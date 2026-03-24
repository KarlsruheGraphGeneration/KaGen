#include "kagen/in_memory_facade.h"

#include "kagen/comm/comm.h"
#include "kagen/comm/comm_types.h"
#include "kagen/comm/hybrid_comm.h"
#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/factories.h"
#include "kagen/generators/generator.h"
#include "kagen/io.h"
#include "kagen/tools/statistics.h"
#include "kagen/tools/validator.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

namespace kagen {
void GenerateInMemoryToDisk(PGeneratorConfig config, Comm& comm) {
    PEID size = comm.Size();
    PEID rank = comm.Rank();

    auto graph = GenerateInMemory(config, GraphRepresentation::EDGE_LIST, comm);

    const auto t_start_io = comm.Wtime();

    const std::string base_filename = config.output_graph.filename;
    for (const FileFormat& format: config.output_graph.formats) {
        const auto& factory = GetGraphFormatFactory(format);

        const std::string filename   = (config.output_graph.extension && !factory->DefaultExtensions().empty())
                                           ? base_filename + "." + factory->DefaultExtensions().front()
                                           : base_filename;
        config.output_graph.filename = filename;

        GraphInfo info(graph, comm);
        auto      writer = factory->CreateWriter(config.output_graph, graph, info, rank, size);
        if (writer != nullptr) {
            WriteGraph(*writer.get(), config.output_graph, rank == ROOT && !config.quiet, comm);
        } else if (!config.quiet && rank == ROOT) {
            std::cout << "Warning: invalid file format " << format << " for writing; skipping\n";
        }
    }

    const auto t_end_io = comm.Wtime();

    if (!config.quiet && rank == ROOT) {
        std::cout << "IO took " << std::fixed << std::setprecision(3) << t_end_io - t_start_io << " seconds"
                  << std::endl;
    }
}

Graph GenerateInMemory(const PGeneratorConfig& config_template, GraphRepresentation representation, Comm& comm) {
    // Threaded path: spawn num_threads virtual PEs within this real PE
    if (config_template.num_threads > 1) {
        const int  T            = config_template.num_threads;
        const PEID real_rank = comm.Rank();

        // Per-thread config: suppress output/stats/validation (done once on merged result)
        PGeneratorConfig thread_config    = config_template;
        thread_config.quiet               = true;
        thread_config.statistics_level    = StatisticsLevel::NONE;
        thread_config.validate_simple_graph = false;
        thread_config.num_threads         = 1; // Prevent recursive threading

        HybridCommShared         shared(T);
        std::vector<Graph>       thread_graphs(T);
        std::vector<std::thread> workers;
        workers.reserve(T);
        for (int t = 0; t < T; ++t) {
            workers.emplace_back([&, t]() {
                HybridComm hcomm(t, comm, shared);
                thread_graphs[t] = GenerateInMemory(thread_config, representation, hcomm);
            });
        }
        for (auto& w : workers) {
            w.join();
        }

        Graph merged = std::move(thread_graphs[0]);
        for (int t = 1; t < T; ++t) {
            merged.Merge(std::move(thread_graphs[t]));
        }

        // Run post-generation steps (validation, statistics) on the merged graph using the real comm
        if (config_template.validate_simple_graph) {
            const bool output_error = real_rank == ROOT;
            const bool output_info  = real_rank == ROOT && !config_template.quiet;
            if (output_info) {
                std::cout << "Validating graph ... " << std::flush;
            }
            bool success = ValidateGraph(merged, config_template.self_loops, config_template.directed, false, comm);
            comm.Allreduce(COMM_IN_PLACE, &success, 1, CommDatatype::C_BOOL, CommOp::LOR);
            if (!success) {
                if (output_error) {
                    std::cerr << "Error: graph validation failed\n";
                }
                comm.Abort(1);
            } else if (output_info) {
                std::cout << "OK" << std::endl;
            }
        }
        if (!config_template.quiet) {
            if (representation == GraphRepresentation::EDGE_LIST) {
                if (config_template.statistics_level >= StatisticsLevel::BASIC) {
                    PrintBasicStatistics(merged.edges, merged.vertex_range, real_rank == ROOT, comm);
                }
                if (config_template.statistics_level >= StatisticsLevel::ADVANCED) {
                    PrintAdvancedStatistics(merged.edges, merged.vertex_range, real_rank == ROOT, comm);
                }
            } else {
                if (config_template.statistics_level >= StatisticsLevel::BASIC) {
                    PrintBasicStatistics(merged.xadj, merged.adjncy, merged.vertex_range, real_rank == ROOT, comm);
                }
            }
        }

        return merged;
    }

    PEID rank = comm.Rank();
    PEID size = comm.Size();

    const bool output_error = rank == ROOT;
    const bool output_info  = rank == ROOT && !config_template.quiet;

    if (output_info && config_template.print_header) {
        PrintHeader(config_template);
    }

    auto             factory = CreateGeneratorFactory(config_template.generator);
    PGeneratorConfig config;
    try {
        config = factory->NormalizeParameters(config_template, rank, size, output_info);
    } catch (const kagen::ConfigurationError& ex) {
        if (output_error) {
            std::cerr << "Error: " << ex.what() << "\n";
        }
        comm.Barrier();
        comm.Abort(1);
    }

    if (output_info) {
        std::cout << "Generating graph ... " << std::flush;
    }

    const auto t_start_graphgen = comm.Wtime();

    auto generator = factory->Create(config, rank, size);
    generator->Generate(representation);
    comm.Barrier();

    if (output_info) {
        std::cout << "OK" << std::endl;
    }

    const SInt num_edges_before_finalize = generator->GetNumberOfEdges();
    if (output_info) {
        std::cout << "Finalizing graph ... " << std::flush;
    }
    if (!config.skip_postprocessing) {
        generator->Finalize(comm);
        comm.Barrier();
    }
    if (output_info) {
        std::cout << "OK" << std::endl;
    }
    const SInt num_edges_after_finalize = generator->GetNumberOfEdges();

    if (output_info) {
        std::cout << "Generating weights ... " << std::flush;
    }
    generator->GenerateEdgeWeights(config.edge_weights, comm);
    generator->GenerateVertexWeights(config.vertex_weights, comm);
    if (output_info) {
        std::cout << "OK" << std::endl;
    }

    const auto t_end_graphgen = comm.Wtime();

    if (!config.skip_postprocessing && !config.quiet) {
        SInt num_global_edges_before = 0, num_global_edges_after = 0;
        comm.Reduce(
            &num_edges_before_finalize, &num_global_edges_before, 1, CommDatatype::UNSIGNED_LONG_LONG, CommOp::SUM,
            ROOT);
        comm.Reduce(
            &num_edges_after_finalize, &num_global_edges_after, 1, CommDatatype::UNSIGNED_LONG_LONG, CommOp::SUM,
            ROOT);

        if (num_global_edges_before != num_global_edges_after && output_info) {
            std::cout << "The number of edges changed from " << num_global_edges_before << " to "
                      << num_global_edges_after << " during finalization (= by "
                      << std::abs(
                             static_cast<SSInt>(num_global_edges_after) - static_cast<SSInt>(num_global_edges_before))
                      << ")" << std::endl;
        }
    }
    if (config.permute) {
        generator->PermuteVertices(config, comm);
    }

    auto graph = generator->Take();

    // Validation
    if (config.validate_simple_graph) {
        if (output_info) {
            std::cout << "Validating graph ... " << std::flush;
        }

        bool success = ValidateGraph(graph, config.self_loops, config.directed, false, comm);
        comm.Allreduce(COMM_IN_PLACE, &success, 1, CommDatatype::C_BOOL, CommOp::LOR);
        if (!success) {
            if (output_error) {
                std::cerr << "Error: graph validation failed\n";
            }
            comm.Abort(1);
        } else if (output_info) {
            std::cout << "OK" << std::endl;
        }
    }

    // Statistics
    if (!config.quiet) {
        if (output_info) {
            std::cout << "Generation took " << std::fixed << std::setprecision(3) << t_end_graphgen - t_start_graphgen
                      << " seconds" << std::endl;
            std::cout << "-------------------------------------------------------------------------------" << std::endl;
        }

        if (representation == GraphRepresentation::EDGE_LIST) {
            if (config.statistics_level >= StatisticsLevel::BASIC) {
                PrintBasicStatistics(graph.edges, graph.vertex_range, rank == ROOT, comm);
            }
            if (config.statistics_level >= StatisticsLevel::ADVANCED) {
                PrintAdvancedStatistics(graph.edges, graph.vertex_range, rank == ROOT, comm);
            }
        } else { // CSR
            if (config.statistics_level >= StatisticsLevel::BASIC) {
                PrintBasicStatistics(graph.xadj, graph.adjncy, graph.vertex_range, rank == ROOT, comm);
            }
            if (config.statistics_level >= StatisticsLevel::ADVANCED) {
                std::cout << "Advanced statistics are not available when generating the graph in CSR representation"
                          << std::endl;
            }
        }
        if (output_info && config.statistics_level != StatisticsLevel::NONE) {
            std::cout << "-------------------------------------------------------------------------------" << std::endl;
        }
    }

    return graph;
}
} // namespace kagen
