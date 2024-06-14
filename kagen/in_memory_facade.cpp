#include "kagen/in_memory_facade.h"

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/factories.h"
#include "kagen/generators/generator.h"
#include "kagen/io.h"
#include "kagen/tools/statistics.h"
#include "kagen/tools/validator.h"

#include <mpi.h>

#include <cmath>
#include <iomanip>
#include <iostream>

namespace kagen {
void GenerateInMemoryToDisk(PGeneratorConfig config, MPI_Comm comm) {
    PEID size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    auto graph = GenerateInMemory(config, GraphRepresentation::EDGE_LIST, comm);

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
}

Graph GenerateInMemory(const PGeneratorConfig& config_template, GraphRepresentation representation, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const bool output_error = rank == ROOT;
    const bool output_info  = rank == ROOT && !config_template.quiet;

    if (output_info && config_template.print_header) {
        PrintHeader(config_template);
    }

    auto             factory = CreateGeneratorFactory(config_template.generator);
    PGeneratorConfig config;
    try {
        config = factory->NormalizeParameters(config_template, rank, size, output_info);
    } catch (ConfigurationError& ex) {
        if (output_error) {
            std::cerr << "Error: " << ex.what() << "\n";
        }
        MPI_Barrier(comm);
        MPI_Abort(comm, 1);
    }

    if (output_info) {
        std::cout << "Generating graph ... " << std::flush;
    }

    const auto t_start_graphgen = MPI_Wtime();

    auto generator = factory->Create(config, rank, size);
    generator->Generate(representation);
    MPI_Barrier(comm);

    if (output_info) {
        std::cout << "OK" << std::endl;
    }

    generator->GenerateEdgeWeights(config.edge_weights, comm);

    const SInt num_edges_before_finalize = generator->GetNumberOfEdges();
    if (output_info) {
        std::cout << "Finalizing graph ... " << std::flush;
    }
    if (!config.skip_postprocessing) {
        generator->Finalize(comm);
        MPI_Barrier(comm);
    }
    if (output_info) {
        std::cout << "OK" << std::endl;
    }
    const SInt num_edges_after_finalize = generator->GetNumberOfEdges();

    const auto t_end_graphgen = MPI_Wtime();

    if (!config.skip_postprocessing && !config.quiet) {
        SInt num_global_edges_before, num_global_edges_after;
        MPI_Reduce(&num_edges_before_finalize, &num_global_edges_before, 1, KAGEN_MPI_SINT, MPI_SUM, ROOT, comm);
        MPI_Reduce(&num_edges_after_finalize, &num_global_edges_after, 1, KAGEN_MPI_SINT, MPI_SUM, ROOT, comm);

        if (num_global_edges_before != num_global_edges_after && output_info) {
            std::cout << "The number of edges changed from " << num_global_edges_before << " to "
                      << num_global_edges_after << " during finalization (= by "
                      << std::abs(
                             static_cast<SSInt>(num_global_edges_after) - static_cast<SSInt>(num_global_edges_before))
                      << ")" << std::endl;
        }
    }

    auto graph = generator->Take();

    // Validation
    if (config.validate_simple_graph) {
        if (output_info) {
            std::cout << "Validating graph ... " << std::flush;
        }

        bool success = ValidateGraph(graph, config.self_loops, config.directed, false, comm);
        MPI_Allreduce(MPI_IN_PLACE, &success, 1, MPI_C_BOOL, MPI_LOR, comm);
        if (!success) {
            if (output_error) {
                std::cerr << "Error: graph validation failed\n";
            }
            MPI_Abort(comm, 1);
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
