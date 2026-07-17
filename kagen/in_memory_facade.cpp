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
#include <sstream>

namespace kagen {
namespace {
// Output formats that require a single PE to own a vertex's whole adjacency (i.e. they group edges by vertex on
// one PE), incompatible with a graph that has split vertices (see the compatibility guard in GenerateInMemory).
bool RequiresSingleVertexOwnership(const FileFormat format) {
    switch (format) {
        case FileFormat::METIS:
        case FileFormat::HMETIS:
        case FileFormat::HMETIS_DIRECTED:
        case FileFormat::HMETIS_EP:
        case FileFormat::DOT:
        case FileFormat::DOT_DIRECTED:
        case FileFormat::PARHIP:
        case FileFormat::XTRAPULP:
        case FileFormat::FREIGHT_NETL:
        case FileFormat::FREIGHT_NETL_EP:
        case FileFormat::NETD_ARE:
            return true;
        default:
            return false;
    }
}
} // namespace

void CheckSplitVertexCompatibility(
    const bool any_split, [[maybe_unused]] const GraphRepresentation representation, const PGeneratorConfig& config) {
    if (!any_split) {
        return;
    }
    // Note: CSR representation is *not* rejected here. A split-vertex CSR graph is only ever produced by the
    // direct strict edge-balanced file read (ParhipReader::ReadStrictEdgeRange), which builds valid xadj/adjncy
    // with partial leading/trailing rows (described by partial_vertices). The paths that would build a
    // *corrupt* CSR from a split edge list -- EdgeListOnlyGenerator::FinalizeCSR and
    // FileGraphGenerator::FinalizeCSR's edge-list branch -- throw before reaching this guard, so any CSR+split
    // graph that gets here is the well-formed kind. A consumer iterating local CSR rows must still honor
    // partial_vertices at the boundaries; the output-format check below rejects adjacency-grouped writers.
    for (const FileFormat& format: config.output_graph.formats) {
        if (RequiresSingleVertexOwnership(format)) {
            std::stringstream msg;
            msg << "the generated graph has vertices whose own edges are split across multiple PEs (from "
                   "--redistribution=balance-edges-strict), which is incompatible with the adjacency-grouped "
                   "output format '"
                << format << "'; use an edge-list-shaped format (edgelist, binary-edgelist) instead";
            throw ConfigurationError(msg.str());
        }
    }
    if (config.validate_simple_graph) {
        throw ConfigurationError(
            "--validate-simple-graph is not supported together with a graph that has split vertices (from "
            "--redistribution=balance-edges-strict)");
    }
    if (config.edge_weights.generator_type == EdgeWeightGeneratorType::HASHING_BASED
        || config.edge_weights.generator_type == EdgeWeightGeneratorType::UNIFORM_RANDOM) {
        std::stringstream msg;
        msg << "edge weight generator '" << config.edge_weights.generator_type
            << "' relies on single-PE vertex ownership and is not supported together with a graph that has "
               "split vertices (from --redistribution=balance-edges-strict)";
        throw ConfigurationError(msg.str());
    }
    if (!config.quiet && config.statistics_level != StatisticsLevel::NONE) {
        // Matches the actual gating in GenerateInMemory ("if (!config.quiet) { ... if (statistics_level >=
        // BASIC) ... }"): with --quiet, no statistics are computed regardless of statistics_level.
        throw ConfigurationError(
            "statistics are not supported together with a graph that has split vertices (from "
            "--redistribution=balance-edges-strict)");
    }
}

void GenerateInMemoryToDisk(PGeneratorConfig config, MPI_Comm comm) {
    PEID size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    auto graph = GenerateInMemory(config, GraphRepresentation::EDGE_LIST, comm);

    const auto t_start_io = MPI_Wtime();

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

    const auto t_end_io = MPI_Wtime();

    if (!config.quiet && rank == ROOT) {
        std::cout << "IO took " << std::fixed << std::setprecision(3) << t_end_io - t_start_io << " seconds"
                  << std::endl;
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
    } catch (const kagen::ConfigurationError& ex) {
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

    // Any of the following stages may throw a ConfigurationError for a runtime-data-dependent condition (e.g. a
    // split vertex combined with an incompatible output format) that could not be caught by NormalizeParameters
    // above; handle it the same way here so it fails fast with a clear message instead of an uncaught abort.
    std::unique_ptr<Generator> generator;
    Graph                      graph;
    double                     t_end_graphgen = t_start_graphgen;
    try {
        generator = factory->Create(config, rank, size);
        generator->Generate(representation);
        MPI_Barrier(comm);

        if (output_info) {
            std::cout << "OK" << std::endl;
        }

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

        // Compatibility guard: BALANCE_EDGES_TRUE may split a vertex's own edges across multiple PEs, breaking
        // the invariant that a single PE owns a vertex's whole adjacency. Several consumers below (and further
        // downstream, in WriteGraph) rely on that invariant and would silently misbehave -- not just produce a
        // suboptimal result -- if it doesn't hold, so fail fast here instead. This is the one place all CLI and
        // library entry points funnel through, and the only point where has_split_vertices -- a runtime,
        // data-dependent property -- is known.
        bool any_split = generator->HasSplitVertices();
        MPI_Allreduce(MPI_IN_PLACE, &any_split, 1, MPI_C_BOOL, MPI_LOR, comm);
        CheckSplitVertexCompatibility(any_split, representation, config);

        if (output_info) {
            std::cout << "Generating weights ... " << std::flush;
        }
        generator->GenerateEdgeWeights(config.edge_weights, comm);
        generator->GenerateVertexWeights(config.vertex_weights, comm);
        if (output_info) {
            std::cout << "OK" << std::endl;
        }

        t_end_graphgen = MPI_Wtime();

        if (!config.skip_postprocessing && !config.quiet) {
            SInt num_global_edges_before, num_global_edges_after;
            MPI_Reduce(&num_edges_before_finalize, &num_global_edges_before, 1, KAGEN_MPI_SINT, MPI_SUM, ROOT, comm);
            MPI_Reduce(&num_edges_after_finalize, &num_global_edges_after, 1, KAGEN_MPI_SINT, MPI_SUM, ROOT, comm);

            if (num_global_edges_before != num_global_edges_after && output_info) {
                std::cout << "The number of edges changed from " << num_global_edges_before << " to "
                          << num_global_edges_after << " during finalization (= by "
                          << std::abs(
                                 static_cast<SSInt>(num_global_edges_after)
                                 - static_cast<SSInt>(num_global_edges_before))
                          << ")" << std::endl;
            }
        }
        if (config.permute) {
            generator->PermuteVertices(config, comm);
        }

        graph = generator->Take();
    } catch (const kagen::ConfigurationError& ex) {
        if (output_error) {
            std::cerr << "Error: " << ex.what() << "\n";
        }
        MPI_Barrier(comm);
        MPI_Abort(comm, 1);
    }

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
                PrintAdvancedStatistics(graph.xadj, graph.adjncy, graph.vertex_range, rank == ROOT, comm);
            }
        }
        if (output_info && config.statistics_level != StatisticsLevel::NONE) {
            std::cout << "-------------------------------------------------------------------------------" << std::endl;
        }
    }

    return graph;
}
} // namespace kagen
