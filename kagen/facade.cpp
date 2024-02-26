#include "kagen/facade.h"

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/tools/statistics.h"

#include <mpi.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

// Generators
#include "kagen/generators/barabassi/barabassi.h"
#include "kagen/generators/file/file_graph.h"
#include "kagen/generators/geometric/rgg.h"
#include "kagen/generators/gnm/gnm_directed.h"
#include "kagen/generators/gnm/gnm_undirected.h"
#include "kagen/generators/gnp/gnp_directed.h"
#include "kagen/generators/gnp/gnp_undirected.h"
#include "kagen/generators/grid/grid_2d.h"
#include "kagen/generators/grid/grid_3d.h"
#include "kagen/generators/hyperbolic/hyperbolic.h"
#include "kagen/generators/image/image_mesh.h"
#include "kagen/generators/kronecker/kronecker.h"
#include "kagen/generators/path/path_directed.h"
#include "kagen/generators/rmat/rmat.h"

#ifdef KAGEN_CGAL_FOUND
    #include "kagen/generators/geometric/delaunay.h"
#endif // KAGEN_CGAL_FOUND

#include "kagen/tools/validator.h"

namespace kagen {
std::unique_ptr<GeneratorFactory> CreateGeneratorFactory(const GeneratorType type) {
    switch (type) {
        case GeneratorType::GNM_DIRECTED:
            return std::make_unique<GNMDirectedFactory>();

        case GeneratorType::GNM_UNDIRECTED:
            return std::make_unique<GNMUndirectedFactory>();

        case GeneratorType::GNP_DIRECTED:
            return std::make_unique<GNPDirectedFactory>();

        case GeneratorType::GNP_UNDIRECTED:
            return std::make_unique<GNPUndirectedFactory>();

        case GeneratorType::RGG_2D:
            return std::make_unique<RGG2DFactory>();

        case GeneratorType::RGG_3D:
            return std::make_unique<RGG3DFactory>();

#ifdef KAGEN_CGAL_FOUND
        case GeneratorType::RDG_2D:
            return std::make_unique<Delaunay2DFactory>();

        case GeneratorType::RDG_3D:
            return std::make_unique<Delaunay3DFactory>();
#else  // KAGEN_CGAL_FOUND
        case GeneratorType::RDG_2D:
        case GeneratorType::RDG_3D:
            // throw exception after switch
            break;
#endif // KAGEN_CGAL_FOUND

        case GeneratorType::GRID_2D:
            return std::make_unique<Grid2DFactory>();

        case GeneratorType::GRID_3D:
            return std::make_unique<Grid3DFactory>();

        case GeneratorType::PATH_DIRECTED:
            return std::make_unique<PathDirectedFactory>();

        case GeneratorType::BA:
            return std::make_unique<BarabassiFactory>();

        case GeneratorType::KRONECKER:
            return std::make_unique<KroneckerFactory>();

        case GeneratorType::RHG:
            return std::make_unique<HyperbolicFactory>();

        case GeneratorType::RMAT:
            return std::make_unique<RMATFactory>();

        case GeneratorType::IMAGE_MESH:
            return std::make_unique<ImageMeshFactory>();

        case GeneratorType::FILE:
            return std::make_unique<FileGraphFactory>();
    }

    throw std::runtime_error("invalid graph generator type");
}

namespace {
void PrintHeader(const PGeneratorConfig& config) {
    std::cout << "###############################################################################\n";
    std::cout << "#                         _  __      ____                                     #\n";
    std::cout << "#                        | |/ /__ _ / ___| ___ _ __                           #\n";
    std::cout << "#                        | ' // _` | |  _ / _ \\ '_ \\                          #\n";
    std::cout << "#                        | . \\ (_| | |_| |  __/ | | |                         #\n";
    std::cout << "#                        |_|\\_\\__,_|\\____|\\___|_| |_|                         #\n";
    std::cout << "#                         Karlsruhe Graph Generation                          #\n";
    std::cout << "#                                                                             #\n";
    std::cout << "###############################################################################\n";
    std::cout << config;
}
} // namespace

Graph Generate(const PGeneratorConfig& config_template, GraphRepresentation representation, MPI_Comm comm) {
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

    // Generate graph
    if (output_info) {
        std::cout << "Generating graph ... " << std::flush;
    }

    const auto start_graphgen = MPI_Wtime();

    auto generator = factory->Create(config, rank, size);
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

    const auto end_graphgen = MPI_Wtime();

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
            std::cout << "Generation took " << std::fixed << std::setprecision(3) << end_graphgen - start_graphgen
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
