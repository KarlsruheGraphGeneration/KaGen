#include "kagen/facade.h"

#include <cmath>
#include <iomanip>
#include <iostream>

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/tools/statistics.h"

// Generators
#include "kagen/generators/barabassi/barabassi.h"
#include "kagen/generators/geometric/rgg.h"
#include "kagen/generators/gnm/gnm_directed.h"
#include "kagen/generators/gnm/gnm_undirected.h"
#include "kagen/generators/gnp/gnp_directed.h"
#include "kagen/generators/gnp/gnp_undirected.h"
#include "kagen/generators/grid/grid_2d.h"
#include "kagen/generators/grid/grid_3d.h"
#include "kagen/generators/hyperbolic/hyperbolic.h"
#include "kagen/generators/kronecker/kronecker.h"

#ifdef KAGEN_CGAL_FOUND
    #include "kagen/generators/geometric/delaunay.h"
#endif // KAGEN_CGAL_FOUND

#include "kagen/tools/postprocessor.h"
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
#endif // KAGEN_CGAL_FOUND

        case GeneratorType::GRID_2D:
            return std::make_unique<Grid2DFactory>();

        case GeneratorType::GRID_3D:
            return std::make_unique<Grid3DFactory>();

        case GeneratorType::BA:
            return std::make_unique<BarabassiFactory>();

        case GeneratorType::KRONECKER:
            return std::make_unique<KroneckerFactory>();

        case GeneratorType::RHG:
            return std::make_unique<HyperbolicFactory>();
    }

    __builtin_unreachable();
}

namespace {
bool IsPowerOfTwo(const SInt value) {
    return (value & (value - 1)) == 0;
}

bool IsSquare(const SInt value) {
    const SInt root = std::round(std::sqrt(value));
    return root * root == value;
}

bool IsCubic(const SInt value) {
    const SInt root = std::round(std::cbrt(value));
    return root * root * root == value;
}

SInt FindMultipleSquare(const SInt value) {
    if (IsSquare(value)) {
        return value;
    }
    if (IsPowerOfTwo(value)) {
        return 2 * value; // every 2nd power of two is square
    }

    // find smallest square number that is a multiple of value
    const SInt root = std::sqrt(value);
    for (SInt cur = root; cur < value; ++cur) {
        const SInt squared = cur * cur;
        if (squared % value == 0) {
            return squared;
        }
    }
    return value * value;
}

SInt FindMultipleCube(const SInt value) {
    if (IsCubic(value)) {
        return value;
    }
    if (IsPowerOfTwo(value)) {
        return IsCubic(value * 2) ? value * 2 : value * 4;
    }

    // find smallest cubic number that is a multiple of value
    const SInt root = std::cbrt(value);
    for (SInt cur = root; cur < value; ++cur) {
        const SInt cubed = cur * cur * cur;
        if (cubed % value == 0) {
            return cubed;
        }
    }
    return value * value * value;
}
} // namespace

std::tuple<EdgeList, VertexRange, Coordinates> Generate(const PGeneratorConfig& config_template, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const bool output_error = rank == ROOT;
    const bool output_info  = rank == ROOT && !config_template.quiet;

    auto factory = CreateGeneratorFactory(config_template.generator);
    PGeneratorConfig config;
    try {
        config = factory->NormalizeParameters(config_template, output_info);
    } catch (ConfigurationError& ex) {
        if (output_error) {
            std::cerr << "Error: " << ex.what() << "\n";
        }
        std::exit(1);
    }

    // Get generator requirements
    const int  requirements                      = factory->Requirements();
    const bool require_power_of_two_communicator = requirements & GeneratorRequirement::POWER_OF_TWO_COMMUNICATOR_SIZE;
    const bool require_square_chunks             = requirements & GeneratorRequirement::SQUARE_CHUNKS;
    const bool require_cubic_chunks              = requirements & GeneratorRequirement::CUBIC_CHUNKS;
    const bool require_one_chunk_per_pe          = requirements & GeneratorRequirement::ONE_CHUNK_PER_PE;

    // Validate number of PEs
    if (require_power_of_two_communicator && !IsPowerOfTwo(size)) {
        if (output_error) {
            std::cerr << "Error: generator requires the number of PEs to be a power of two\n";
        }
        std::exit(1);
    }

    // Setup chunk count
    if (config.k == 0) { // automatic selection
        if (require_square_chunks) {
            config.k = FindMultipleSquare(size);
        } else if (require_cubic_chunks) {
            config.k = FindMultipleCube(size);
        } else {
            config.k = size;
        }

        if (output_info) {
            std::cout << "Setting number of chunks to " << config.k << std::endl;
        }
    } else { // only validation
        if (require_square_chunks && !IsSquare(config.k)) {
            if (output_error) {
                std::cerr << "Error: generator requires a square number of chunks\n";
            }
            std::exit(1);
        }
        if (require_cubic_chunks && !IsCubic(config.k)) {
            if (output_error) {
                std::cerr << "Error: generator requires a cubic number of chunks\n";
            }
            std::exit(1);
        }
    }
    if (require_one_chunk_per_pe && config.k != static_cast<SInt>(size)) {
        if (output_error) {
            std::cerr << "Error: generator requires exactly one chunk per PE\n";
        }
        std::exit(1);
    }

    // Generate graph
    if (output_info) {
        std::cout << "Generating graph ..." << std::endl;
    }

    const auto start_graphgen = MPI_Wtime();

    auto generator = factory->Create(config, rank, size);
    generator->Generate();
    auto edges        = generator->TakeEdges();
    auto vertex_range = generator->GetVertexRange();
    auto coordinates  = generator->GetCoordinates();

    // Postprocessing: cut edges
    const SInt num_edges_before_postprocessing = edges.size();
    if (!config.skip_postprocessing && generator->IsAlmostUndirected()) {
        AddReverseEdges(edges, vertex_range, comm);
    }
    const SInt num_edges_after_postprocessing = edges.size();

    const auto end_graphgen = MPI_Wtime();

    // Postprocessing
    if (!config.skip_postprocessing && generator->IsAlmostUndirected()) {
        if (output_info) {
            std::cout << "Postprocessing:" << std::endl;
        }

        SInt num_global_edges_before, num_global_edges_after;
        MPI_Reduce(&num_edges_before_postprocessing, &num_global_edges_before, 1, KAGEN_MPI_SINT, MPI_SUM, ROOT, comm);
        MPI_Reduce(&num_edges_after_postprocessing, &num_global_edges_after, 1, KAGEN_MPI_SINT, MPI_SUM, ROOT, comm);

        if (output_info) {
            std::cout << "  Number of global edges changed from " << num_global_edges_before << " to "
                      << num_global_edges_after << " edges: by "
                      << std::abs(
                             static_cast<SSInt>(num_global_edges_after) - static_cast<SSInt>(num_global_edges_before))
                      << std::endl;
        }
    }

    // Validation
    if (config.validate_simple_graph) {
        if (output_info) {
            std::cout << "Validating simple graph ... " << std::flush;
        }

        bool success = ValidateSimpleGraph(edges, vertex_range, comm);
        MPI_Allreduce(MPI_IN_PLACE, &success, 1, MPI_C_BOOL, MPI_LOR, comm);
        if (!success) {
            if (output_error) {
                std::cerr << "Error: simple graph validation failed\n";
            }
            std::exit(1);
        } else if (output_info) {
            std::cout << "ok" << std::endl;
        }
    }

    // Statistics
    if (!config.quiet) {
        // Running time
        if (output_info) {
            std::cout << "Generation took " << std::fixed << std::setprecision(3) << end_graphgen - start_graphgen
                      << " seconds" << std::endl;
        }

        if (config.statistics_level >= StatisticsLevel::BASIC) {
            // Basic graph statistics
            PrintBasicStatistics(edges, vertex_range, rank == ROOT, comm);
        }
        if (config.statistics_level >= StatisticsLevel::ADVANCED) {
            PrintAdvancedStatistics(edges, vertex_range, rank == ROOT, comm);
        }
    }

    return {std::move(edges), vertex_range, std::move(coordinates)};
}
} // namespace kagen
