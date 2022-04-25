#include "kagen/facade.h"

#include <iomanip>
#include <iostream>

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/tools/statistics.h"

// Generators
#include "kagen/generators/barabassi/barabassi.h"
#include "kagen/generators/geometric/rgg/rgg_2d.h"
#include "kagen/generators/geometric/rgg/rgg_3d.h"
#include "kagen/generators/gnm/gnm_directed.h"
#include "kagen/generators/gnm/gnm_undirected.h"
#include "kagen/generators/gnp/gnp_directed.h"
#include "kagen/generators/gnp/gnp_undirected.h"
#include "kagen/generators/grid/grid_2d.h"
#include "kagen/generators/grid/grid_3d.h"
#include "kagen/generators/hyperbolic/hyperbolic.h"
#include "kagen/generators/kronecker/kronecker.h"

#ifdef KAGEN_CGAL_FOUND
    #include "kagen/generators/geometric/delaunay/delaunay_2d.h"
    #include "kagen/generators/geometric/delaunay/delaunay_3d.h"
#endif // KAGEN_CGAL_FOUND

#include "kagen/tools/postprocessor.h"
#include "kagen/tools/validator.h"

namespace kagen {
std::unique_ptr<Generator> CreateGenerator(const PGeneratorConfig& config, const PEID rank, const PEID size) {
    switch (config.generator) {
        case GeneratorType::GNM_DIRECTED:
            return std::make_unique<GNMDirected>(config, rank, size);

        case GeneratorType::GNM_UNDIRECTED:
            return std::make_unique<GNMUndirected>(config, rank, size);

        case GeneratorType::GNP_DIRECTED:
            return std::make_unique<GNPDirected>(config, rank, size);

        case GeneratorType::GNP_UNDIRECTED:
            return std::make_unique<GNPUndirected>(config, rank, size);

        case GeneratorType::RGG_2D:
            return std::make_unique<RGG2D>(config, rank, size);

        case GeneratorType::RGG_3D:
            return std::make_unique<RGG3D>(config, rank, size);

#ifdef KAGEN_CGAL_FOUND
        case GeneratorType::RDG_2D:
            return std::make_unique<Delaunay2D>(config, rank, size);

        case GeneratorType::RDG_3D:
            return std::make_unique<Delaunay3D>(config, rank, size);
#endif // KAGEN_CGAL_FOUND

        case GeneratorType::GRID_2D:
            return std::make_unique<Grid2D>(config, rank, size);

        case GeneratorType::GRID_3D:
            return std::make_unique<Grid3D>(config, rank, size);

        case GeneratorType::BA:
            return std::make_unique<Barabassi>(config, rank, size);

        case GeneratorType::KRONECKER:
            return std::make_unique<Kronecker>(config, rank, size);

        case GeneratorType::RHG:
            return std::make_unique<Hyperbolic>(config, rank, size);
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

std::pair<EdgeList, VertexRange> Generate(const PGeneratorConfig& config_template, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const bool output = rank == ROOT; // && !config_template.quiet; // Always output errors that crash the program
    auto       config = config_template;

    // Get generator requirements @todo this is ugly
    config.k               = size;
    const int requirements = CreateGenerator(config, rank, size)->Requirements();
    config.k               = config_template.k;

    const bool require_power_of_two_communicator = requirements & GeneratorRequirement::POWER_OF_TWO_COMMUNICATOR_SIZE;
    const bool require_square_chunks             = requirements & GeneratorRequirement::SQUARE_CHUNKS;
    const bool require_cubic_chunks              = requirements & GeneratorRequirement::CUBIC_CHUNKS;
    const bool require_one_chunk_per_pe          = requirements & GeneratorRequirement::ONE_CHUNK_PER_PE;

    // Validate number of PEs
    if (require_power_of_two_communicator && !IsPowerOfTwo(size)) {
        if (output) {
            std::cerr << "Generator requires the number of PEs to be a power of two\n";
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

        if (output) {
            std::cout << "Set number of chunks to " << config.k << std::endl;
        }
    } else { // only validation
        if (require_square_chunks && !IsSquare(config.k)) {
            if (output) {
                std::cerr << "Generator requires a square number of chunks\n";
            }
            std::exit(1);
        }
        if (require_cubic_chunks && !IsCubic(config.k)) {
            if (output) {
                std::cerr << "Generator requires a cubic number of chunks\n";
            }
            std::exit(1);
        }
    }
    if (require_one_chunk_per_pe && config.k != static_cast<SInt>(size)) {
        if (output) {
            std::cerr << "Generator requires exactly one chunk per PE\n";
        }
        std::exit(1);
    }

    // Generate graph
    const auto start_graphgen  = MPI_Wtime();
    auto       generator       = CreateGenerator(config, rank, size);
    auto [edges, vertex_range] = generator->Generate();

    // Postprocessing
    if (generator->AlmostUndirected()) {
        AddReverseEdges(edges, vertex_range, comm);
    }
    const auto end_graphgen = MPI_Wtime();

    // Validation
    if (config.validate_simple_graph) {
        bool success = ValidateSimpleGraph(edges, vertex_range, comm);
        MPI_Allreduce(MPI_IN_PLACE, &success, 1, MPI_C_BOOL, MPI_LOR, comm);
        if (!success) {
            if (output) {
                std::cerr << "Simple graph validation failed\n";
            }
            std::exit(1);
        }
    }

    // Statistics
    if (!config.quiet) {
        if (config.statistics_level >= StatisticsLevel::BASIC) {
            if (output) {
                std::cout << "Generation took " << std::fixed << std::setprecision(3) << end_graphgen - start_graphgen
                          << " seconds" << std::endl;
            }
            PrintBasicStatistics(edges, vertex_range, rank == ROOT, comm);
        }
        if (config.statistics_level >= StatisticsLevel::ADVANCED) {
            PrintAdvancedStatistics(edges, vertex_range, rank == ROOT, comm);
        }
    }

    return {std::move(edges), vertex_range};
}
} // namespace kagen
