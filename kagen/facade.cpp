#include "kagen/facade.h"

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/generators/generator.h"

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

std::pair<EdgeList, VertexRange> Generate(const PGeneratorConfig& config_template) {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return Generate(config_template, rank, size);
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
        return IsCubic(value * 2) ? value * 2 : value * 3;
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

std::pair<EdgeList, VertexRange> Generate(const PGeneratorConfig& config_template, const PEID rank, const PEID size) {
    const bool output = !config_template.quiet && rank == ROOT;
    auto       config = config_template;

    // Get generator requirements @todo this is ugly
    config.k               = size;
    const int requirements = CreateGenerator(config, rank, size)->Requirements();
    config.k               = config_template.k;

    const bool require_power_of_two_communicator = requirements & GeneratorRequirement::POWER_OF_TWO_COMMUNICATOR_SIZE;
    const bool require_square_chunks             = requirements & GeneratorRequirement::SQAURE_CHUNKS;
    const bool require_cubic_chunks              = requirements & GeneratorRequirement::CUBIC_CHUNKS;

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

    // Generate graph
    auto generator             = CreateGenerator(config, rank, size);
    auto [edges, vertex_range] = generator->Generate();

    // Postprocessing
    if (generator->AlmostUndirected()) {
        AddReverseEdges(edges, vertex_range);
    }

    // Validation
    if (config.validate_simple_graph) {
        bool success = ValidateSimpleGraph(edges, vertex_range);
        MPI_Allreduce(MPI_IN_PLACE, &success, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
        if (!success) {
            std::exit(1);
        }
    }

    return {std::move(edges), vertex_range};
}
} // namespace kagen
