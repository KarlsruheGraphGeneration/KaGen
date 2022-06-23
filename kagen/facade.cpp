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
template <typename ApproxRadius, typename ApproxNumNodes>
void SetupRGGParameters(
    PGeneratorConfig& config, ApproxRadius&& approx_radius, ApproxNumNodes&& approx_num_nodes, const bool output_error,
    const bool output_info) {
    // We have three parameters, from which two must be given:
    // - Number of nodes
    // - Number of edges
    // - Radius
    // If number of nodes or radius is missing, compute it from the other two values
    if (config.r == 0.0) {
        if (config.n == 0 || config.m == 0) {
            if (output_error) {
                std::cerr << "At least two parameters out of {n, m, r} must be nonzero\n";
            }
            std::exit(1);
        }

        config.r = approx_radius(config.n, config.m);
        if (output_info) {
            std::cout << "Setting radius to " << config.r << std::endl;
        }
    } else if (config.n == 0) {
        if (config.r == 0.0 || config.m == 0) {
            if (output_error) {
                std::cerr << "At least two parameters out of {n, m, r} must be nonzero\n";
            }
            std::exit(1);
        }

        config.n = approx_num_nodes(config.m, config.r);
        if (output_info) {
            std::cout << "Setting number of nodes to " << config.n << std::endl;
        }
    }
}

void SetupRHGParameters(PGeneratorConfig& config, const bool output_error, const bool output_info) {
    if (config.avg_degree == 0) {
        if (config.m == 0 || config.n == 0) {
            if (output_error) {
                std::cerr << "At least two parameters out of {n, m, d} must be nonzero\n";
            }
            std::exit(1);
        }

        config.avg_degree = 1.0 * config.m / config.n;
        if (output_info) {
            std::cout << "Setting average degree to " << config.avg_degree << std::endl;
        }
    } else if (config.n == 0) {
        if (config.avg_degree == 0 || config.m == 0) {
            if (output_error) {
                std::cerr << "At least two parameters out of {n, m, d} must be nonzero\n";
            }
            std::exit(1);
        }

        config.n = static_cast<SInt>(config.m / config.avg_degree);
        if (output_info) {
            std::cout << "Setting number of nodes to " << config.n << std::endl;
        }
    }
}

#ifdef KAGEN_CGAL_FOUND
void SetupRDGParameters(
    PGeneratorConfig& config, const double factor, const bool output_error, const bool output_info) {
    if (config.n == 0) {
        if (config.m == 0) {
            if (output_error) {
                std::cerr << "At least one parameter out of {n, m} must be nonzero\n";
            }
            std::exit(1);
        }

        config.n = config.m / factor + 2;
        if (output_info) {
            std::cout << "Setting number of nodes to " << config.n << std::endl;
        }
    }
}
#endif // KAGEN_CGAL_FOUND

void ApproxMissingParameters(PGeneratorConfig& config, const double output_error, const double output_info) {
    switch (config.generator) {
        case GeneratorType::RGG_2D:
            SetupRGGParameters(config, &RGG2D::ApproxRadius, &RGG2D::ApproxNumNodes, output_error, output_info);
            break;

        case GeneratorType::RGG_3D:
            SetupRGGParameters(config, &RGG3D::ApproxRadius, &RGG3D::ApproxNumNodes, output_error, output_info);
            break;

        case GeneratorType::RHG:
            SetupRHGParameters(config, output_error, output_info);
            break;

#ifdef KAGEN_CGAL_FOUND
        case GeneratorType::RDG_2D:
            SetupRDGParameters(config, 3, output_error, output_info);
            break;

        case GeneratorType::RDG_3D:
            SetupRDGParameters(config, 15.53, output_error, output_info);
            break;
#endif // KAGEN_CGAL_FOUND

        default:
            // do nothing
            break;
    }
}

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
    auto       config       = config_template;

    // Some parameters can option or can be deduced from other parameters; do that now, before
    // instantiating the generators to obtain their chunk / PE requirements
    ApproxMissingParameters(config, output_error, output_info);

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
        if (output_error) {
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

        if (output_info) {
            std::cout << "Setting number of chunks to " << config.k << std::endl;
        }
    } else { // only validation
        if (require_square_chunks && !IsSquare(config.k)) {
            if (output_error) {
                std::cerr << "Generator requires a square number of chunks\n";
            }
            std::exit(1);
        }
        if (require_cubic_chunks && !IsCubic(config.k)) {
            if (output_error) {
                std::cerr << "Generator requires a cubic number of chunks\n";
            }
            std::exit(1);
        }
    }
    if (require_one_chunk_per_pe && config.k != static_cast<SInt>(size)) {
        if (output_error) {
            std::cerr << "Generator requires exactly one chunk per PE\n";
        }
        std::exit(1);
    }

    // Generate graph
    if (output_info) {
        std::cout << "Generating graph ..." << std::endl;
    }

    const auto start_graphgen = MPI_Wtime();

    auto generator = CreateGenerator(config, rank, size);
    generator->Generate();

    // Postprocessing: vertex ranges
    const VertexRange old_vertex_range = generator->GetVertexRange();
    if (generator->InvalidVertexRangeIfEmpty()) {
        generator->SetVertexRange(RecomputeVertexRanges(generator->GetVertexRange(), comm));
    }

    auto edges        = generator->TakeEdges();
    auto vertex_range = generator->GetVertexRange();
    auto coordinates  = generator->GetCoordinates();

    // Postprocessing: cut edges
    const SInt num_edges_before_postprocessing = edges.size();
    if (!config.skip_postprocessing && generator->AlmostUndirected()) {
        AddReverseEdges(edges, vertex_range, comm);
    }
    const SInt num_edges_after_postprocessing = edges.size();

    const auto end_graphgen = MPI_Wtime();

    // Postprocessing
    if (!config.skip_postprocessing && output_error
        && (generator->AlmostUndirected() || generator->InvalidVertexRangeIfEmpty())) {
        std::cout << "Postprocessing:" << std::endl;
    }
    if (!config.skip_postprocessing && generator->AlmostUndirected()) {
        SInt num_global_edges_before, num_global_edges_after;
        MPI_Reduce(&num_edges_before_postprocessing, &num_global_edges_before, 1, KAGEN_MPI_SINT, MPI_SUM, ROOT, comm);
        MPI_Reduce(&num_edges_after_postprocessing, &num_global_edges_after, 1, KAGEN_MPI_SINT, MPI_SUM, ROOT, comm);

        if (output_error) {
            std::cout << "  Number of global edges changed from " << num_global_edges_before << " to "
                      << num_global_edges_after << " edges: by "
                      << std::abs(
                             static_cast<SSInt>(num_global_edges_after) - static_cast<SSInt>(num_global_edges_before))
                      << std::endl;
        }
    }
    if (!config.skip_postprocessing && generator->InvalidVertexRangeIfEmpty()) {
        std::vector<SInt> from_before(size, old_vertex_range.first);
        std::vector<SInt> from_after(size, vertex_range.first);
        MPI_Gather(from_before.data(), 1, KAGEN_MPI_SINT, from_before.data(), 1, KAGEN_MPI_SINT, ROOT, comm);
        MPI_Gather(from_after.data(), 1, KAGEN_MPI_SINT, from_after.data(), 1, KAGEN_MPI_SINT, ROOT, comm);

        if (output_error) {
            bool nothing_changed = true;
            for (PEID pe = 0; pe < size; ++pe) {
                if (from_before[pe] != from_after[pe]) {
                    std::cout << "  Changed start of vertex range for PE " << std::setw(std::log10(size) + 1) << pe
                              << " from " << from_before[pe] << " to " << from_after[pe] << std::endl;
                    nothing_changed = false;
                }
            }
            if (nothing_changed) {
                std::cout << "  Vertex ranges where already correct on all PEs" << std::endl;
            }
        }
    }

    // Validation
    if (config.validate_simple_graph) {
        bool success = ValidateSimpleGraph(edges, vertex_range, comm);
        MPI_Allreduce(MPI_IN_PLACE, &success, 1, MPI_C_BOOL, MPI_LOR, comm);
        if (!success) {
            if (output_error) {
                std::cerr << "Simple graph validation failed\n";
            }
            std::exit(1);
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
