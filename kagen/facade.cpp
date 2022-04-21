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

std::pair<EdgeList, VertexRange> Generate(const PGeneratorConfig& config) {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return Generate(config, rank, size);
}

std::pair<EdgeList, VertexRange> Generate(const PGeneratorConfig& config, const PEID rank, const PEID size) {
    auto generator             = CreateGenerator(config, rank, size);
    auto [edges, vertex_range] = generator->Generate();

    if ((generator->Features() & GeneratorFeature::ALMOST_UNDIRECTED) != GeneratorFeature::NONE) {
        AddReverseEdges(edges, vertex_range);
    }

    if (config.validate_simple_graph) {
        ValidateSimpleGraph(edges, vertex_range);
    }

    return {std::move(edges), vertex_range};
}
} // namespace kagen
