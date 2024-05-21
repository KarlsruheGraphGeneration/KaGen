#include "kagen/factories.h"

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
} // namespace kagen
