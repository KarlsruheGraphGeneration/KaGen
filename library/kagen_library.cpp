#include "kagen_library.h"

#include "kagen/generator_config.h"
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
#include "kagen/postprocessing.h"

#ifdef KAGEN_CGAL_FOUND
    #include "kagen/generators/geometric/delaunay/delaunay_2d.h"
    #include "kagen/generators/geometric/delaunay/delaunay_3d.h"
#endif

namespace kagen {
KaGen::KaGen(const PEID rank, const PEID size)
    : rank_(rank),
      size_(size),
      config_(std::make_unique<PGeneratorConfig>()),
      validate_undirected_graph_(false) {
    SetDefaults();
}

KaGen::~KaGen() = default;

void KaGen::SetSeed(const int seed) {
    config_->seed = seed;
}

void KaGen::EnableUndirectedGraphVerification() {
    validate_undirected_graph_ = true;
}

KaGenResult KaGen::GenerateDirectedGMM(const SInt n, const SInt m, const SInt k, const bool self_loops) {
    // Update config
    config_->n          = n;
    config_->m          = m;
    config_->k          = (k == 0 ? config_->k : k);
    config_->self_loops = self_loops;

    // Init and run generator
    GNMDirected gen(*config_, rank_, size_);
    gen.Generate();

    return {std::move(gen.IO().GetEdges()), gen.GetVertexRange()};
}

// Normalized output format
KaGenResult KaGen::GenerateUndirectedGNM(const SInt n, const SInt m, const SInt k, const bool self_loops) {
    // Update config
    config_->n          = n;
    config_->m          = m;
    config_->k          = (k == 0 ? config_->k : k);
    config_->self_loops = self_loops;

    // Init and run generator
    GNMUndirected gen(*config_, rank_, size_);
    gen.Generate();

    if (validate_undirected_graph_) {
        Postprocess(Postprocessing::VALIDATE_UNDIRECTED, gen);
    }

    return {std::move(gen.IO().GetEdges()), gen.GetVertexRange()};
}

KaGenResult KaGen::GenerateDirectedGNP(const SInt n, const LPFloat p, const SInt k, const bool self_loops) {
    // Update config
    config_->n          = n;
    config_->p          = p;
    config_->k          = (k == 0 ? config_->k : k);
    config_->self_loops = self_loops;

    // Init and run generator
    GNPDirected gen(*config_, rank_, size_);
    gen.Generate();

    return {std::move(gen.IO().GetEdges()), gen.GetVertexRange()};
}

KaGenResult KaGen::GenerateUndirectedGNP(const SInt n, const LPFloat p, const SInt k, const bool self_loops) {
    // Update config
    config_->n          = n;
    config_->p          = p;
    config_->k          = (k == 0 ? config_->k : k);
    config_->self_loops = self_loops;

    // Init and run generator
    GNPUndirected gen(*config_, rank_, size_);
    gen.Generate();

    if (validate_undirected_graph_) {
        Postprocess(Postprocessing::VALIDATE_UNDIRECTED, gen);
    }

    return {std::move(gen.IO().GetEdges()), gen.GetVertexRange()};
}

// Normalized output format
KaGenResult KaGen::Generate2DRGG(const SInt n, const LPFloat r, const SInt k) {
    // Update config
    config_->n = n;
    config_->r = r;
    config_->k = (k == 0 ? config_->k : k);

    // Init and run generator
    RGG2D gen(*config_, rank_, size_);
    gen.Generate();

    if (validate_undirected_graph_) {
        Postprocess(Postprocessing::VALIDATE_UNDIRECTED, gen);
    }

    return {std::move(gen.IO().GetEdges()), gen.GetVertexRange()};
}

// Normalized output format
KaGenResult KaGen::Generate3DRGG(const SInt n, const LPFloat r, const SInt k) {
    // Update config
    config_->n = n;
    config_->r = r;
    config_->k = (k == 0 ? config_->k : k);

    // Init and run generator
    RGG3D gen(*config_, rank_, size_);
    gen.Generate();

    if (validate_undirected_graph_) {
        Postprocess(Postprocessing::VALIDATE_UNDIRECTED, gen);
    }

    return {std::move(gen.IO().GetEdges()), gen.GetVertexRange()};
}

#ifdef KAGEN_CGAL_FOUND
KaGenResult KaGen::Generate2DRDG(const SInt n, const SInt k) {
    // Update config
    config_->n = n;
    config_->k = (k == 0 ? config_->k : k);

    // Init and run generator
    Delaunay2D gen(*config_, rank_, size_);
    gen.Generate();

    // Fix broken edge list
    Postprocess(Postprocessing::FIX_UNDIRECTED_EDGE_LIST, gen);

    if (validate_undirected_graph_) {
        Postprocess(Postprocessing::VALIDATE_UNDIRECTED, gen);
    }

    return {std::move(gen.IO().GetEdges()), gen.GetVertexRange()};
}

KaGenResult KaGen::Generate3DRDG(const SInt n, const SInt k) {
    // Update config
    config_->n = n;
    config_->k = (k == 0 ? config_->k : k);

    // Init and run generator
    Delaunay3D gen(*config_, rank_, size_);
    gen.Generate();

    if (validate_undirected_graph_) {
        Postprocess(Postprocessing::VALIDATE_UNDIRECTED, gen);
    }

    return {std::move(gen.IO().GetEdges()), gen.GetVertexRange()};
}
#else  // KAGEN_CGAL_FOUND
KaGenResult KaGen::Generate2DRDG(SInt, SInt) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}

KaGenResult KaGen::Generate3DRDG(SInt, SInt) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}
#endif // KAGEN_CGAL_FOUND

// Normalized output format
KaGenResult KaGen::GenerateBA(const SInt n, const SInt d, const SInt k) {
    // Update config
    config_->n          = n;
    config_->min_degree = d;
    config_->k          = (k == 0 ? config_->k : k);

    // Init and run generator
    Barabassi gen(*config_, rank_, size_);
    gen.Generate();

    // Redistribute graph such that each PE owns a consecutive set of vertices
    Postprocess(Postprocessing::REDISTRIBUTE_GRAPH, gen);

    if (validate_undirected_graph_) {
        Postprocess(Postprocessing::VALIDATE_UNDIRECTED, gen);
    }

    return {std::move(gen.IO().GetEdges()), gen.GetVertexRange()};
}

// Normalized output format
KaGenResult KaGen::GenerateRHG(const SInt n, const LPFloat gamma, const SInt d, const SInt k) {
    // Update config
    config_->n          = n;
    config_->plexp      = gamma;
    config_->avg_degree = d;
    config_->k          = (k == 0 ? config_->k : k);

    // Init and run generator
    Hyperbolic gen(*config_, rank_, size_);
    gen.Generate();

    // Fix broken edge list
    Postprocess(Postprocessing::FIX_UNDIRECTED_EDGE_LIST, gen);

    if (validate_undirected_graph_) {
        Postprocess(Postprocessing::VALIDATE_UNDIRECTED, gen);
    }

    return {std::move(gen.IO().GetEdges()), gen.GetVertexRange()};
}

KaGenResult KaGen::Generate2DGrid(const SInt n, const LPFloat p, const bool periodic, const SInt k) {
    SInt sqrt_n = std::sqrt(n);
    return Generate2DGrid(sqrt_n, sqrt_n, p, periodic, k);
}

KaGenResult
KaGen::Generate2DGrid(const SInt grid_x, const SInt grid_y, const LPFloat p, const bool periodic, const SInt k) {
    // Update config
    config_->grid_x   = grid_x;
    config_->grid_y   = grid_y;
    config_->p        = p;
    config_->periodic = periodic;
    config_->k        = (k == 0 ? config_->k : k);

    // Init and run generator
    Grid2D gen(*config_, rank_, size_);
    gen.Generate();

    if (validate_undirected_graph_) {
        Postprocess(Postprocessing::VALIDATE_UNDIRECTED, gen);
    }

    return {std::move(gen.IO().GetEdges()), gen.GetVertexRange()};
}

KaGenResult KaGen::Generate3DGrid(
    const SInt grid_x, const SInt grid_y, const SInt grid_z, const LPFloat p, const bool periodic, const SInt k) {
    // Update config
    config_->grid_x   = grid_x;
    config_->grid_y   = grid_y;
    config_->grid_z   = grid_z;
    config_->p        = p;
    config_->periodic = periodic;
    config_->k        = (k == 0 ? config_->k : k);

    // Init and run generator
    Grid3D gen(*config_, rank_, size_);
    gen.Generate();

    if (validate_undirected_graph_) {
        Postprocess(Postprocessing::VALIDATE_UNDIRECTED, gen);
    }

    return {std::move(gen.IO().GetEdges()), gen.GetVertexRange()};
}

// @todo broken generator
KaGenResult KaGen::GenerateKronecker(const SInt n, const SInt m, const SInt k) {
    // Update config
    config_->n = n;
    config_->m = m;
    config_->k = (k == 0 ? config_->k : k);

    // Init and run generator
    Kronecker gen(*config_, rank_, size_);
    gen.Generate();

    if (validate_undirected_graph_) {
        Postprocess(Postprocessing::VALIDATE_UNDIRECTED, gen);
    }

    return {std::move(gen.IO().GetEdges()), gen.GetVertexRange()};
}

void KaGen::SetDefaults() {
    config_->n            = 100;
    config_->m            = 0;
    config_->k            = size_;
    config_->seed         = 1;
    config_->hash_sample  = false;
    config_->use_binom    = false;
    config_->output_file  = "out";
    config_->debug_output = "dbg";
    config_->dist_size    = 10;
    config_->p            = 0.0;
    config_->self_loops   = false;
    config_->r            = 0.125;
    config_->avg_degree   = 5.0;
    config_->plexp        = 2.6;
    config_->thres        = 0;
    config_->query_both   = false;
    config_->min_degree   = 4;
    config_->precision    = 32;
    config_->base_size    = (SInt)1 << 8;
    config_->hyp_base     = (SInt)1 << 8;
    config_->iterations   = 1;
}
} // namespace kagen
