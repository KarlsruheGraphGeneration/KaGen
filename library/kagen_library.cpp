#include "kagen_library.h"

#include <cmath>

#include "kagen/facade.h"
#include "kagen/generator_config.h"

namespace kagen {
KaGen::KaGen(const PEID rank, const PEID size)
    : rank_(rank),
      size_(size),
      config_(std::make_unique<PGeneratorConfig>()) {
    SetDefaults();
}

KaGen::~KaGen() = default;

void KaGen::SetSeed(const int seed) {
    config_->seed = seed;
}

void KaGen::EnableUndirectedGraphVerification() {
    config_->validate_undirected_graph = true;
}

KaGenResult KaGen::GenerateDirectedGMM(const SInt n, const SInt m, const SInt k, const bool self_loops) {
    config_->generator  = GeneratorType::GNM_DIRECTED;
    config_->n          = n;
    config_->m          = m;
    config_->k          = (k == 0 ? config_->k : k);
    config_->self_loops = self_loops;
    return Generate(*config_, rank_, size_);
}

KaGenResult KaGen::GenerateUndirectedGNM(const SInt n, const SInt m, const SInt k, const bool self_loops) {
    config_->generator  = GeneratorType::GNM_UNDIRECTED;
    config_->n          = n;
    config_->m          = m;
    config_->k          = (k == 0 ? config_->k : k);
    config_->self_loops = self_loops;
    return Generate(*config_, rank_, size_);
}

KaGenResult KaGen::GenerateDirectedGNP(const SInt n, const LPFloat p, const SInt k, const bool self_loops) {
    config_->generator  = GeneratorType::GNP_DIRECTED;
    config_->n          = n;
    config_->p          = p;
    config_->k          = (k == 0 ? config_->k : k);
    config_->self_loops = self_loops;
    return Generate(*config_, rank_, size_);
}

KaGenResult KaGen::GenerateUndirectedGNP(const SInt n, const LPFloat p, const SInt k, const bool self_loops) {
    config_->generator  = GeneratorType::GNP_UNDIRECTED;
    config_->n          = n;
    config_->p          = p;
    config_->k          = (k == 0 ? config_->k : k);
    config_->self_loops = self_loops;
    return Generate(*config_, rank_, size_);
}

KaGenResult KaGen::Generate2DRGG(const SInt n, const LPFloat r, const SInt k) {
    config_->generator = GeneratorType::RGG_2D;
    config_->n         = n;
    config_->r         = r;
    config_->k         = (k == 0 ? config_->k : k);
    return Generate(*config_, rank_, size_);
}

KaGenResult KaGen::Generate3DRGG(const SInt n, const LPFloat r, const SInt k) {
    config_->generator = GeneratorType::RGG_3D;
    config_->n         = n;
    config_->r         = r;
    config_->k         = (k == 0 ? config_->k : k);
    return Generate(*config_, rank_, size_);
}

#ifdef KAGEN_CGAL_FOUND
KaGenResult KaGen::Generate2DRDG(const SInt n, const SInt k) {
    config_->generator = GeneratorType::RDG_2D;
    config_->n         = n;
    config_->k         = (k == 0 ? config_->k : k);
    return Generate(*config_, rank_, size_);
}

KaGenResult KaGen::Generate3DRDG(const SInt n, const SInt k) {
    config_->generator = GeneratorType::RDG_3D;
    config_->n         = n;
    config_->k         = (k == 0 ? config_->k : k);
    return Generate(*config_, rank_, size_);
}
#else  // KAGEN_CGAL_FOUND
KaGenResult KaGen::Generate2DRDG(SInt, SInt) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}

KaGenResult KaGen::Generate3DRDG(SInt, SInt) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}
#endif // KAGEN_CGAL_FOUND

KaGenResult KaGen::GenerateBA(const SInt n, const SInt d, const SInt k) {
    config_->generator  = GeneratorType::BA;
    config_->n          = n;
    config_->min_degree = d;
    config_->k          = (k == 0 ? config_->k : k);
    return Generate(*config_, rank_, size_);
}

// Normalized output format
KaGenResult KaGen::GenerateRHG(const SInt n, const LPFloat gamma, const SInt d, const SInt k) {
    config_->generator  = GeneratorType::RHG;
    config_->n          = n;
    config_->plexp      = gamma;
    config_->avg_degree = d;
    config_->k          = (k == 0 ? config_->k : k);
    return Generate(*config_, rank_, size_);
}

KaGenResult KaGen::Generate2DGrid(const SInt n, const LPFloat p, const bool periodic, const SInt k) {
    SInt sqrt_n = std::sqrt(n);
    return Generate2DGrid(sqrt_n, sqrt_n, p, periodic, k);
}

KaGenResult
KaGen::Generate2DGrid(const SInt grid_x, const SInt grid_y, const LPFloat p, const bool periodic, const SInt k) {
    config_->generator = GeneratorType::GRID_2D;
    config_->grid_x    = grid_x;
    config_->grid_y    = grid_y;
    config_->p         = p;
    config_->periodic  = periodic;
    config_->k         = (k == 0 ? config_->k : k);
    return Generate(*config_, rank_, size_);
}

KaGenResult KaGen::Generate3DGrid(
    const SInt grid_x, const SInt grid_y, const SInt grid_z, const LPFloat p, const bool periodic, const SInt k) {
    config_->generator = GeneratorType::GRID_3D;
    config_->grid_x    = grid_x;
    config_->grid_y    = grid_y;
    config_->grid_z    = grid_z;
    config_->p         = p;
    config_->periodic  = periodic;
    config_->k         = (k == 0 ? config_->k : k);
    return Generate(*config_, rank_, size_);
}

// @todo broken generator
KaGenResult KaGen::GenerateKronecker(const SInt n, const SInt m, const SInt k) {
    config_->generator = GeneratorType::KRONECKER;
    config_->n         = n;
    config_->m         = m;
    config_->k         = (k == 0 ? config_->k : k);
    return Generate(*config_, rank_, size_);
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
