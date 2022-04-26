#include "kagen_library.h"

#include <cmath>

#include "kagen/context.h"
#include "kagen/facade.h"

namespace kagen {
namespace {}

KaGen::KaGen(MPI_Comm comm) : comm_(comm), config_(std::make_unique<PGeneratorConfig>()) {
    SetDefaults();
}

KaGen::~KaGen() = default;

void KaGen::SetSeed(const int seed) {
    config_->seed = seed;
}

void KaGen::EnableUndirectedGraphVerification() {
    config_->validate_simple_graph = true;
}

void KaGen::EnableBasicStatistics() {
    config_->statistics_level = StatisticsLevel::BASIC;
    config_->quiet            = false;
}

void KaGen::EnableAdvancedStatistics() {
    config_->statistics_level = StatisticsLevel::ADVANCED;
    config_->quiet            = false;
}

void KaGen::SetNumberOfChunks(const SInt k) {
    config_->k = k;
}

KaGenResult KaGen::GenerateDirectedGMM(const SInt n, const SInt m, const bool self_loops) {
    config_->generator  = GeneratorType::GNM_DIRECTED;
    config_->n          = n;
    config_->m          = m;
    config_->self_loops = self_loops;
    return Generate(*config_, comm_);
}

KaGenResult KaGen::GenerateUndirectedGNM(const SInt n, const SInt m, const bool self_loops) {
    config_->generator  = GeneratorType::GNM_UNDIRECTED;
    config_->n          = n;
    config_->m          = m;
    config_->self_loops = self_loops;
    return Generate(*config_, comm_);
}

KaGenResult KaGen::GenerateDirectedGNP(const SInt n, const LPFloat p, const bool self_loops) {
    config_->generator  = GeneratorType::GNP_DIRECTED;
    config_->n          = n;
    config_->p          = p;
    config_->self_loops = self_loops;
    return Generate(*config_, comm_);
}

KaGenResult KaGen::GenerateUndirectedGNP(const SInt n, const LPFloat p, const bool self_loops) {
    config_->generator  = GeneratorType::GNP_UNDIRECTED;
    config_->n          = n;
    config_->p          = p;
    config_->self_loops = self_loops;
    return Generate(*config_, comm_);
}

KaGenResult KaGen::Generate2DRGG(const SInt n, const LPFloat r) {
    config_->generator = GeneratorType::RGG_2D;
    config_->n         = n;
    config_->r         = r;
    return Generate(*config_, comm_);
}

KaGenResult KaGen::Generate3DRGG(const SInt n, const LPFloat r) {
    config_->generator = GeneratorType::RGG_3D;
    config_->n         = n;
    config_->r         = r;
    return Generate(*config_, comm_);
}

#ifdef KAGEN_CGAL_FOUND
KaGenResult KaGen::Generate2DRDG(const SInt n) {
    config_->generator = GeneratorType::RDG_2D;
    config_->n         = n;
    return Generate(*config_, comm_);
}

KaGenResult KaGen::Generate3DRDG(const SInt n) {
    config_->generator = GeneratorType::RDG_3D;
    config_->n         = n;
    return Generate(*config_, comm_);
}
#else  // KAGEN_CGAL_FOUND
KaGenResult KaGen::Generate2DRDG(SInt) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}

KaGenResult KaGen::Generate3DRDG(SInt) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}
#endif // KAGEN_CGAL_FOUND

// @todo broken generator
KaGenResult KaGen::GenerateBA(const SInt n, const SInt d) {
    config_->generator  = GeneratorType::BA;
    config_->n          = n;
    config_->min_degree = d;
    return Generate(*config_, comm_);
}

KaGenResult KaGen::GenerateRHG(const SInt n, const LPFloat gamma, const SInt d) {
    config_->generator  = GeneratorType::RHG;
    config_->n          = n;
    config_->plexp      = gamma;
    config_->avg_degree = d;
    return Generate(*config_, comm_);
}

KaGenResult KaGen::Generate2DGrid(const SInt n, const LPFloat p, const bool periodic) {
    const SInt sqrt_n = std::sqrt(n);
    return Generate2DGrid(sqrt_n, sqrt_n, p, periodic);
}

KaGenResult KaGen::Generate2DGrid(const SInt grid_x, const SInt grid_y, const LPFloat p, const bool periodic) {
    config_->generator = GeneratorType::GRID_2D;
    config_->grid_x    = grid_x;
    config_->grid_y    = grid_y;
    config_->p         = p;
    config_->periodic  = periodic;
    return Generate(*config_, comm_);
}

KaGenResult KaGen::Generate3DGrid(const SInt n, const LPFloat p, const bool periodic) {
    const SInt cbrt_n = std::cbrt(n);
    return Generate3DGrid(cbrt_n, cbrt_n, cbrt_n, p, periodic);
}

KaGenResult
KaGen::Generate3DGrid(const SInt grid_x, const SInt grid_y, const SInt grid_z, const LPFloat p, const bool periodic) {
    config_->generator = GeneratorType::GRID_3D;
    config_->grid_x    = grid_x;
    config_->grid_y    = grid_y;
    config_->grid_z    = grid_z;
    config_->p         = p;
    config_->periodic  = periodic;
    return Generate(*config_, comm_);
}

// @todo broken generator
KaGenResult KaGen::GenerateKronecker(const SInt n, const SInt m) {
    config_->generator = GeneratorType::KRONECKER;
    config_->n         = n;
    config_->m         = m;
    return Generate(*config_, comm_);
}

void KaGen::SetDefaults() {
    config_->quiet = true;
    // keep all other defaults
}
} // namespace kagen
