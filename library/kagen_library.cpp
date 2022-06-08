#include "kagen_library.h"

#include <cmath>

#include "kagen/context.h"
#include "kagen/facade.h"

namespace kagen {
namespace {}

KaGen::KaGen(MPI_Comm comm) : comm_(comm), config_(std::make_unique<PGeneratorConfig>()) {
    SetDefaults();
}

KaGen::KaGen(KaGen&&) noexcept = default;
KaGen& KaGen::operator=(KaGen&&) noexcept = default;

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

KaGenResult KaGen::GenerateRGG2D(const SInt n, const LPFloat r) {
    config_->generator = GeneratorType::RGG_2D;
    config_->n         = n;
    config_->r         = r;
    return Generate(*config_, comm_);
}

KaGenResult KaGen::GenerateRGG2DNM(const SInt n, const SInt m) {
    config_->generator = GeneratorType::RGG_2D;
    config_->n         = n;
    config_->m         = m;
    return Generate(*config_, comm_);
}

KaGenResult KaGen::GenerateRGG2DMR(const SInt m, const LPFloat r) {
    config_->generator = GeneratorType::RGG_2D;
    config_->m         = m;
    config_->r         = r;
    return Generate(*config_, comm_);
}

KaGenResult KaGen::GenerateRGG3D(const SInt n, const LPFloat r) {
    config_->generator = GeneratorType::RGG_3D;
    config_->n         = n;
    config_->r         = r;
    return Generate(*config_, comm_);
}

KaGenResult KaGen::GenerateRGG3DNM(const SInt n, const SInt m) {
    config_->generator = GeneratorType::RGG_3D;
    config_->n         = n;
    config_->m         = m;
    return Generate(*config_, comm_);
}

KaGenResult KaGen::GenerateRGG3DMR(const SInt m, const LPFloat r) {
    config_->generator = GeneratorType::RGG_3D;
    config_->m         = m;
    config_->r         = r;
    return Generate(*config_, comm_);
}

#ifdef KAGEN_CGAL_FOUND
KaGenResult KaGen::GenerateRDG2D(const SInt n) {
    config_->generator = GeneratorType::RDG_2D;
    config_->n         = n;
    return Generate(*config_, comm_);
}

KaGenResult KaGen::GenerateRDG3D(const SInt n) {
    config_->generator = GeneratorType::RDG_3D;
    config_->n         = n;
    return Generate(*config_, comm_);
}
#else  // KAGEN_CGAL_FOUND
KaGenResult KaGen::GenerateRDG2D(SInt) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}

KaGenResult KaGen::GenerateRDG3D(SInt) {
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

namespace {
KaGenResult2D GenerateRHGImpl(
    PGeneratorConfig& config, const SInt n, const LPFloat gamma, const SInt d, const bool coordinates, MPI_Comm comm) {
    config.generator   = GeneratorType::RHG;
    config.n           = n;
    config.plexp       = gamma;
    config.avg_degree  = d;
    config.coordinates = coordinates;
    return Generate(config, comm);
}
} // namespace

KaGenResult KaGen::GenerateRHG(const SInt n, const LPFloat gamma, const SInt d) {
    return GenerateRHGImpl(*config_, n, gamma, d, false, comm_);
}

KaGenResult2D KaGen::GenerateRHGWithCoordinates(const SInt n, const LPFloat gamma, const SInt d) {
    return GenerateRHGImpl(*config_, n, gamma, d, true, comm_);
}

KaGenResult KaGen::GenerateGrid2D(const SInt n, const LPFloat p, const bool periodic) {
    const SInt sqrt_n = std::sqrt(n);
    return GenerateGrid2D(sqrt_n, sqrt_n, p, periodic);
}

KaGenResult KaGen::GenerateGrid2D(const SInt grid_x, const SInt grid_y, const LPFloat p, const bool periodic) {
    config_->generator = GeneratorType::GRID_2D;
    config_->grid_x    = grid_x;
    config_->grid_y    = grid_y;
    config_->p         = p;
    config_->periodic  = periodic;
    return Generate(*config_, comm_);
}

KaGenResult KaGen::GenerateGrid3D(const SInt n, const LPFloat p, const bool periodic) {
    const SInt cbrt_n = std::cbrt(n);
    return GenerateGrid3D(cbrt_n, cbrt_n, cbrt_n, p, periodic);
}

KaGenResult
KaGen::GenerateGrid3D(const SInt grid_x, const SInt grid_y, const SInt grid_z, const LPFloat p, const bool periodic) {
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
