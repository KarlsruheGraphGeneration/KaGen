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

namespace {
KaGenResult2D GenerateRGG2D_Impl(
    PGeneratorConfig& config, const SInt n, const SInt m, const LPFloat r, const bool coordinates, MPI_Comm comm) {
    config.generator   = GeneratorType::RGG_2D;
    config.n           = n;
    config.m           = m;
    config.r           = r;
    config.coordinates = coordinates;
    return Generate(config, comm);
}
} // namespace

KaGenResult KaGen::GenerateRGG2D(const SInt n, const LPFloat r) {
    return GenerateRGG2D_Impl(*config_, n, 0, r, false, comm_);
}

KaGenResult KaGen::GenerateRGG2D_NM(const SInt n, const SInt m) {
    return GenerateRGG2D_Impl(*config_, n, m, 0.0, false, comm_);
}

KaGenResult KaGen::GenerateRGG2D_MR(const SInt m, const LPFloat r) {
    return GenerateRGG2D_Impl(*config_, 0, m, r, false, comm_);
}

KaGenResult2D KaGen::GenerateRGG2D_Coordinates(const SInt n, const LPFloat r) {
    return GenerateRGG2D_Impl(*config_, n, 0, r, true, comm_);
}

namespace {
KaGenResult3D GenerateRGG3D_Impl(
    PGeneratorConfig& config, const SInt n, const SInt m, const LPFloat r, const bool coordinates, MPI_Comm comm) {
    config.generator   = GeneratorType::RGG_3D;
    config.m           = m;
    config.n           = n;
    config.r           = r;
    config.coordinates = coordinates;
    return Generate(config, comm);
}
} // namespace

KaGenResult KaGen::GenerateRGG3D(const SInt n, const LPFloat r) {
    return GenerateRGG3D_Impl(*config_, n, 0, r, false, comm_);
}

KaGenResult KaGen::GenerateRGG3D_NM(const SInt n, const SInt m) {
    return GenerateRGG3D_Impl(*config_, n, m, 0.0, false, comm_);
}

KaGenResult KaGen::GenerateRGG3D_MR(const SInt m, const LPFloat r) {
    return GenerateRGG3D_Impl(*config_, 0, m, r, false, comm_);
}

KaGenResult3D KaGen::GenerateRGG3D_Coordinates(const SInt n, const LPFloat r) {
    return GenerateRGG3D_Impl(*config_, n, 0, r, true, comm_);
}

#ifdef KAGEN_CGAL_FOUND
namespace {
KaGenResult2D GenerateRDG2D_Impl(
    PGeneratorConfig& config, const SInt n, const SInt m, const bool periodic, const bool coordinates, MPI_Comm comm) {
    config.generator   = GeneratorType::RDG_2D;
    config.n           = n;
    config.m           = m;
    config.periodic    = periodic;
    config.coordinates = coordinates;
    return Generate(config, comm);
}
} // namespace

KaGenResult KaGen::GenerateRDG2D(const SInt n, const bool periodic) {
    return GenerateRDG2D_Impl(*config_, n, 0, periodic, false, comm_);
}

KaGenResult KaGen::GenerateRDG2D_M(const SInt m, const bool periodic) {
    return GenerateRDG2D_Impl(*config_, 0, m, periodic, false, comm_);
}

KaGenResult2D KaGen::GenerateRDG2D_Coordinates(const SInt n, const bool periodic) {
    return GenerateRDG2D_Impl(*config_, n, 0, periodic, true, comm_);
}

namespace {
KaGenResult3D
GenerateRDG3D_Impl(PGeneratorConfig& config, const SInt n, const SInt m, const bool coordinates, MPI_Comm comm) {
    config.generator   = GeneratorType::RDG_3D;
    config.n           = n;
    config.m           = m;
    config.coordinates = coordinates;
    return Generate(config, comm);
}
} // namespace

KaGenResult KaGen::GenerateRDG3D(const SInt n) {
    return GenerateRDG3D_Impl(*config_, n, 0, false, comm_);
}

KaGenResult KaGen::GenerateRDG3D_M(const SInt m) {
    return GenerateRDG3D_Impl(*config_, 0, m, false, comm_);
}

KaGenResult3D KaGen::GenerateRDG3D_Coordinates(const SInt n) {
    return GenerateRDG3D_Impl(*config_, n, 0, true, comm_);
}
#else  // KAGEN_CGAL_FOUND
KaGenResult KaGen::GenerateRDG2D(SInt) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}

KaGenResult KaGen::GenerateRDG2D_M(SInt) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}

KaGenResult2D KaGen::GenerateRDG2D_Coordinates(SInt) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}

KaGenResult KaGen::GenerateRDG3D(SInt) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}

KaGenResult KaGen::GenerateRDG3D_M(SInt) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}

KaGenResult3D KaGen::GenerateRDG3D_Coordinates(SInt) {
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
KaGenResult2D GenerateRHG_Impl(
    PGeneratorConfig& config, const LPFloat gamma, const SInt n, const SInt m, const LPFloat d, const bool coordinates,
    MPI_Comm comm) {
    config.generator   = GeneratorType::RHG;
    config.n           = n;
    config.m           = m;
    config.avg_degree  = d;
    config.plexp       = gamma;
    config.coordinates = coordinates;
    return Generate(config, comm);
}
} // namespace

KaGenResult KaGen::GenerateRHG(const LPFloat gamma, const SInt n, const LPFloat d) {
    return GenerateRHG_Impl(*config_, gamma, n, 0, d, false, comm_);
}

KaGenResult KaGen::GenerateRHG_NM(const LPFloat gamma, const SInt n, const SInt m) {
    return GenerateRHG_Impl(*config_, gamma, n, m, 0.0, false, comm_);
}

KaGenResult KaGen::GenerateRHG_MD(const LPFloat gamma, const SInt m, const LPFloat d) {
    return GenerateRHG_Impl(*config_, gamma, 0, m, d, false, comm_);
}

KaGenResult2D KaGen::GenerateRHG_Coordinates(const LPFloat gamma, const SInt n, const LPFloat d) {
    return GenerateRHG_Impl(*config_, gamma, n, 0, d, true, comm_);
}

namespace {
KaGenResult2D GenerateGrid2D_Impl(
    PGeneratorConfig& config, const SInt grid_x, const SInt grid_y, const LPFloat p, const bool periodic,
    const bool coordinates, MPI_Comm comm) {
    config.generator   = GeneratorType::GRID_2D;
    config.grid_x      = grid_x;
    config.grid_y      = grid_y;
    config.p           = p;
    config.periodic    = periodic;
    config.coordinates = coordinates;
    return Generate(config, comm);
}
} // namespace

KaGenResult KaGen::GenerateGrid2D(const SInt grid_x, const SInt grid_y, const LPFloat p, const bool periodic) {
    return GenerateGrid2D_Impl(*config_, grid_x, grid_y, p, periodic, false, comm_);
}

KaGenResult KaGen::GenerateGrid2D_N(const SInt n, const LPFloat p, const bool periodic) {
    const SInt sqrt_n = std::sqrt(n);
    return GenerateGrid2D(sqrt_n, sqrt_n, p, periodic);
}

KaGenResult2D
KaGen::GenerateGrid2D_Coordinates(const SInt grid_x, const SInt grid_y, const LPFloat p, const bool periodic) {
    return GenerateGrid2D_Impl(*config_, grid_x, grid_y, p, periodic, true, comm_);
}

namespace {
KaGenResult3D GenerateGrid3D_Impl(
    PGeneratorConfig& config, const SInt grid_x, const SInt grid_y, const SInt grid_z, const LPFloat p,
    const bool periodic, const bool coordinates, MPI_Comm comm) {
    config.generator   = GeneratorType::GRID_3D;
    config.grid_x      = grid_x;
    config.grid_y      = grid_y;
    config.grid_z      = grid_z;
    config.p           = p;
    config.periodic    = periodic;
    config.coordinates = coordinates;
    return Generate(config, comm);
}
} // namespace

KaGenResult
KaGen::GenerateGrid3D(const SInt grid_x, const SInt grid_y, const SInt grid_z, const LPFloat p, const bool periodic) {
    return GenerateGrid3D_Impl(*config_, grid_x, grid_y, grid_z, p, periodic, false, comm_);
}

KaGenResult KaGen::GenerateGrid3D_N(const SInt n, const LPFloat p, const bool periodic) {
    const SInt cbrt_n = std::cbrt(n);
    return GenerateGrid3D(cbrt_n, cbrt_n, cbrt_n, p, periodic);
}

KaGenResult3D KaGen::GenerateGrid3D_Coordinates(
    const SInt grid_x, const SInt grid_y, const SInt grid_z, const LPFloat p, const bool periodic) {
    return GenerateGrid3D_Impl(*config_, grid_x, grid_y, grid_z, p, periodic, true, comm_);
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
