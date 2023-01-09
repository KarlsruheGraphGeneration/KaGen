#include "kagen.h"

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include "kagen/context.h"
#include "kagen/facade.h"

namespace kagen {
namespace {}

KaGen::KaGen(MPI_Comm comm) : comm_(comm), config_(std::make_unique<PGeneratorConfig>()) {
    SetDefaults();
}

KaGen::KaGen(KaGen&&) noexcept            = default;
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

void KaGen::EnableOutput(const bool header) {
    config_->quiet        = false;
    config_->print_header = header;
}

void KaGen::UseHPFloats(const bool state) {
    config_->hp_floats = state ? 1 : -1;
}

void KaGen::SetNumberOfChunks(const SInt k) {
    config_->k = k;
}

namespace {
auto GenericGenerateFromOptionString(const std::string& options_str, PGeneratorConfig base_config, MPI_Comm comm) {
    return Generate(CreateConfigFromString(options_str, base_config), comm);
}
} // namespace

KaGenResult KaGen::GenerateFromOptionString(const std::string& options) {
    return GenericGenerateFromOptionString(options, *config_, comm_).tuple();
}

KaGenResult KaGen::GenerateDirectedGNM(const SInt n, const SInt m, const bool self_loops) {
    config_->generator  = GeneratorType::GNM_DIRECTED;
    config_->n          = n;
    config_->m          = m;
    config_->self_loops = self_loops;
    return Generate(*config_, comm_).tuple();
}

KaGenResult KaGen::GenerateUndirectedGNM(const SInt n, const SInt m, const bool self_loops) {
    config_->generator  = GeneratorType::GNM_UNDIRECTED;
    config_->n          = n;
    config_->m          = m;
    config_->self_loops = self_loops;
    return Generate(*config_, comm_).tuple();
}

KaGenResult KaGen::GenerateDirectedGNP(const SInt n, const LPFloat p, const bool self_loops) {
    config_->generator  = GeneratorType::GNP_DIRECTED;
    config_->n          = n;
    config_->p          = p;
    config_->self_loops = self_loops;
    return Generate(*config_, comm_).tuple();
}

KaGenResult KaGen::GenerateUndirectedGNP(const SInt n, const LPFloat p, const bool self_loops) {
    config_->generator  = GeneratorType::GNP_UNDIRECTED;
    config_->n          = n;
    config_->p          = p;
    config_->self_loops = self_loops;
    return Generate(*config_, comm_).tuple();
}

namespace {
KaGenResult GenerateRGG2D_Impl(
    PGeneratorConfig& config, const SInt n, const SInt m, const LPFloat r, const bool coordinates, MPI_Comm comm) {
    config.generator   = GeneratorType::RGG_2D;
    config.n           = n;
    config.m           = m;
    config.r           = r;
    config.coordinates = coordinates;
    return Generate(config, comm).tuple();
}
} // namespace

KaGenResult KaGen::GenerateRGG2D(const SInt n, const LPFloat r, const bool coordinates) {
    return GenerateRGG2D_Impl(*config_, n, 0, r, coordinates, comm_);
}

KaGenResult KaGen::GenerateRGG2D_NM(const SInt n, const SInt m, const bool coordinates) {
    return GenerateRGG2D_Impl(*config_, n, m, 0.0, coordinates, comm_);
}

KaGenResult KaGen::GenerateRGG2D_MR(const SInt m, const LPFloat r, const bool coordinates) {
    return GenerateRGG2D_Impl(*config_, 0, m, r, coordinates, comm_);
}

namespace {
KaGenResult GenerateRGG3D_Impl(
    PGeneratorConfig& config, const SInt n, const SInt m, const LPFloat r, const bool coordinates, MPI_Comm comm) {
    config.generator   = GeneratorType::RGG_3D;
    config.m           = m;
    config.n           = n;
    config.r           = r;
    config.coordinates = coordinates;
    return Generate(config, comm).tuple();
}
} // namespace

KaGenResult KaGen::GenerateRGG3D(const SInt n, const LPFloat r, const bool coordinates) {
    return GenerateRGG3D_Impl(*config_, n, 0, r, coordinates, comm_);
}

KaGenResult KaGen::GenerateRGG3D_NM(const SInt n, const SInt m, const bool coordinates) {
    return GenerateRGG3D_Impl(*config_, n, m, 0.0, coordinates, comm_);
}

KaGenResult KaGen::GenerateRGG3D_MR(const SInt m, const LPFloat r, const bool coordinates) {
    return GenerateRGG3D_Impl(*config_, 0, m, r, coordinates, comm_);
}

#ifdef KAGEN_CGAL_FOUND
namespace {
KaGenResult GenerateRDG2D_Impl(
    PGeneratorConfig& config, const SInt n, const SInt m, const bool periodic, const bool coordinates, MPI_Comm comm) {
    config.generator   = GeneratorType::RDG_2D;
    config.n           = n;
    config.m           = m;
    config.periodic    = periodic;
    config.coordinates = coordinates;
    return Generate(config, comm).tuple();
}
} // namespace

KaGenResult KaGen::GenerateRDG2D(const SInt n, const bool periodic, const bool coordinates) {
    return GenerateRDG2D_Impl(*config_, n, 0, periodic, coordinates, comm_);
}

KaGenResult KaGen::GenerateRDG2D_M(const SInt m, const bool periodic, const bool coordinates) {
    return GenerateRDG2D_Impl(*config_, 0, m, periodic, coordinates, comm_);
}

namespace {
KaGenResult
GenerateRDG3D_Impl(PGeneratorConfig& config, const SInt n, const SInt m, const bool coordinates, MPI_Comm comm) {
    config.generator   = GeneratorType::RDG_3D;
    config.n           = n;
    config.m           = m;
    config.coordinates = coordinates;
    return Generate(config, comm).tuple();
}
} // namespace

KaGenResult KaGen::GenerateRDG3D(const SInt n, const bool coordinates) {
    return GenerateRDG3D_Impl(*config_, n, 0, coordinates, comm_);
}

KaGenResult KaGen::GenerateRDG3D_M(const SInt m, const bool coordinates) {
    return GenerateRDG3D_Impl(*config_, 0, m, coordinates, comm_);
}
#else  // KAGEN_CGAL_FOUND
KaGenResult KaGen::GenerateRDG2D(SInt, bool, bool) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}

KaGenResult KaGen::GenerateRDG2D_M(SInt, bool, bool) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}

KaGenResult KaGen::GenerateRDG3D(SInt, bool) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}

KaGenResult KaGen::GenerateRDG3D_M(SInt, bool) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}
#endif // KAGEN_CGAL_FOUND

namespace {
KaGenResult GenerateBA_Impl(
    PGeneratorConfig& config, const SInt n, const SInt m, const SInt d, const bool directed, const bool self_loops,
    MPI_Comm comm) {
    config.generator  = GeneratorType::BA;
    config.m          = m;
    config.n          = n;
    config.min_degree = d;
    config.self_loops = self_loops;
    config.directed   = directed;
    return Generate(config, comm).tuple();
}
} // namespace

KaGenResult KaGen::GenerateBA(const SInt n, const SInt d, const bool directed, const bool self_loops) {
    return GenerateBA_Impl(*config_, n, 0, d, directed, self_loops, comm_);
}

KaGenResult KaGen::GenerateBA_NM(const SInt n, const SInt m, const bool directed, const bool self_loops) {
    return GenerateBA_Impl(*config_, n, m, 0.0, directed, self_loops, comm_);
}

KaGenResult KaGen::GenerateBA_MD(const SInt m, const SInt d, const bool directed, const bool self_loops) {
    return GenerateBA_Impl(*config_, 0, m, d, directed, self_loops, comm_);
}

namespace {
KaGenResult GenerateRHG_Impl(
    PGeneratorConfig& config, const LPFloat gamma, const SInt n, const SInt m, const LPFloat d, const bool coordinates,
    MPI_Comm comm) {
    config.generator   = GeneratorType::RHG;
    config.n           = n;
    config.m           = m;
    config.avg_degree  = d;
    config.plexp       = gamma;
    config.coordinates = coordinates;
    return Generate(config, comm).tuple();
}
} // namespace

KaGenResult KaGen::GenerateRHG(const LPFloat gamma, const SInt n, const LPFloat d, const bool coordinates) {
    return GenerateRHG_Impl(*config_, gamma, n, 0, d, coordinates, comm_);
}

KaGenResult KaGen::GenerateRHG_NM(const LPFloat gamma, const SInt n, const SInt m, const bool coordinates) {
    return GenerateRHG_Impl(*config_, gamma, n, m, 0.0, coordinates, comm_);
}

KaGenResult KaGen::GenerateRHG_MD(const LPFloat gamma, const SInt m, const LPFloat d, const bool coordinates) {
    return GenerateRHG_Impl(*config_, gamma, 0, m, d, coordinates, comm_);
}

namespace {
KaGenResult GenerateGrid2D_Impl(
    PGeneratorConfig& config, const SInt grid_x, const SInt grid_y, const LPFloat p, const SInt m, const bool periodic,
    const bool coordinates, MPI_Comm comm) {
    config.generator   = GeneratorType::GRID_2D;
    config.grid_x      = grid_x;
    config.grid_y      = grid_y;
    config.p           = p;
    config.m           = m;
    config.periodic    = periodic;
    config.coordinates = coordinates;
    return Generate(config, comm).tuple();
}
} // namespace

KaGenResult KaGen::GenerateGrid2D(
    const SInt grid_x, const SInt grid_y, const LPFloat p, const bool periodic, const bool coordinates) {
    return GenerateGrid2D_Impl(*config_, grid_x, grid_y, p, 0, periodic, coordinates, comm_);
}

KaGenResult KaGen::GenerateGrid2D_N(const SInt n, const LPFloat p, const bool periodic, const bool coordinates) {
    const SInt sqrt_n = std::sqrt(n);
    return GenerateGrid2D(sqrt_n, sqrt_n, p, periodic, coordinates);
}

KaGenResult KaGen::GenerateGrid2D_NM(const SInt n, const SInt m, const bool periodic, const bool coordinates) {
    const SInt sqrt_n = std::sqrt(n);
    return GenerateGrid2D_Impl(*config_, sqrt_n, sqrt_n, 0.0, m, periodic, coordinates, comm_);
}

namespace {
KaGenResult GenerateGrid3D_Impl(
    PGeneratorConfig& config, const SInt grid_x, const SInt grid_y, const SInt grid_z, const LPFloat p, const SInt m,
    const bool periodic, const bool coordinates, MPI_Comm comm) {
    config.generator   = GeneratorType::GRID_3D;
    config.grid_x      = grid_x;
    config.grid_y      = grid_y;
    config.grid_z      = grid_z;
    config.p           = p;
    config.m           = m;
    config.periodic    = periodic;
    config.coordinates = coordinates;
    return Generate(config, comm).tuple();
}
} // namespace

KaGenResult KaGen::GenerateGrid3D(
    const SInt grid_x, const SInt grid_y, const SInt grid_z, const LPFloat p, const bool periodic,
    const bool coordinates) {
    return GenerateGrid3D_Impl(*config_, grid_x, grid_y, grid_z, p, 0, periodic, coordinates, comm_);
}

KaGenResult KaGen::GenerateGrid3D_N(const SInt n, const LPFloat p, const bool periodic, const bool coordinates) {
    const SInt cbrt_n = std::cbrt(n);
    return GenerateGrid3D(cbrt_n, cbrt_n, cbrt_n, p, periodic, coordinates);
}

KaGenResult KaGen::GenerateGrid3D_NM(const SInt n, const SInt m, const bool periodic, const bool coordinates) {
    const SInt cbrt_n = std::cbrt(n);
    return GenerateGrid3D_Impl(*config_, cbrt_n, cbrt_n, cbrt_n, 0.0, m, periodic, coordinates, comm_);
}

KaGenResult KaGen::GenerateKronecker(const SInt n, const SInt m, const bool directed, const bool self_loops) {
    config_->generator  = GeneratorType::KRONECKER;
    config_->n          = n;
    config_->m          = m;
    config_->directed   = directed;
    config_->self_loops = self_loops;
    return Generate(*config_, comm_).tuple();
}

KaGenResult KaGen::GenerateRMAT(
    const SInt n, const SInt m, const LPFloat a, const LPFloat b, const LPFloat c, const bool directed,
    const bool self_loops) {
    config_->generator  = GeneratorType::RMAT;
    config_->n          = n;
    config_->m          = m;
    config_->rmat_a     = a;
    config_->rmat_b     = b;
    config_->rmat_c     = c;
    config_->directed   = directed;
    config_->self_loops = self_loops;
    return Generate(*config_, comm_).tuple();
}

void KaGen::SetDefaults() {
    config_->quiet         = true;
    config_->output_format = GraphFormat::NONE; // ignored anyways
    // keep all other defaults
}
} // namespace kagen
