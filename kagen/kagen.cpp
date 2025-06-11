#include "kagen.h"

#include "kagen/context.h"
#include "kagen/in_memory_facade.h"
#include "kagen/streaming_facade.h"

#include <cmath>
#include <numeric>
#include <sstream>
#include <unordered_map>

namespace kagen {
std::string BuildDescription() {
    std::stringstream ss;
    ss << "v" << KAGEN_VERSION_MAJOR << "." << KAGEN_VERSION_MINOR << "." << KAGEN_VERSION_PATCH;
    ss << " (";
#ifdef KAGEN_CGAL_FOUND
    ss << "CGal=yes";
#else  // KAGEN_CGAL_FOUND
    ss << "CGal=no";
#endif // KAGEN_CGAL_FOUND
    ss << ", ";
#ifdef KAGEN_SPARSEHASH_FOUND
    ss << "Sparsehash=yes";
#else  // KAGEN_SPARSEHASH_FOUND
    ss << "Sparsehash=no";
#endif // KAGEN_SPARSEHASH_FOUND
    ss << ", ";
#ifdef KAGEN_XXHASH_FOUND
    ss << "xxHash=yes";
#else  // KAGEN_XXHASH_FOUND
    ss << "xxHash=no";
#endif // KAGEN_XXHASH_FOUND
    ss << ", ";
#ifdef KAGEN_MKL_FOUND
    ss << "MKL=yes";
#else  // KAGEN_MKL_FOUND
    ss << "MKL=no";
#endif // KAGEN_MKL_FOUND
    ss << ")";
    return ss.str();
}

std::unordered_map<std::string, FileFormat> GetOutputFormatMap() {
    return {
        {"noop", FileFormat::NOOP},
        {"edgelist", FileFormat::EDGE_LIST},
        {"edgelist-undirected", FileFormat::EDGE_LIST_UNDIRECTED},
        {"binary-edgelist", FileFormat::BINARY_EDGE_LIST},
        {"binary-edgelist-undirected", FileFormat::BINARY_EDGE_LIST_UNDIRECTED},
        {"plain-edgelist", FileFormat::PLAIN_EDGE_LIST},
        {"metis", FileFormat::METIS},
        {"hmetis", FileFormat::HMETIS},
        {"hmetis-directed", FileFormat::HMETIS_DIRECTED},
        {"dot", FileFormat::DOT},
        {"dot-directed", FileFormat::DOT_DIRECTED},
        {"coordinates", FileFormat::COORDINATES},
        {"parhip", FileFormat::PARHIP},
        {"xtrapulp", FileFormat::XTRAPULP},
        {"experimental/hmetis-ep", FileFormat::HMETIS_EP},
        {"experimental/freight-netl", FileFormat::FREIGHT_NETL},
        {"experimental/freight-netl-ep", FileFormat::FREIGHT_NETL_EP},
        {"experimental/weighted-binary-edgelist", FileFormat::WEIGHTED_BINARY_EDGE_LIST},
    };
}

std::unordered_map<std::string, FileFormat> GetInputFormatMap() {
    return {
        {"extension", FileFormat::EXTENSION},
        {"metis", FileFormat::METIS},
        {"parhip", FileFormat::PARHIP},
        {"plain-edgelist", FileFormat::PLAIN_EDGE_LIST},
        {"experimental/weighted-binary-edgelist", FileFormat::WEIGHTED_BINARY_EDGE_LIST},
        {"experimental/netd-are", FileFormat::NETD_ARE},
    };
}

std::ostream& operator<<(std::ostream& out, FileFormat output_format) {
    switch (output_format) {
        case FileFormat::NOOP:
            return out << "noop";

        case FileFormat::EXTENSION:
            return out << "extension";

        case FileFormat::EDGE_LIST:
            return out << "edgelist";

        case FileFormat::EDGE_LIST_UNDIRECTED:
            return out << "edgelist-undirected";

        case FileFormat::BINARY_EDGE_LIST:
            return out << "binary-edgelist";

        case FileFormat::BINARY_EDGE_LIST_UNDIRECTED:
            return out << "binary-edgelist-undirected";

        case FileFormat::PLAIN_EDGE_LIST:
            return out << "plain-edgelist";

        case FileFormat::METIS:
            return out << "metis";

        case FileFormat::HMETIS:
            return out << "hmetis";

        case FileFormat::HMETIS_DIRECTED:
            return out << "hmetis-directed";

        case FileFormat::HMETIS_EP:
            return out << "experimental/hmetis-ep";

        case FileFormat::DOT:
            return out << "dot";

        case FileFormat::DOT_DIRECTED:
            return out << "dot-directed";

        case FileFormat::COORDINATES:
            return out << "coordinates";

        case FileFormat::PARHIP:
            return out << "parhip";

        case FileFormat::XTRAPULP:
            return out << "xtrapulp";

        case FileFormat::FREIGHT_NETL:
            return out << "experimental/freight-netl";

        case FileFormat::FREIGHT_NETL_EP:
            return out << "experimental/freight-netl-ep";

        case FileFormat::WEIGHTED_BINARY_EDGE_LIST:
            return out << "experimental/weighted-binary-edge-list";

        case FileFormat::NETD_ARE:
            return out << "experimental/netd-are";
    }

    return out << "<invalid>";
}

std::unordered_map<std::string, GeneratorType> GetGeneratorTypeMap() {
    return {
        {"gnm-directed", GeneratorType::GNM_DIRECTED},
        {"gnm_directed", GeneratorType::GNM_DIRECTED},
        {"gnm-undirected", GeneratorType::GNM_UNDIRECTED},
        {"gnm_undirected", GeneratorType::GNM_UNDIRECTED},
        {"gnp-directed", GeneratorType::GNP_DIRECTED},
        {"gnp_directed", GeneratorType::GNP_DIRECTED},
        {"gnp-undirected", GeneratorType::GNP_UNDIRECTED},
        {"gnp-undirected", GeneratorType::GNP_UNDIRECTED},
        {"rgg2d", GeneratorType::RGG_2D},
        {"rgg3d", GeneratorType::RGG_3D},
#ifdef KAGEN_CGAL_FOUND
        {"rdg2d", GeneratorType::RDG_2D},
        {"rdg3d", GeneratorType::RDG_3D},
#endif // KAGEN_CGAL_FOUND
        {"grid2d", GeneratorType::GRID_2D},
        {"grid3d", GeneratorType::GRID_3D},
        {"path", GeneratorType::PATH_DIRECTED},
        {"ba", GeneratorType::BA},
        {"kronecker", GeneratorType::KRONECKER},
        {"rhg", GeneratorType::RHG},
        {"rmat", GeneratorType::RMAT},
        {"image", GeneratorType::IMAGE_MESH},
        {"imagemesh", GeneratorType::IMAGE_MESH},
        {"image-mesh", GeneratorType::IMAGE_MESH},
        {"file", GeneratorType::FILE},
        {"static", GeneratorType::FILE}, // @deprecated
    };
}

std::ostream& operator<<(std::ostream& out, GeneratorType generator_type) {
    switch (generator_type) {
        case GeneratorType::GNM_DIRECTED:
            return out << "gnm-directed";

        case GeneratorType::GNM_UNDIRECTED:
            return out << "gnm-undirected";

        case GeneratorType::GNP_DIRECTED:
            return out << "gnp-directed";

        case GeneratorType::GNP_UNDIRECTED:
            return out << "gnp-undirected";

        case GeneratorType::RGG_2D:
            return out << "rgg2d";

        case GeneratorType::RGG_3D:
            return out << "rgg3d";

        case GeneratorType::RDG_2D:
            return out << "rdg2d";

        case GeneratorType::RDG_3D:
            return out << "rdg3d";

        case GeneratorType::GRID_2D:
            return out << "grid2d";

        case GeneratorType::GRID_3D:
            return out << "grid3d";

        case GeneratorType::PATH_DIRECTED:
            return out << "path-directed";

        case GeneratorType::BA:
            return out << "ba";

        case GeneratorType::KRONECKER:
            return out << "kronecker";

        case GeneratorType::RHG:
            return out << "rhg";

        case GeneratorType::RMAT:
            return out << "rmat";

        case GeneratorType::IMAGE_MESH:
            return out << "image-mesh";

        case GeneratorType::FILE:
            return out << "file";
    }

    return out << "<invalid>";
}

bool operator<=(StatisticsLevel a, StatisticsLevel b) {
    return static_cast<std::uint8_t>(a) <= static_cast<std::uint8_t>(b);
}

std::unordered_map<std::string, StatisticsLevel> GetStatisticsLevelMap() {
    return {
        {"none", StatisticsLevel::NONE},
        {"basic", StatisticsLevel::BASIC},
        {"advanced", StatisticsLevel::ADVANCED},
    };
}

std::ostream& operator<<(std::ostream& out, StatisticsLevel statistics_level) {
    switch (statistics_level) {
        case StatisticsLevel::NONE:
            return out << "none";

        case StatisticsLevel::BASIC:
            return out << "basic";

        case StatisticsLevel::ADVANCED:
            return out << "advanced";
    }

    return out << "<invalid>";
}

std::unordered_map<std::string, ImageMeshWeightModel> GetImageMeshWeightModelMap() {
    return {
        {"l2", ImageMeshWeightModel::L2},          {"inv-l2", ImageMeshWeightModel::INV_L2},
        {"ratio", ImageMeshWeightModel::RATIO},    {"inv-ratio", ImageMeshWeightModel::INV_RATIO},
        {"sim", ImageMeshWeightModel::SIMILARITY}, {"similarity", ImageMeshWeightModel::SIMILARITY},
    };
}

std::ostream& operator<<(std::ostream& out, ImageMeshWeightModel weight_model) {
    switch (weight_model) {
        case ImageMeshWeightModel::L2:
            return out << "l2";
        case ImageMeshWeightModel::INV_L2:
            return out << "inv-l2";
        case ImageMeshWeightModel::RATIO:
            return out << "ratio";
        case ImageMeshWeightModel::INV_RATIO:
            return out << "inv-ratio";
        case ImageMeshWeightModel::SIMILARITY:
            return out << "similarity";
    }

    return out << "<invalid>";
}

std::unordered_map<std::string, GraphDistribution> GetGraphDistributionMap() {
    return {
        {"root", GraphDistribution::ROOT},
        {"balance-vertices", GraphDistribution::BALANCE_VERTICES},
        {"balance-edges", GraphDistribution::BALANCE_EDGES},
        {"explicit", GraphDistribution::EXPLICIT},
    };
}

std::ostream& operator<<(std::ostream& out, GraphDistribution distribution) {
    switch (distribution) {
        case GraphDistribution::ROOT:
            return out << "root";
        case GraphDistribution::BALANCE_VERTICES:
            return out << "balance-vertices";
        case GraphDistribution::BALANCE_EDGES:
            return out << "balance-edges";
        case GraphDistribution::EXPLICIT:
            return out << "explicit";
    }

    return out << "<invalid>";
}

std::unordered_map<std::string, VertexWeightGeneratorType> GetVertexWeightGeneratorTypeMap() {
    return {
        {"default", VertexWeightGeneratorType::DEFAULT},
        {"voiding", VertexWeightGeneratorType::VOIDING},
        {"uniform_random", VertexWeightGeneratorType::UNIFORM_RANDOM}};
}

std::ostream& operator<<(std::ostream& out, VertexWeightGeneratorType generator) {
    switch (generator) {
        case kagen::VertexWeightGeneratorType::DEFAULT:
            return out << "default";
        case kagen::VertexWeightGeneratorType::VOIDING:
            return out << "voiding";
        case kagen::VertexWeightGeneratorType::UNIFORM_RANDOM:
            return out << "uniform_random";
    }

    return out << "<invalid>";
}

std::unordered_map<std::string, EdgeWeightGeneratorType> GetEdgeWeightGeneratorTypeMap() {
    return {
        {"default", EdgeWeightGeneratorType::DEFAULT},
        {"voiding", EdgeWeightGeneratorType::VOIDING},
        {"hashing_based", EdgeWeightGeneratorType::HASHING_BASED},
        {"euclidean_distance", EdgeWeightGeneratorType::HASHING_BASED},
        {"uniform_random", EdgeWeightGeneratorType::UNIFORM_RANDOM}};
}

std::ostream& operator<<(std::ostream& out, EdgeWeightGeneratorType generator) {
    switch (generator) {
        case kagen::EdgeWeightGeneratorType::DEFAULT:
            return out << "default";
        case kagen::EdgeWeightGeneratorType::VOIDING:
            return out << "voiding";
        case kagen::EdgeWeightGeneratorType::HASHING_BASED:
            return out << "hashing_based";
        case kagen::EdgeWeightGeneratorType::UNIFORM_RANDOM:
            return out << "uniform_random";
        case kagen::EdgeWeightGeneratorType::EUCLIDEAN_DISTANCE:
            return out << "euclidean_distance";
    }

    return out << "<invalid>";
}

std::unordered_map<std::string, GraphRepresentation> GetGraphRepresentationMap() {
    return {
        {"edge-list", GraphRepresentation::EDGE_LIST},
        {"csr", GraphRepresentation::CSR},
    };
}

std::ostream& operator<<(std::ostream& out, const GraphRepresentation representation) {
    switch (representation) {
        case GraphRepresentation::EDGE_LIST:
            return out << "edge-list";

        case GraphRepresentation::CSR:
            return out << "csr";
    }

    return out << "<invalid>";
}

SInt Graph::NumberOfLocalVertices() const {
    return vertex_range.second - vertex_range.first;
}

SInt Graph::NumberOfGlobalVertices() const {
    SInt number_vertices = NumberOfLocalVertices();
    MPI_Allreduce(MPI_IN_PLACE, &number_vertices, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    return number_vertices;
}

SInt Graph::NumberOfLocalEdges() const {
    switch (representation) {
        case GraphRepresentation::EDGE_LIST:
            return edges.size();

        case GraphRepresentation::CSR:
            return adjncy.size();
    }

    __builtin_unreachable();
}

SInt Graph::NumberOfGlobalEdges() const {
    SInt number_edges = NumberOfLocalEdges();
    MPI_Allreduce(MPI_IN_PLACE, &number_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    return number_edges;
}

void Graph::Clear() {
    edges.clear();
    xadj.clear();
    adjncy.clear();
    vertex_weights.clear();
    edge_weights.clear();
    coordinates.first.clear();
    coordinates.second.clear();
}

void Graph::FreeEdgelist() {
    [[maybe_unused]] auto free_edgelist = std::move(edges);
}

void Graph::FreeCSR() {
    [[maybe_unused]] auto free_xadj   = std::move(xadj);
    [[maybe_unused]] auto free_adjncy = std::move(adjncy);
}

void Graph::SortEdgelist() {
    auto cmp_from = [](const auto& lhs, const auto& rhs) {
        return std::get<0>(lhs) < std::get<0>(rhs);
    };

    if (!std::is_sorted(edges.begin(), edges.end(), cmp_from)) {
        if (!edge_weights.empty()) {
            std::vector<SSInt> indices(edges.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&](const auto& lhs, const auto& rhs) {
                return cmp_from(edges[lhs], edges[rhs]);
            });
            for (std::size_t e = 0; e < edges.size(); ++e) {
                indices[e] = edge_weights[indices[e]];
            }
            std::swap(edge_weights, indices);
        }

        std::sort(edges.begin(), edges.end(), cmp_from);
    }
}

KaGen::KaGen(MPI_Comm comm)
    : comm_(comm),
      config_(std::make_unique<PGeneratorConfig>()),
      representation_(GraphRepresentation::EDGE_LIST) {
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

void KaGen::EnableVertexPermutation() {
    config_->permute = true;
}

void KaGen::ConfigureEdgeWeightGeneration(
    EdgeWeightGeneratorType generator, SInt weight_range_begin, SInt weight_range_end) {
    config_->edge_weights.generator_type     = generator;
    config_->edge_weights.weight_range_begin = weight_range_begin;
    config_->edge_weights.weight_range_end   = weight_range_end;
}

void KaGen::ConfigureVertexWeightGeneration(
    VertexWeightGeneratorType generator, SInt weight_range_begin, SInt weight_range_end) {
    config_->vertex_weights.generator_type     = generator;
    config_->vertex_weights.weight_range_begin = weight_range_begin;
    config_->vertex_weights.weight_range_end   = weight_range_end;
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

void KaGen::UseEdgeListRepresentation() {
    representation_ = GraphRepresentation::EDGE_LIST;
}

void KaGen::UseCSRRepresentation() {
    representation_ = GraphRepresentation::CSR;
}

namespace {
auto GenericGenerateFromOptionString(
    const std::string& options_str, PGeneratorConfig base_config, const GraphRepresentation representation,
    MPI_Comm comm) {
    return GenerateInMemory(CreateConfigFromString(options_str, base_config), representation, comm);
}
} // namespace

Graph KaGen::GenerateFromOptionString(const std::string& options) {
    return GenericGenerateFromOptionString(options, *config_, representation_, comm_);
}

Graph KaGen::GenerateDirectedGNM(const SInt n, const SInt m, const bool self_loops) {
    config_->generator  = GeneratorType::GNM_DIRECTED;
    config_->n          = n;
    config_->m          = m;
    config_->self_loops = self_loops;
    return GenerateInMemory(*config_, representation_, comm_);
}

Graph KaGen::GenerateUndirectedGNM(const SInt n, const SInt m, const bool self_loops) {
    config_->generator  = GeneratorType::GNM_UNDIRECTED;
    config_->n          = n;
    config_->m          = m;
    config_->self_loops = self_loops;
    return GenerateInMemory(*config_, representation_, comm_);
}

Graph KaGen::GenerateDirectedGNP(const SInt n, const LPFloat p, const bool self_loops) {
    config_->generator  = GeneratorType::GNP_DIRECTED;
    config_->n          = n;
    config_->p          = p;
    config_->self_loops = self_loops;
    return GenerateInMemory(*config_, representation_, comm_);
}

Graph KaGen::GenerateUndirectedGNP(const SInt n, const LPFloat p, const bool self_loops) {
    config_->generator  = GeneratorType::GNP_UNDIRECTED;
    config_->n          = n;
    config_->p          = p;
    config_->self_loops = self_loops;
    return GenerateInMemory(*config_, representation_, comm_);
}

namespace {
Graph GenerateRGG2D_Impl(
    PGeneratorConfig& config, const SInt n, const SInt m, const LPFloat r, const bool coordinates,
    const GraphRepresentation representation, MPI_Comm comm) {
    config.generator   = GeneratorType::RGG_2D;
    config.n           = n;
    config.m           = m;
    config.r           = r;
    config.coordinates = coordinates;
    return GenerateInMemory(config, representation, comm);
}
} // namespace

Graph KaGen::GenerateRGG2D(const SInt n, const LPFloat r, const bool coordinates) {
    return GenerateRGG2D_Impl(*config_, n, 0, r, coordinates, representation_, comm_);
}

Graph KaGen::GenerateRGG2D_NM(const SInt n, const SInt m, const bool coordinates) {
    return GenerateRGG2D_Impl(*config_, n, m, 0.0, coordinates, representation_, comm_);
}

Graph KaGen::GenerateRGG2D_MR(const SInt m, const LPFloat r, const bool coordinates) {
    return GenerateRGG2D_Impl(*config_, 0, m, r, coordinates, representation_, comm_);
}

namespace {
Graph GenerateRGG3D_Impl(
    PGeneratorConfig& config, const SInt n, const SInt m, const LPFloat r, const bool coordinates,
    const GraphRepresentation representation, MPI_Comm comm) {
    config.generator   = GeneratorType::RGG_3D;
    config.m           = m;
    config.n           = n;
    config.r           = r;
    config.coordinates = coordinates;
    return GenerateInMemory(config, representation, comm);
}
} // namespace

Graph KaGen::GenerateRGG3D(const SInt n, const LPFloat r, const bool coordinates) {
    return GenerateRGG3D_Impl(*config_, n, 0, r, coordinates, representation_, comm_);
}

Graph KaGen::GenerateRGG3D_NM(const SInt n, const SInt m, const bool coordinates) {
    return GenerateRGG3D_Impl(*config_, n, m, 0.0, coordinates, representation_, comm_);
}

Graph KaGen::GenerateRGG3D_MR(const SInt m, const LPFloat r, const bool coordinates) {
    return GenerateRGG3D_Impl(*config_, 0, m, r, coordinates, representation_, comm_);
}

#ifdef KAGEN_CGAL_FOUND
namespace {
Graph GenerateRDG2D_Impl(
    PGeneratorConfig& config, const SInt n, const SInt m, const bool periodic, const bool coordinates,
    const GraphRepresentation representation, MPI_Comm comm) {
    config.generator   = GeneratorType::RDG_2D;
    config.n           = n;
    config.m           = m;
    config.periodic    = periodic;
    config.coordinates = coordinates;
    return GenerateInMemory(config, representation, comm);
}
} // namespace

Graph KaGen::GenerateRDG2D(const SInt n, const bool periodic, const bool coordinates) {
    return GenerateRDG2D_Impl(*config_, n, 0, periodic, coordinates, representation_, comm_);
}

Graph KaGen::GenerateRDG2D_M(const SInt m, const bool periodic, const bool coordinates) {
    return GenerateRDG2D_Impl(*config_, 0, m, periodic, coordinates, representation_, comm_);
}

namespace {
Graph GenerateRDG3D_Impl(
    PGeneratorConfig& config, const SInt n, const SInt m, const bool coordinates,
    const GraphRepresentation representation, MPI_Comm comm) {
    config.generator   = GeneratorType::RDG_3D;
    config.n           = n;
    config.m           = m;
    config.coordinates = coordinates;
    return GenerateInMemory(config, representation, comm);
}
} // namespace

Graph KaGen::GenerateRDG3D(const SInt n, const bool coordinates) {
    return GenerateRDG3D_Impl(*config_, n, 0, coordinates, representation_, comm_);
}

Graph KaGen::GenerateRDG3D_M(const SInt m, const bool coordinates) {
    return GenerateRDG3D_Impl(*config_, 0, m, coordinates, representation_, comm_);
}
#else  // KAGEN_CGAL_FOUND
Graph KaGen::GenerateRDG2D(SInt, bool, bool) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}

Graph KaGen::GenerateRDG2D_M(SInt, bool, bool) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}

Graph KaGen::GenerateRDG3D(SInt, bool) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}

Graph KaGen::GenerateRDG3D_M(SInt, bool) {
    throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}
#endif // KAGEN_CGAL_FOUND

namespace {
Graph GenerateBA_Impl(
    PGeneratorConfig& config, const SInt n, const SInt m, const SInt d, const bool directed, const bool self_loops,
    const GraphRepresentation representation, MPI_Comm comm) {
    config.generator  = GeneratorType::BA;
    config.m          = m;
    config.n          = n;
    config.min_degree = d;
    config.self_loops = self_loops;
    config.directed   = directed;
    return GenerateInMemory(config, representation, comm);
}
} // namespace

Graph KaGen::GenerateBA(const SInt n, const SInt d, const bool directed, const bool self_loops) {
    return GenerateBA_Impl(*config_, n, 0, d, directed, self_loops, representation_, comm_);
}

Graph KaGen::GenerateBA_NM(const SInt n, const SInt m, const bool directed, const bool self_loops) {
    return GenerateBA_Impl(*config_, n, m, 0.0, directed, self_loops, representation_, comm_);
}

Graph KaGen::GenerateBA_MD(const SInt m, const SInt d, const bool directed, const bool self_loops) {
    return GenerateBA_Impl(*config_, 0, m, d, directed, self_loops, representation_, comm_);
}

namespace {
Graph GenerateRHG_Impl(
    PGeneratorConfig& config, const LPFloat gamma, const SInt n, const SInt m, const LPFloat d, const bool coordinates,
    const GraphRepresentation representation, MPI_Comm comm) {
    config.generator   = GeneratorType::RHG;
    config.n           = n;
    config.m           = m;
    config.avg_degree  = d;
    config.plexp       = gamma;
    config.coordinates = coordinates;
    return GenerateInMemory(config, representation, comm);
}
} // namespace

Graph KaGen::GenerateRHG(const LPFloat gamma, const SInt n, const LPFloat d, const bool coordinates) {
    return GenerateRHG_Impl(*config_, gamma, n, 0, d, coordinates, representation_, comm_);
}

Graph KaGen::GenerateRHG_NM(const LPFloat gamma, const SInt n, const SInt m, const bool coordinates) {
    return GenerateRHG_Impl(*config_, gamma, n, m, 0.0, coordinates, representation_, comm_);
}

Graph KaGen::GenerateRHG_MD(const LPFloat gamma, const SInt m, const LPFloat d, const bool coordinates) {
    return GenerateRHG_Impl(*config_, gamma, 0, m, d, coordinates, representation_, comm_);
}

namespace {
Graph GenerateGrid2D_Impl(
    PGeneratorConfig& config, const SInt n, const SInt grid_x, const SInt grid_y, const LPFloat p, const SInt m,
    const bool periodic, const bool coordinates, const GraphRepresentation representation, MPI_Comm comm) {
    config.generator   = GeneratorType::GRID_2D;
    config.grid_x      = grid_x;
    config.grid_y      = grid_y;
    config.p           = p;
    config.n           = n;
    config.m           = m;
    config.periodic    = periodic;
    config.coordinates = coordinates;
    return GenerateInMemory(config, representation, comm);
}
} // namespace

Graph KaGen::GenerateGrid2D(
    const SInt grid_x, const SInt grid_y, const LPFloat p, const bool periodic, const bool coordinates) {
    return GenerateGrid2D_Impl(*config_, 0, grid_x, grid_y, p, 0, periodic, coordinates, representation_, comm_);
}

Graph KaGen::GenerateGrid2D_N(const SInt n, const LPFloat p, const bool periodic, const bool coordinates) {
    return GenerateGrid2D_Impl(*config_, n, 0, 0, p, 0, periodic, coordinates, representation_, comm_);
}

Graph KaGen::GenerateGrid2D_NM(const SInt n, const SInt m, const bool periodic, const bool coordinates) {
    return GenerateGrid2D_Impl(*config_, n, 0, 0, 0.0, m, periodic, coordinates, representation_, comm_);
}

namespace {
Graph GenerateGrid3D_Impl(
    PGeneratorConfig& config, const SInt grid_x, const SInt grid_y, const SInt grid_z, const LPFloat p, const SInt m,
    const bool periodic, const bool coordinates, const GraphRepresentation representation, MPI_Comm comm) {
    config.generator   = GeneratorType::GRID_3D;
    config.grid_x      = grid_x;
    config.grid_y      = grid_y;
    config.grid_z      = grid_z;
    config.p           = p;
    config.m           = m;
    config.periodic    = periodic;
    config.coordinates = coordinates;
    return GenerateInMemory(config, representation, comm);
}
} // namespace

Graph KaGen::GenerateGrid3D(
    const SInt grid_x, const SInt grid_y, const SInt grid_z, const LPFloat p, const bool periodic,
    const bool coordinates) {
    return GenerateGrid3D_Impl(*config_, grid_x, grid_y, grid_z, p, 0, periodic, coordinates, representation_, comm_);
}

Graph KaGen::GenerateGrid3D_N(const SInt n, const LPFloat p, const bool periodic, const bool coordinates) {
    const SInt cbrt_n = std::cbrt(n);
    return GenerateGrid3D(cbrt_n, cbrt_n, cbrt_n, p, periodic, coordinates);
}

Graph KaGen::GenerateGrid3D_NM(const SInt n, const SInt m, const bool periodic, const bool coordinates) {
    const SInt cbrt_n = std::cbrt(n);
    return GenerateGrid3D_Impl(*config_, cbrt_n, cbrt_n, cbrt_n, 0.0, m, periodic, coordinates, representation_, comm_);
}

Graph KaGen::GenerateDirectedPath(unsigned long long n, bool permute, bool periodic) {
    config_->generator = GeneratorType::PATH_DIRECTED;
    config_->n         = n;
    config_->permute   = permute;
    config_->periodic  = periodic;
    return GenerateInMemory(*config_, representation_, comm_);
}

Graph KaGen::GenerateKronecker(const SInt n, const SInt m, const bool directed, const bool self_loops) {
    config_->generator  = GeneratorType::KRONECKER;
    config_->n          = n;
    config_->m          = m;
    config_->directed   = directed;
    config_->self_loops = self_loops;
    return GenerateInMemory(*config_, representation_, comm_);
}

Graph KaGen::GenerateRMAT(
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
    return GenerateInMemory(*config_, representation_, comm_);
}

Graph KaGen::ReadFromFile(std::string const& filename, const FileFormat format, const GraphDistribution distribution) {
    config_->generator                = GeneratorType::FILE;
    config_->input_graph.filename     = filename;
    config_->input_graph.format       = format;
    config_->input_graph.distribution = distribution;
    return GenerateInMemory(*config_, representation_, comm_);
}

void KaGen::SetDefaults() {
    config_->quiet = true;
    config_->output_graph.formats.clear();
    // (keep other defaults)
}

//
// Streaming interface
//

[[nodiscard]] SInt StreamedGraph::NumberOfLocalVertices() const {
    return vertex_range.second - vertex_range.first;
}

[[nodiscard]] SInt StreamedGraph::NumberOfLocalEdges() const {
    return primary_edges.size() + secondary_edges.size();
}

sKaGen::sKaGen(const std::string& options, PEID chunks_per_pe, MPI_Comm comm)
    : generator_(std::make_unique<StreamingGenerator>(options, chunks_per_pe, comm)) {}

sKaGen::~sKaGen() = default;

VertexRange sKaGen::EstimateVertexRange(const PEID pe) const {
    return generator_->EstimateVertexRange(pe);
}

void sKaGen::Initialize() {
    generator_->Initialize();
}

[[nodiscard]] StreamedGraph sKaGen::Next() {
    return generator_->Next();
}

[[nodiscard]] bool sKaGen::Continue() {
    return generator_->Continue();
}
} // namespace kagen
