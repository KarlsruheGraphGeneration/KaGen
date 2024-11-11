#include "kagen/generators/generator.h"

#include "kagen/context.h"
#include "kagen/edgeweight_generators/default_generator.h"
#include "kagen/edgeweight_generators/edge_weight_generator.h"
#include "kagen/edgeweight_generators/hashing_based_generator.h"
#include "kagen/edgeweight_generators/uniform_random_generator.h"
#include "kagen/edgeweight_generators/voiding_generator.h"
#include "kagen/kagen.h"
#include "kagen/tools/converter.h"
#include "kagen/tools/random_permutation.h"
#include "kagen/vertexweight_generators/default_generator.h"
#include "kagen/vertexweight_generators/uniform_random_generator.h"
#include "kagen/vertexweight_generators/vertex_weight_generator.h"
#include "kagen/vertexweight_generators/voiding_generator.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>

namespace kagen {
Generator::~Generator() = default;

Generator* Generator::Generate(const GraphRepresentation representation) {
    Reset();
    desired_representation_ = representation;

    switch (desired_representation_) {
        case GraphRepresentation::EDGE_LIST:
            GenerateEdgeList();
            break;

        case GraphRepresentation::CSR:
            GenerateCSR();
            break;
    }

    return this;
}

Generator* Generator::Finalize(MPI_Comm comm) {
    switch (desired_representation_) {
        case GraphRepresentation::EDGE_LIST:
            FinalizeEdgeList(comm);
            break;

        case GraphRepresentation::CSR:
            FinalizeCSR(comm);
            break;
    }

    graph_.representation = desired_representation_;

    return this;
}

std::unique_ptr<kagen::EdgeWeightGenerator>
CreateEdgeWeightGenerator(const EdgeWeightConfig weight_config, MPI_Comm comm, const VertexRange vertex_range) {
    switch (weight_config.generator_type) {
        case EdgeWeightGeneratorType::DEFAULT:
            return std::make_unique<DefaultEdgeWeightGenerator>(weight_config);
        case EdgeWeightGeneratorType::VOIDING:
            return std::make_unique<VoidingEdgeWeightGenerator>(weight_config);
        case EdgeWeightGeneratorType::HASHING_BASED:
            return std::make_unique<HashingBasedEdgeWeightGenerator>(weight_config);
        case EdgeWeightGeneratorType::UNIFORM_RANDOM:
            return std::make_unique<UniformRandomEdgeWeightGenerator>(weight_config, comm, vertex_range);
    }

    throw std::runtime_error("invalid weight generator type");
}

void Generator::GenerateEdgeWeights(EdgeWeightConfig weight_config, MPI_Comm comm) {
    std::unique_ptr<kagen::EdgeWeightGenerator> edge_weight_generator =
        CreateEdgeWeightGenerator(weight_config, comm, graph_.vertex_range);

    switch (desired_representation_) {
        case GraphRepresentation::EDGE_LIST:
            edge_weight_generator->GenerateEdgeWeights(graph_.edges, graph_.edge_weights);
            break;
        case GraphRepresentation::CSR:
            edge_weight_generator->GenerateEdgeWeights(graph_.xadj, graph_.adjncy, graph_.edge_weights);
            break;
    }
}

namespace {
template <typename Permutator>
auto ApplyPermutationAndComputeSendBuffersEdgeList(
    const Graph& graph, const std::vector<VertexRange>& recv_ranges, Permutator&& permute) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Edgelist edges = graph.edges;
    for (auto& [src, dst]: edges) {
        src = permute(src);
        dst = permute(dst);
    }
    std::unordered_map<PEID, std::vector<SInt>>  send_buffers;
    std::unordered_map<PEID, std::vector<SSInt>> edge_weights_send_buffers;

    bool has_edge_weights = graph.NumberOfLocalEdges() == 0 || !graph.edge_weights.empty();

    for (std::size_t i = 0; i < edges.size(); ++i) {
        const auto& [src, dst]               = edges[i];
        const PEID          target_pe        = FindPEInRange(src, recv_ranges);
        std::vector<SInt>&  send_buf         = send_buffers[target_pe];
        std::vector<SSInt>& weights_send_buf = edge_weights_send_buffers[target_pe];
        send_buf.push_back(src);
        send_buf.push_back(dst);
        // [Permuted_Src_Id, Degree, EdgeWeights]
        if (has_edge_weights) {
            weights_send_buf.push_back(graph.edge_weights[i]);
        }
    }
    return std::make_tuple(std::move(send_buffers), std::move(edge_weights_send_buffers));
}
template <typename Permutator>
auto ApplyPermutationAndComputeSendBuffersCSR(
    const Graph& graph, const std::vector<VertexRange>& recv_ranges, Permutator&& permute) {
    AdjncyArray permuted_adjncy = graph.adjncy;

    for (auto& edge: permuted_adjncy) {
        edge = permute(edge);
    }

    std::unordered_map<PEID, std::vector<SInt>>  send_buffers;
    std::unordered_map<PEID, std::vector<SSInt>> edge_weights_send_buffers;

    bool has_edge_weights = graph.NumberOfLocalEdges() == 0 || !graph.edge_weights.empty();

    for (std::size_t i = 0; i + 1 < graph.xadj.size(); ++i) {
        const SInt          degree             = graph.xadj[i + 1] - graph.xadj[i];
        const SInt          global_id          = graph.vertex_range.first + i;
        const SInt          permuted_global_id = permute(global_id);
        const PEID          target_pe          = FindPEInRange(permuted_global_id, recv_ranges);
        std::vector<SInt>&  send_buf           = send_buffers[target_pe];
        std::vector<SSInt>& weights_send_buf   = edge_weights_send_buffers[target_pe];
        auto                edge_begin_offset  = graph.xadj[i];
        auto                edge_end_offset    = graph.xadj[i + 1];
        // [Permuted_Src_Id, Degree, [Permuted_Dst_Ids]] with #Permuted_Dst_Ids = Degree
        send_buf.push_back(permuted_global_id);
        send_buf.push_back(degree);
        send_buf.insert(
            send_buf.end(), permuted_adjncy.begin() + edge_begin_offset, permuted_adjncy.begin() + edge_end_offset);
        // [Permuted_Src_Id, Degree, EdgeWeights]
        if (has_edge_weights) {
            weights_send_buf.push_back(permuted_global_id);
            weights_send_buf.push_back(degree);
            weights_send_buf.insert(
                weights_send_buf.end(), graph.edge_weights.begin() + edge_begin_offset,
                graph.edge_weights.begin() + edge_end_offset);
        }
    }
    return std::make_tuple(std::move(send_buffers), std::move(edge_weights_send_buffers));
}

template <typename Permutator>
auto ApplyPermutationAndComputeSendBuffers(
    const Graph& graph, const std::vector<VertexRange>& recv_ranges, Permutator&& permutator) {
    switch (graph.representation) {
        case GraphRepresentation::EDGE_LIST:
            return ApplyPermutationAndComputeSendBuffersEdgeList(
                graph, recv_ranges, std::forward<Permutator>(permutator));
        case GraphRepresentation::CSR:
            return ApplyPermutationAndComputeSendBuffersCSR(graph, recv_ranges, std::forward<Permutator>(permutator));
        default:
            throw std::runtime_error("Unexpected graph representation type.");
    }
}

inline auto ConstructPermutedGraphCSR(
    VertexRange recv_range, const std::vector<SInt>& recv_edges, const std::vector<SSInt>& recv_edge_weights) {
    std::size_t       num_local_vertices = recv_range.second - recv_range.first;
    std::vector<SInt> degree_count(num_local_vertices, 0);

    // scan received data for degrees
    for (std::size_t cur_pos = 0; cur_pos < recv_edges.size();) {
        const SInt src_id                       = recv_edges[cur_pos];
        const SInt degree                       = recv_edges[cur_pos + 1];
        degree_count[src_id - recv_range.first] = degree;
        // skip edges
        cur_pos += 1 + degree + 1;
    }
    XadjArray xadj(num_local_vertices + 1, 0);
    // compute xadj array for received graph
    std::exclusive_scan(degree_count.begin(), degree_count.end(), xadj.begin(), SInt{0});
    xadj.back() = degree_count.back() + xadj[num_local_vertices - 1];

    // compute adjncy for received graph
    const std::size_t num_local_edges = xadj.back();
    XadjArray         adjncy(num_local_edges);
    for (std::size_t cur_pos = 0; cur_pos < recv_edges.size();) {
        const SInt global_src_id = recv_edges[cur_pos];
        const SInt degree        = recv_edges[cur_pos + 1];
        const SInt local_src_id  = global_src_id - recv_range.first;
        std::copy_n(recv_edges.begin() + cur_pos + 2, degree, adjncy.begin() + xadj[local_src_id]);
        // forward to next received src vertex
        cur_pos += 1 + degree + 1;
    }
    // compute edge weights for received graph
    EdgeWeights edge_weights(num_local_edges);
    if (!recv_edge_weights.empty()) {
        for (std::size_t cur_pos = 0; cur_pos < recv_edge_weights.size();) {
            const SInt global_src_id = static_cast<SInt>(recv_edge_weights[cur_pos]);
            const SInt degree        = static_cast<SInt>(recv_edge_weights[cur_pos + 1]);
            const SInt local_src_id  = global_src_id - recv_range.first;
            std::copy_n(recv_edge_weights.begin() + cur_pos + 2, degree, edge_weights.begin() + xadj[local_src_id]);
            // forward to next received src vertex
            cur_pos += 1 + degree + 1;
        }
    }

    // TODO handle vertex weights
    return std::make_tuple(std::move(xadj), std::move(adjncy), std::move(edge_weights));
}
inline auto ConstructPermutedGraphEdgeList(
    VertexRange recv_range, const std::vector<SInt>& recv_edges, const std::vector<SSInt>& recv_edge_weights) {
    std::size_t       num_local_vertices = recv_range.second - recv_range.first;
    std::vector<SInt> degree(num_local_vertices, 0);
    int               rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for (std::size_t i = 0; i < recv_edges.size(); i += 2) {
        const auto src = recv_edges[i];
        ++degree[src - recv_range.first];
    }
    std::vector<SInt> write_idx(num_local_vertices);
    std::exclusive_scan(degree.begin(), degree.end(), write_idx.begin(), SInt{0});
    Edgelist    edgelist(recv_edges.size() / 2); // edges are sent flat - not as pairs
    EdgeWeights edge_weights(recv_edge_weights.size());
    for (std::size_t i = 0; i < recv_edges.size(); i += 2) {
        const auto src          = recv_edges[i];
        const auto dst          = recv_edges[i + 1];
        const auto local_src_id = src - recv_range.first;
        const auto idx          = write_idx[local_src_id];
        edgelist[idx]           = std::make_pair(src, dst);
        if (!edge_weights.empty()) {
            edge_weights[idx] = recv_edge_weights[i / 2];
        }
        ++write_idx[local_src_id];
    }
    return std::make_tuple(std::move(edgelist), std::move(edge_weights));
}
} // namespace

void Generator::PermuteVertices(const PGeneratorConfig& config, MPI_Comm comm) {
    if (!graph_.vertex_weights.empty())
        throw std::runtime_error(
            "Graph is vertex weight but this is not yet supported by the vertex permutation routine!");
    int size = -1;
    int rank = -1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    auto permutator = random_permutation::FeistelPseudoRandomPermutation::buildPermutation(config.n - 1, 0);
    auto permute    = [&permutator](SInt v) {
        return permutator.f(v);
    };

    // all PE get n / size vertices
    // the first n modulo size PEs obtain one additional vertices.
    const SInt vertices_per_pe             = config.n / size;
    const PEID num_pe_with_additional_node = config.n % size;
    const bool has_pe_additional_node      = rank < num_pe_with_additional_node;
    const SInt begin_vertices              = std::min(num_pe_with_additional_node, rank) + rank * vertices_per_pe;
    const SInt end_vertices                = begin_vertices + vertices_per_pe + has_pe_additional_node;

    VertexRange              recv_range{begin_vertices, end_vertices};
    std::vector<VertexRange> recv_ranges = AllgatherVertexRange(recv_range, comm);

    auto [send_buffers, edge_weight_send_buffers] = ApplyPermutationAndComputeSendBuffers(graph_, recv_ranges, permute);
    auto recv_edges        = ExchangeMessageBuffers(std::move(send_buffers), KAGEN_MPI_SINT, comm);
    auto recv_edge_weights = ExchangeMessageBuffers(std::move(edge_weight_send_buffers), KAGEN_MPI_SSINT, comm);

    switch (desired_representation_) {
        case GraphRepresentation::EDGE_LIST: {
            auto [permuted_edgelist, permuted_edge_weights] =
                ConstructPermutedGraphEdgeList(recv_range, recv_edges, recv_edge_weights);
            graph_.edges        = std::move(permuted_edgelist);
            graph_.edge_weights = std::move(permuted_edge_weights);
            break;
        }
        case GraphRepresentation::CSR: {
            auto [permuted_xadj, permuted_adjncy, permuted_edge_weights] =
                ConstructPermutedGraphCSR(recv_range, recv_edges, recv_edge_weights);

            graph_.xadj         = std::move(permuted_xadj);
            graph_.adjncy       = std::move(permuted_adjncy);
            graph_.edge_weights = std::move(permuted_edge_weights);
            break;
        }
    }
    SetVertexRange(recv_range);
}

std::unique_ptr<kagen::VertexWeightGenerator>
CreateVertexWeightGenerator(const VertexWeightConfig weight_config, MPI_Comm comm) {
    switch (weight_config.generator_type) {
        case VertexWeightGeneratorType::DEFAULT:
            return std::make_unique<DefaultVertexWeightGenerator>(weight_config);
        case VertexWeightGeneratorType::VOIDING:
            return std::make_unique<VoidingVertexWeightGenerator>(weight_config);
        case VertexWeightGeneratorType::UNIFORM_RANDOM:
            return std::make_unique<UniformRandomVertexWeightGenerator>(weight_config, comm);
    }

    throw std::runtime_error("invalid weight generator type");
}

void Generator::GenerateVertexWeights(VertexWeightConfig weight_config, MPI_Comm comm) {
    std::unique_ptr<kagen::VertexWeightGenerator> vertex_weight_generator =
        CreateVertexWeightGenerator(weight_config, comm);

    switch (desired_representation_) {
        case GraphRepresentation::EDGE_LIST:
            vertex_weight_generator->GenerateVertexWeights(graph_.vertex_range, graph_.edges, graph_.vertex_weights);
            break;
        case GraphRepresentation::CSR:
            vertex_weight_generator->GenerateVertexWeights(
                graph_.vertex_range, graph_.xadj, graph_.adjncy, graph_.vertex_weights);

            break;
    }
}

void Generator::FinalizeEdgeList(MPI_Comm) {}

void Generator::FinalizeCSR(MPI_Comm) {}

void CSROnlyGenerator::GenerateEdgeList() {
    GenerateCSR();
}

void CSROnlyGenerator::FinalizeEdgeList(MPI_Comm comm) {
    if (graph_.xadj.empty()) {
        return;
    }

    // Otherwise, we have generated the graph in CSR representation, but
    // actually want edge list representation -> transform graph
    FinalizeCSR(comm);
    graph_.edges = BuildEdgeListFromCSR(graph_.vertex_range, graph_.xadj, graph_.adjncy);
    {
        XadjArray tmp;
        std::swap(graph_.xadj, tmp);
    }
    {
        AdjncyArray tmp;
        std::swap(graph_.adjncy, tmp);
    }
}

void EdgeListOnlyGenerator::GenerateCSR() {
    GenerateEdgeList();
}

void EdgeListOnlyGenerator::FinalizeCSR(MPI_Comm comm) {
    if (!graph_.xadj.empty()) {
        return;
    }

    // Otherwise, we have generated the graph in edge list representation, but
    // actually want CSR format --> transform graph
    FinalizeEdgeList(comm);
    std::tie(graph_.xadj, graph_.adjncy) = BuildCSRFromEdgeList(graph_.vertex_range, graph_.edges, graph_.edge_weights);
    {
        Edgelist tmp;
        std::swap(graph_.edges, tmp);
    }
}

SInt Generator::GetNumberOfEdges() const {
    return std::max(graph_.adjncy.size(), graph_.edges.size());
}

Graph Generator::Take() {
    return std::move(graph_);
}

void Generator::SetVertexRange(const VertexRange vertex_range) {
    graph_.vertex_range = vertex_range;
}

void Generator::FilterDuplicateEdges() {
    std::sort(graph_.edges.begin(), graph_.edges.end());
    auto it = std::unique(graph_.edges.begin(), graph_.edges.end());
    graph_.edges.erase(it, graph_.edges.end());
}

void Generator::Reset() {
    graph_.Clear();
}

GeneratorFactory::~GeneratorFactory() = default;

PGeneratorConfig
GeneratorFactory::NormalizeParameters(PGeneratorConfig config, PEID, const PEID size, const bool output) const {
    if (config.k == 0) {
        config.k = static_cast<SInt>(size);
        if (output) {
            std::cout << "Setting number of chunks to " << config.k << std::endl;
        }
    }
    return config;
}

namespace {
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
} // namespace

void GeneratorFactory::EnsureSquarePowerOfTwoChunkSize(
    PGeneratorConfig& config, const PEID size, const bool output) const {
    if (config.k == 0) {
        if (IsSquare(size) && IsPowerOfTwo(size)) {
            config.k = static_cast<SInt>(size);
        } else {
            const SInt l = std::ceil(std::log2(size));
            config.k     = 1 << l;
            if (!IsSquare(config.k)) {
                config.k *= 2;
            }
            while (std::ceil(1.0 * config.k / size) > (1.0 + config.max_vertex_imbalance) * config.k / size) {
                config.k <<= 2;
            }
        }
        if (output) {
            std::cout << "Setting number of chunks to " << config.k << std::endl;
        }
    } else if (config.k < static_cast<SInt>(size) || !IsSquare(config.k) || !IsPowerOfTwo(config.k)) {
        throw ConfigurationError("number of chunks must be square power of two and larger than number of PEs");
    }
}

void GeneratorFactory::EnsureCubicPowerOfTwoChunkSize(
    PGeneratorConfig& config, const PEID size, const bool output) const {
    if (config.k == 0) {
        if (IsCubic(size) && IsPowerOfTwo(size)) {
            config.k = static_cast<SInt>(size);
        } else {
            const SInt l = std::ceil(std::log2(size));
            config.k     = 1 << l;
            if (!IsCubic(config.k)) {
                config.k *= 2;
            }
            if (!IsCubic(config.k)) {
                config.k *= 2;
            }

            while (std::ceil(1.0 * config.k / size) > (1.0 + config.max_vertex_imbalance) * config.k / size) {
                config.k <<= 3;
            }
        }
        if (output) {
            std::cout << "Setting number of chunks to " << config.k << std::endl;
        }
    } else if (config.k < static_cast<SInt>(size) || !IsCubic(config.k)) {
        throw ConfigurationError("number of chunks must be cubic and larger than the number of PEs");
    }
}

void GeneratorFactory::EnsureOneChunkPerPE(PGeneratorConfig& config, const PEID size) const {
    if (config.k != static_cast<SInt>(size)) {
        throw ConfigurationError("number of chunks must match the number of PEs");
    }
}
} // namespace kagen
