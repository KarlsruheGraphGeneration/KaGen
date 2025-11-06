#include "kagen/edgeweight_generators/uniform_random_generator.h"

#include "kagen/context.h"
#include "kagen/edge_range.h"
#include "kagen/kagen.h"
#include "kagen/tools/converter.h"
#include "kagen/tools/utils.h"

#include <cassert>
#include <random>
#include <unordered_map>
#include <utility>

namespace kagen {

namespace {

template <typename T>
void dump(T&& t) {
    t.clear();
    t.shrink_to_fit();
    T temp{std::move(t)};
}

auto get_random_generator(int rank) {
    std::mt19937 gen((rank + 42) * 3);
    return gen;
}

struct EdgeHasher {
    using Edge = std::pair<SInt, SInt>;
    inline void hash_combine(std::size_t& seed, std::size_t v) const noexcept {
        seed ^= v + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    }
    size_t operator()(const Edge& edge) const noexcept {
        std::size_t h1   = std::hash<std::uint64_t>{}(static_cast<uint64_t>(edge.first));
        std::size_t h2   = std::hash<std::uint64_t>{}(static_cast<uint64_t>(edge.second));
        std::size_t seed = h1;
        hash_combine(seed, h2);
	return seed;
    }
};

std::pair<SInt, SInt> get_canonical_order(std::pair<SInt, SInt> const& edge) {
    auto u = std::get<0>(edge);
    auto v = std::get<1>(edge);
    if (u < v) {
        return std::make_pair(u, v);
    } else {
        return std::make_pair(v, u);
    }
}

bool is_canonically_ordered(SInt u, SInt v) {
    return u < v;
}

template <typename Edge>
bool is_canonically_ordered(Edge const& edge) {
    return is_canonically_ordered(std::get<0>(edge), std::get<1>(edge));
}

struct EdgeData {
    SInt  u;
    SInt  v;
    SSInt weight;

    EdgeData() = default;

    EdgeData(SInt u_param, SInt v_param, SSInt weight_param) : u{u_param}, v{v_param}, weight{weight_param} {}

    EdgeData(std::pair<SInt, SInt> const& edge, SSInt weight)
        : EdgeData(std::get<0>(edge), std::get<1>(edge), weight) {}
};

struct SrcDstOrder {
    bool operator()(const EdgeData& lhs, const EdgeData& rhs) const {
        return std::tie(lhs.u, lhs.v) < std::tie(rhs.u, rhs.v);
    }
};

struct SrcDstEqual {
    bool operator()(const EdgeData& lhs, const EdgeData& rhs) const {
        return std::tie(lhs.u, lhs.v) == std::tie(rhs.u, rhs.v);
    }
};

void GenerateEdgeWeightsImpl(
    EdgeRange edge_range, EdgeWeights& weights, const EdgeWeightConfig& config, VertexRange vertex_range,
    MPI_Comm comm_) {
    PEID rank;
    MPI_Comm_rank(comm_, &rank);
    std::mt19937                         gen = get_random_generator(rank);
    std::uniform_int_distribution<SSInt> weight_dist(config.weight_range_begin, config.weight_range_end - 1);

    const auto ranges                  = AllgatherVertexRange(vertex_range, comm_);
    const auto& [local_from, local_to] = ranges[rank];
    using Edge                         = typename EdgeRange::Edge;

    std::unordered_map<PEID, std::vector<EdgeData>> message_buffers;
    std::unordered_map<Edge, SInt, EdgeHasher>      edge_to_weight_map;

    for (const auto& edge: edge_range) {
        if (is_canonically_ordered(edge)) {
            auto it = edge_to_weight_map.find(edge);
            if (it == edge_to_weight_map.end()) {
                const SSInt weight = weight_dist(gen);
                const SInt  head   = std::get<1>(edge);
                edge_to_weight_map.emplace(edge, weight);
                if ((head < local_from || head >= local_to)) {
                    const SInt pe = static_cast<SInt>(FindPEInRange(head, ranges));
                    message_buffers[pe].emplace_back(edge, weight);
                }
            }
        }
    }

    MPI_Datatype edgedata_mpi_type;
    MPI_Type_contiguous(sizeof(EdgeData), MPI_BYTE, &edgedata_mpi_type);
    MPI_Type_commit(&edgedata_mpi_type);
    auto recv_buf = ExchangeMessageBuffers(std::move(message_buffers), edgedata_mpi_type, comm_);
    MPI_Type_free(&edgedata_mpi_type);

    // add received cut edges into edge_weight storage
    for (const auto& [u, v, weight]: recv_buf) {
        const auto key = std::make_pair(u, v);
        assert(is_canonically_ordered(u, v));
        edge_to_weight_map.emplace(key, weight);
    }
    dump(std::move(recv_buf));
    weights.resize(edge_range.size());
    for (auto it = edge_range.begin(); it != edge_range.end(); ++it) {
        auto find_result_it = edge_to_weight_map.find(get_canonical_order(*it));
        assert(find_result_it != edge_to_weight_map.end());
        weights[it.edge_index()] = find_result_it->second;
    }
}

} // namespace

UniformRandomEdgeWeightGenerator::UniformRandomEdgeWeightGenerator(
    EdgeWeightConfig config, MPI_Comm comm, VertexRange vertex_range)
    : config_(config),
      comm_(comm),
      vertex_range_(vertex_range) {
    if (config_.weight_range_begin >= config_.weight_range_end) {
        throw std::runtime_error("Weight causes undefined behavior, need weight_range_begin < weight_range_end.");
    }
}

void UniformRandomEdgeWeightGenerator::GenerateEdgeWeights(
    const XadjArray& xadj, const AdjncyArray& adjncy, EdgeWeights& weights) {
    EdgeRange edge_range(xadj, adjncy, vertex_range_);
    GenerateEdgeWeightsImpl(edge_range, weights, config_, vertex_range_, comm_);
}

void UniformRandomEdgeWeightGenerator::GenerateEdgeWeights(const Edgelist& edgelist, EdgeWeights& weights) {
    EdgeRange edge_range(edgelist);
    GenerateEdgeWeightsImpl(edge_range, weights, config_, vertex_range_, comm_);
}
} // namespace kagen
