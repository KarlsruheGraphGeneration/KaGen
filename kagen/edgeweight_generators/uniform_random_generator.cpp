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
    size_t operator()(const Edge& edge) const noexcept {
        return edge.first ^ (edge.second << 1);
    }
};

std::pair<SInt, SInt> get_canonical_order(SInt u, SInt v) {
    if (u < v) {
        return std::make_pair(u, v);
    } else {
        return std::make_pair(v, u);
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
};

bool is_canonically_ordered(SInt u, SInt v) {
    return u < v;
};

template <typename Edge>
bool is_canonically_ordered(Edge const& edge) {
    return is_canonically_ordered(std::get<0>(edge), std::get<1>(edge));
};

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

// void UniformRandomEdgeWeightGenerator::GenerateEdgeWeights(
//     // TODO add more native implementation
//     const XadjArray& xadj, const AdjncyArray& adjncy, EdgeWeights& weights) {
//     int size;
//     MPI_Comm_size(comm_, &size);
//     // if (size == 1) {
//     //     GenerateEdgeWeightsForSingleRanks(xadj, adjncy, weights, config_);
//     // } else {
//     const auto edge_list = BuildEdgeListFromCSR(vertex_range_, xadj, adjncy);
//     GenerateEdgeWeights(edge_list, weights);
//     //}
// }

void UniformRandomEdgeWeightGenerator::GenerateEdgeWeights(
    const XadjArray& xadj, const AdjncyArray& adjncy, EdgeWeights& weights) {
    PEID rank;
    MPI_Comm_rank(comm_, &rank);
    std::mt19937                         gen = get_random_generator(rank);
    std::uniform_int_distribution<SSInt> weight_dist(config_.weight_range_begin, config_.weight_range_end - 1);

    const auto ranges                  = AllgatherVertexRange(vertex_range_, comm_);
    const auto& [local_from, local_to] = ranges[rank];
    const auto num_local_edges         = adjncy.size();
    const auto num_local_vertices      = vertex_range_.second - vertex_range_.first;
    using Edge                         = Edgelist::value_type;

    std::unordered_map<PEID, std::vector<EdgeData>> message_buffers;
    std::unordered_map<Edge, SInt, EdgeHasher>      edge_to_weight_map;

    EdgeRange range(xadj, adjncy, vertex_range_);
    for (const auto& edge: range) {
        if (is_canonically_ordered(edge)) {
            // auto edge = std::make_pair(u_global, v_global);
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
    // for (std::size_t u_local = 0; u_local < num_local_vertices; ++u_local) {
    //     const SInt        u_global = static_cast<SInt>(local_from + static_cast<SInt>(u_local));
    //     const std::size_t start    = static_cast<std::size_t>(xadj[u_local]);
    //     const std::size_t end      = static_cast<std::size_t>(xadj[u_local + 1]);

    //    for (std::size_t offset = start; offset < end; ++offset) {
    //        const SInt v_global = static_cast<SInt>(adjncy[offset]);
    //        if (is_canonically_ordered(u_global, v_global)) {
    //            auto edge = std::make_pair(u_global, v_global);
    //            auto it   = edge_to_weight_map.find(edge);
    //            if (it == edge_to_weight_map.end()) {
    //                const SSInt weight = weight_dist(gen);
    //                const SInt  head   = std::get<1>(edge);
    //                edge_to_weight_map.emplace(edge, weight);
    //                if ((head < local_from || head >= local_to)) {
    //                    const SInt pe = static_cast<SInt>(FindPEInRange(head, ranges));
    //                    message_buffers[pe].emplace_back(edge, weight);
    //                }
    //            }
    //        }
    //    }
    //}
    // exchange messages
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
    // weights.resize(num_local_edges);
    // for (std::size_t u_local = 0; u_local < num_local_vertices; ++u_local) {
    //     const SInt        u_global = static_cast<SInt>(local_from + static_cast<SInt>(u_local));
    //     const std::size_t start    = static_cast<std::size_t>(xadj[u_local]);
    //     const std::size_t end      = static_cast<std::size_t>(xadj[u_local + 1]);

    //    for (std::size_t offset = start; offset < end; ++offset) {
    //        const SInt v_global = static_cast<SInt>(adjncy[offset]);
    //        auto       it       = edge_to_weight_map.find(get_canonical_order(u_global, v_global));

    //        assert(it != edge_to_weight_map.end());
    //        weights[offset] = it->second;
    //    }
    //}
    weights.resize(num_local_edges);
    for (std::size_t u_local = 0; u_local < num_local_vertices; ++u_local) {
        const SInt        u_global = static_cast<SInt>(local_from + static_cast<SInt>(u_local));
        const std::size_t start    = static_cast<std::size_t>(xadj[u_local]);
        const std::size_t end      = static_cast<std::size_t>(xadj[u_local + 1]);

        for (std::size_t offset = start; offset < end; ++offset) {
            const SInt v_global = static_cast<SInt>(adjncy[offset]);
            auto       it       = edge_to_weight_map.find(get_canonical_order(u_global, v_global));

            assert(it != edge_to_weight_map.end());
            weights[offset] = it->second;
        }
    }
}

void UniformRandomEdgeWeightGenerator::GenerateEdgeWeights(const Edgelist& edgelist, EdgeWeights& weights) {
    using Edge = std::pair<SInt, SInt>;
    static_assert(std::is_same_v<typename Edgelist::value_type, Edge>);

    PEID rank;
    MPI_Comm_rank(comm_, &rank);
    std::mt19937                         gen = get_random_generator(rank);
    std::uniform_int_distribution<SSInt> weight_dist(config_.weight_range_begin, config_.weight_range_end - 1);

    const auto ranges                  = AllgatherVertexRange(vertex_range_, comm_);
    const auto& [local_from, local_to] = ranges[rank];
    std::unordered_map<PEID, std::vector<EdgeData>> message_buffers;
    std::unordered_map<Edge, SInt, EdgeHasher>      edge_to_weight_map;

    for (size_t i = 0; i < edgelist.size(); ++i) {
        const auto& edge = edgelist[i];
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

    // exchange messages
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

    // look up the actual weight
    weights.resize(edgelist.size());
    for (size_t i = 0; i < edgelist.size(); ++i) {
        const auto& edge = edgelist[i];
        auto        it   = edge_to_weight_map.find(get_canonical_order(edge));
        assert(it != edge_to_weight_map.end());
        weights[i] = it->second;
    }
}
} // namespace kagen
