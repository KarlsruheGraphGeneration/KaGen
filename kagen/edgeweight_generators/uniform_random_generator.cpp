#include "kagen/edgeweight_generators/uniform_random_generator.h"

#include "kagen/context.h"
#include "kagen/kagen.h"
#include "kagen/tools/converter.h"
#include "kagen/tools/utils.h"

#include "xxhash.h"
#include <unordered_map>

namespace kagen {

namespace {

struct EdgeData {
    SInt  u;
    SInt  v;
    SSInt randomness;
    SSInt weight;
    EdgeData() = default;
    EdgeData(SInt u_param, SInt v_param, SSInt randomness_param, SSInt weight_param)
        : u{u_param},
          v{v_param},
          randomness{randomness_param},
          weight{weight_param} {}
};

// Stores edge (u,v) together with an associated random integer and an edge weight.
// Since an edge weight for each (directed) edge (u,v) and (v,u) is chosen uniformly at random, one has to agree on a
// common weight for the undirected edge {u,v}. This is done choosing the weight with the smaller random integer to
// avoid biases. The random integer is also used to break ties in case KaGen generates duplicate edges.
class EdgeWeightStorage {
    using Edge           = typename Edgelist::value_type;
    using RandInt_Weight = std::pair<SSInt, SSInt>;

public:
    void InsertOrReplace(const Edge& key, const RandInt_Weight& value) {
        auto it = edge_to_weightdata.find(key);
        if (it != edge_to_weightdata.end() && value < it->second) {
            // if edge is already present (due to duplicate edges) store edge weight with smaller (rand_int, weight)
            // pair.
            it->second = value;
        } else {
            edge_to_weightdata.emplace(key, value);
        }
    }

    // Use weight with smaller associated random integer.
    SSInt AgreeOnEdgeWeight(const Edge& edge) {
        const auto& reversed_edge = std::make_pair(edge.second, edge.first);
        auto        edge_it       = edge_to_weightdata.find(edge);
        if (edge_it == edge_to_weightdata.end()) {
            throw std::runtime_error("edge should have been inserted into map!");
        }
        const auto& edge_value        = edge_it->second;
        SSInt       weight            = edge_value.second; // use this value as default
        auto        reveresed_edge_it = edge_to_weightdata.find(reversed_edge);
        if (reveresed_edge_it != edge_to_weightdata.end() && reveresed_edge_it->second < edge_value) {
            const auto& reversed_edge_value = reveresed_edge_it->second;
            weight                          = reversed_edge_value.second;
        }
        return weight;
    }

private:
    struct EdgeHasher {
        size_t operator()(const Edge& edge) const noexcept {
            return edge.first ^ (edge.second << 1);
        }
    };
    std::unordered_map<Edge, std::pair<SSInt, SSInt>, EdgeHasher> edge_to_weightdata;
};
} // namespace

UniformRandomEdgeWeightGenerator::UniformRandomEdgeWeightGenerator(
    EdgeWeightConfig config, MPI_Comm comm, VertexRange vertex_range)
    : config_(config),
      comm_{comm},
      vertex_range_{vertex_range} {
    if (config_.weight_range_begin >= config_.weight_range_end) {
        throw std::runtime_error("Weight causes undefined behavior, need weight_range_begin > weight_range_end.");
    }
}

EdgeWeights UniformRandomEdgeWeightGenerator::GenerateEdgeWeights(const XadjArray& xadj, const AdjncyArray& adjncy) {
    const auto edge_list = BuildEdgeListFromCSR(vertex_range_, xadj, adjncy);
    return GenerateEdgeWeights(edge_list);
}

EdgeWeights UniformRandomEdgeWeightGenerator::GenerateEdgeWeights(const Edgelist& edgelist) {
    PEID rank;
    MPI_Comm_rank(comm_, &rank);
    std::mt19937                         gen((rank + 42) * 3);
    std::uniform_int_distribution<SSInt> weight_dist(config_.weight_range_begin, config_.weight_range_end - 1);

    const auto ranges                  = AllgatherVertexRange(vertex_range_, comm_);
    const auto& [local_from, local_to] = ranges[rank];
    std::unordered_map<PEID, std::vector<EdgeData>> message_buffers;
    EdgeWeightStorage                               edge_weight_storage;

    // store (rand_int, weight) for each edge and exchange this data for cut edges
    for (size_t i = 0; i < edgelist.size(); ++i) {
        const auto& edge         = edgelist[i];
        const auto& [tail, head] = edge;
        const SSInt randomness   = weight_dist(gen);
        const SSInt weight       = weight_dist(gen);
        // Each PE gets the edges that we have to that PE
        if ((head < local_from || head >= local_to)) {
            const SInt pe = static_cast<SInt>(FindPEInRange(head, ranges));
            message_buffers[pe].emplace_back(tail, head, randomness, weight);
        }
        auto value = std::make_pair(randomness, weight);
        edge_weight_storage.InsertOrReplace(edge, value);
    }

    {
        // exchange messages
        MPI_Datatype edgedata_mpi_type;
        MPI_Type_contiguous(sizeof(EdgeData), MPI_BYTE, &edgedata_mpi_type);
        MPI_Type_commit(&edgedata_mpi_type);
        auto recv_buf = ExchangeMessageBuffers(message_buffers, edgedata_mpi_type, comm_);
        MPI_Type_free(&edgedata_mpi_type);

        // add received cut edges into edge_weight storage
        for (const auto& [u, v, randomness, weight]: recv_buf) {
            const auto key   = std::make_pair(u, v);
            const auto value = std::make_pair(randomness, weight);
            edge_weight_storage.InsertOrReplace(key, value);
        }
    }

    // agree on which weight to choose
    EdgeWeights weights(edgelist.size());
    for (size_t i = 0; i < edgelist.size(); ++i) {
        const auto& edge = edgelist[i];
        weights[i]       = edge_weight_storage.AgreeOnEdgeWeight(edge);
    }
    return weights;
}
} // namespace kagen
