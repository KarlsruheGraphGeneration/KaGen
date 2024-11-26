
#include "kagen/kagen.h"

namespace kagen::testing {
using SrcDstEdgeWeight = std::tuple<SInt, SInt, SSInt>;
inline std::vector<SrcDstEdgeWeight> ConvertToWeightedEdgelist(const Graph& graph) {
    std::vector<SrcDstEdgeWeight> weighted_edges;
    bool                          is_edge_weighted = !graph.edge_weights.empty();
    switch (graph.representation) {
        case kagen::GraphRepresentation::EDGE_LIST: {
            for (std::size_t i = 0; i < graph.edges.size(); ++i) {
                const auto& [src, dst]  = graph.edges[i];
                const SSInt edge_weight = is_edge_weighted ? graph.edge_weights[i] : 0;
                weighted_edges.emplace_back(src, dst, edge_weight);
            }
        }
        case kagen::GraphRepresentation::CSR: {
            const SInt  global_v_offset = graph.vertex_range.first;
            std::size_t j               = 0;
            for (std::size_t i = 0; i + 1 < graph.xadj.size(); ++i) {
                const SInt src = global_v_offset + i;
                for (; j < graph.xadj[i + 1]; ++j) {
                    const SInt  dst         = graph.adjncy[j];
                    const SSInt edge_weight = is_edge_weighted ? graph.edge_weights[j] : 0;
                    weighted_edges.emplace_back(src, dst, edge_weight);
                }
            }
        };
    }
    return weighted_edges;
}
}
