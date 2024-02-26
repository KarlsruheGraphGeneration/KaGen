#pragma once

#include "kagen/kagen.h"

#include <algorithm>
#include <numeric>

namespace kagen {
inline std::pair<std::vector<SInt>, std::vector<SInt>>
BuildCSRFromEdgeList(VertexRange vertex_range, Edgelist& edges, EdgeWeights& edge_weights) {
    const SInt num_local_nodes = vertex_range.second - vertex_range.first;
    const SInt num_local_edges = edges.size();

    // Edges must be sorted by from node
    auto cmp_from = [](const auto& lhs, const auto& rhs) {
        return std::get<0>(lhs) < std::get<0>(rhs);
    };

    if (!std::is_sorted(edges.begin(), edges.end(), cmp_from)) {
        // If we have edge weights, sort them the same way as the edges
        // This not very efficient; ideally, we should probably implement some kind of zip iterator to sort edges
        // and edge weights without extra allocation / expensive permutation step (@todo)
        if (!edge_weights.empty()) {
            std::vector<EdgeWeights::value_type> indices(num_local_edges);
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&](const auto& lhs, const auto& rhs) {
                return cmp_from(edges[lhs], edges[rhs]);
            });
            for (std::size_t e = 0; e < num_local_edges; ++e) {
                indices[e] = edge_weights[indices[e]];
            }
            std::swap(edge_weights, indices);
        }

        std::sort(edges.begin(), edges.end(), cmp_from);
    }

    std::vector<SInt> xadj(num_local_nodes + 1);
    std::vector<SInt> adjncy(num_local_edges);

    SInt cur_vertex = 0;
    SInt cur_edge   = 0;

    for (const auto& [from, to]: edges) {
        while (from - vertex_range.first > cur_vertex) {
            xadj[++cur_vertex] = cur_edge;
        }
        adjncy[cur_edge++] = to;
    }
    while (cur_vertex < num_local_nodes) {
        xadj[++cur_vertex] = cur_edge;
    }

    return {std::move(xadj), std::move(adjncy)};
}

inline Edgelist BuildEdgeListFromCSR(VertexRange vertex_range, const XadjArray& xadj, const AdjncyArray& adjncy) {
    Edgelist edge_list;
    edge_list.reserve(adjncy.size());

    for (SInt u = 0; u + 1 < xadj.size(); ++u) {
        for (SInt e = xadj[u]; e < xadj[u + 1]; ++e) {
            edge_list.emplace_back(u + vertex_range.first, adjncy[e]);
        }
    }

    return edge_list;
}
} // namespace kagen
