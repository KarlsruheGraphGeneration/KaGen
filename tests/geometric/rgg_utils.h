#pragma once

#include <gtest/gtest.h>

#include <mpi.h>

#include <cmath>
#include <numeric>
#include <utility>

#include "kagen/context.h"

namespace kagen::testing {

inline void f () {

}

inline Edgelist CreateExpectedRGG2DEdges(PGeneratorConfig config, const Graph& graph) {
    std::vector<std::pair<SInt, SInt>> edge_List;
    for (SInt i = 0; i < config.n; i++) {
        auto [x1, y1] = graph.coordinates.first[i];
        for (SInt j = 0; j < config.n; j++) {
            auto [x2, y2] = graph.coordinates.first[j];
            // Comparing all coordinates
            if (i != j && std::hypot(x1 - x2, y1 - y2) <= config.r) {
                edge_List.emplace_back(i, j);
            }
        }
    }
    return edge_List;
}

inline Edgelist CreateExpectedRGG3DEdges(PGeneratorConfig config, const Graph& graph) {
    std::vector<std::pair<SInt, SInt>> edge_List;
    for (SInt i = 0; i < config.n; i++) {
        auto [x1, y1, z1] = graph.coordinates.second[i];
        for (SInt j = 0; j < config.n; j++) {
            auto [x2, y2, z2] = graph.coordinates.second[j];
            // Comparing all coordinates
            if (i != j && std::hypot(x1 - x2, y1 - y2, z1 - z2) < config.r) {
                edge_List.emplace_back(i, j);
            }
        }
    }
    return edge_List;
}

} // namespace kagen::testing
