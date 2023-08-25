#pragma once

#include <gtest/gtest.h>

#include <mpi.h>

#include <cmath>
#include <numeric>
#include <utility>

#include "kagen/context.h"

namespace kagen::testing {

    inline void CheckReverseEdges(Edgelist edges, EdgeWeights weights) {
        int left, right, id;
        for (int i = 0; i < edges.size(); i++) {
            SInt tail = edges[i].first;
            SInt head = edges[i].second;

            if (tail < head) {
                left = 0;
                right = edges.size() - 1;

                while (left <= right) {
                    id = left + (right - left) / 2;

                    if (edges[id].first == head && edges[id].second == tail) {
                        break;
                    } else if (edges[id].first < head ||
                               (edges[id].first == head && edges[id].second < tail)) {
                        left = id + 1;
                    } else {
                        right = id - 1;
                    }
                }

                std::cout << "Edge i (" << edges[i].first << ", " << edges[i].second << ") " << weights[i] << " Edge id ("
                          << edges[id].first << ", " << edges[id].second << ") " << weights[id] << std::endl;


                ASSERT_EQ(weights[i], weights[id]);
            }
        }


//        std::vector<std::pair<SInt, SInt>> edge_List;
//        for (SInt i = 0; i < config.n; i++) {
//            auto [x1, y1] = graph.coordinates.first[i];
//            for (SInt j = 0; j < config.n; j++) {
//                auto [x2, y2] = graph.coordinates.first[j];
//                // Comparing all coordinates
//                if (i != j && std::hypot(x1 - x2, y1 - y2) <= config.r) {
//                    edge_List.emplace_back(i, j);
//                }
//            }
//        }
//        return edge_List;
    }

} // namespace kagen::testing
