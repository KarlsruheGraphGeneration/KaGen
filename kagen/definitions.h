#pragma once

#include <ostream>
#include <tuple>
#include <vector>

#include "kagen/kagen.h"

namespace kagen {
constexpr PEID ROOT = 0;

enum Direction { Up, Down, Left, Right, Front, Back };

constexpr std::size_t NEWTON_MAX_ITERS = 10000;
constexpr double      NEWTON_EPS       = 0.001;

struct Graph {
    VertexRange         vertex_range;
    GraphRepresentation representation;

    // Edge list representation only
    EdgeList edges;

    // CSR representation only
    XadjArray   xadj;
    AdjncyArray adjncy;

    VertexWeights vertex_weights;
    EdgeWeights   edge_weights;
    Coordinates   coordinates;

    std::tuple<VertexRange, EdgeList, XadjArray, AdjncyArray, VertexWeights, EdgeWeights, Coordinates> tuple() && {
        return std::make_tuple(
            vertex_range, std::move(edges), std::move(xadj), std::move(adjncy), std::move(vertex_weights),
            std::move(edge_weights), std::move(coordinates));
    }
};
} // namespace kagen
