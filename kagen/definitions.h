/*******************************************************************************
 * include/definitions.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include <ostream>
#include <tuple>
#include <vector>

#include "kagen/interface_definitions.h"

namespace kagen {

constexpr PEID ROOT = 0;
#define KAGEN_MPI_LONG    MPI_LONG_LONG
#define KAGEN_MPI_ULONG   MPI_UNSIGNED_LONG_LONG
#define KAGEN_MPI_INT     MPI_INT
#define KAGEN_MPI_UINT    MPI_UNSIGNED
#define KAGEN_MPI_PEID    MPI_INT
#define KAGEN_MPI_HPFLOAT MPI_LONG_DOUBLE
#define KAGEN_MPI_LPFLOAT MPI_DOUBLE
#define KAGEN_MPI_SINT    MPI_UNSIGNED_LONG_LONG
#define KAGEN_MPI_SSINT   MPI_LONG_LONG

enum Direction { Up, Down, Left, Right, Front, Back };

constexpr std::size_t NEWTON_MAX_ITERS = 10000;
constexpr double      NEWTON_EPS       = 0.001;

inline std::ostream& operator<<(std::ostream& out, const GraphRepresentation representation) {
    switch (representation) {
        case GraphRepresentation::EDGE_LIST:
            return out << "edge-list";

        case GraphRepresentation::CSR:
            return out << "csr";
    }

    return out << "<invalid>";
}

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
