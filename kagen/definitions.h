/*******************************************************************************
 * include/definitions.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include <tuple>
#include <vector>

namespace kagen {
// Constants
using LONG  = long long;
using ULONG = unsigned long long;
using INT   = int;
using UINT  = unsigned int;
using PEID  = int;

#ifdef KAGEN_32BIT
using SInt  = UINT;
using SSInt = INT;
#else
using SInt  = ULONG;
using SSInt = LONG;
#endif

// High/low prec
using HPFloat = long double;
using LPFloat = double;

constexpr PEID ROOT = 0;
#define KAGEN_MPI_LONG    MPI_LONG_LONG
#define KAGEN_MPI_ULONG   MPI_UNSIGNED_LONG_LONG
#define KAGEN_MPI_INT     MPI_INT
#define KAGEN_MPI_UINT    MPI_UNSIGNED
#define KAGEN_MPI_PEID    MPI_INT
#define KAGEN_MPI_HPFLOAT MPI_LONG_DOUBLE
#define KAGEN_MPI_LPFLOAT MPI_DOUBLE

#ifdef KAGEN_32BIT
    #define KAGEN_MPI_SINT  MPI_UNSIGNED_INT
    #define KAGEN_MPI_SSINT MPI_INT
#else
    #define KAGEN_MPI_SINT  MPI_UNSIGNED_LONG_LONG
    #define KAGEN_MPI_SSINT MPI_LONG_LONG
#endif

enum Direction { Up, Down, Left, Right, Front, Back };

using EdgeList    = std::vector<std::tuple<SInt, SInt>>;
using VertexRange = std::pair<SInt, SInt>;

using Coordinates2D = std::vector<std::tuple<HPFloat, HPFloat>>;
using Coordinates3D = std::vector<std::tuple<HPFloat, HPFloat, HPFloat>>;
using Coordinates   = std::pair<Coordinates2D, Coordinates3D>;

using VertexWeights = std::vector<SSInt>;
using EdgeWeights   = std::vector<SSInt>;

using XadjArray   = std::vector<SInt>;
using AdjncyArray = std::vector<SInt>;

constexpr std::size_t NEWTON_MAX_ITERS = 10000;
constexpr double      NEWTON_EPS       = 0.001;

enum class GraphRepresentation {
    EDGE_LIST,
    CSR,
};

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
