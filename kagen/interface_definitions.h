/*******************************************************************************
 * kagen/interface_definitions.h
 *
 * Copyright (C) 2023 Tim Niklas Uhl <uhl@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include <tuple>
#include <vector>

namespace kagen {
enum class GraphRepresentation {
    EDGE_LIST,
    CSR,
};
enum class StaticGraphFormat {
    METIS,
    BINARY_PARHIP,
};

enum class StaticGraphDistribution {
    BALANCE_VERTICES,
    BALANCE_EDGES,
};
using SInt          = unsigned long long;
using SSInt         = long long;
using EdgeList      = std::vector<std::pair<SInt, SInt>>;
using VertexRange   = std::pair<SInt, SInt>;
using PEID          = int;
using HPFloat       = long double;
using LPFloat       = double;
using Coordinates2D = std::vector<std::tuple<HPFloat, HPFloat>>;
using Coordinates3D = std::vector<std::tuple<HPFloat, HPFloat, HPFloat>>;
using Coordinates   = std::pair<Coordinates2D, Coordinates3D>;
using VertexWeights = std::vector<SSInt>;
using EdgeWeights   = std::vector<SSInt>;
using XadjArray     = std::vector<SInt>;
using AdjncyArray   = std::vector<SInt>;
}
