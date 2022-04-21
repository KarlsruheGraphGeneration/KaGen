/*******************************************************************************
 interface/kagen_interface.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace kagen {
struct PGeneratorConfig;

using SInt     = unsigned long long;
using SSInt    = long long;
using EdgeList = std::vector<std::tuple<SInt, SInt>>;
using PEID     = int;
using HPFloat  = long double;
using LPFloat  = double;

struct KaGenResult {
    EdgeList              edges;
    std::pair<SInt, SInt> vertex_range;
};

class KaGen {
public:
    KaGen(PEID rank, PEID size);
    ~KaGen();

    void SetSeed(int seed);

    void EnableUndirectedGraphVerification();

    KaGenResult GenerateDirectedGMM(SInt n, SInt m, SInt k = 0, bool self_loops = false);

    KaGenResult GenerateUndirectedGNM(SInt n, SInt m, SInt k = 0, bool self_loops = false);

    KaGenResult GenerateDirectedGNP(SInt n, LPFloat p, SInt k = 0, bool self_loops = false);

    KaGenResult GenerateUndirectedGNP(SInt n, LPFloat p, SInt k = 0, bool self_loops = false);

    KaGenResult Generate2DRGG(SInt n, LPFloat r, SInt k = 0);

    KaGenResult Generate3DRGG(SInt n, LPFloat r, SInt k = 0);

    KaGenResult Generate2DRDG(SInt n, SInt k = 0);

    KaGenResult Generate3DRDG(SInt n, SInt k = 0);

    KaGenResult GenerateBA(SInt n, SInt d, SInt k = 0);

    KaGenResult GenerateRHG(SInt n, LPFloat gamma, SInt d, SInt k = 0);

    KaGenResult Generate2DGrid(SInt n, LPFloat p, bool periodic, SInt k = 0);

    KaGenResult Generate2DGrid(SInt grid_x, SInt grid_y, LPFloat p, bool periodic, SInt k = 0);

    KaGenResult Generate3DGrid(SInt grid_x, SInt grid_y, SInt grid_z, LPFloat p, bool periodic, SInt k = 0);

    KaGenResult GenerateKronecker(SInt n, SInt m, SInt k = 0);

private:
    void SetDefaults();

    PEID                              rank_, size_;
    std::unique_ptr<PGeneratorConfig> config_;
    bool                              validate_undirected_graph_;
};
} // namespace kagen
