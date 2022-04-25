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

#include <mpi.h>

namespace kagen {
struct PGeneratorConfig;

using SInt        = unsigned long long;
using SSInt       = long long;
using EdgeList    = std::vector<std::tuple<SInt, SInt>>;
using VertexRange = std::pair<SInt, SInt>;
using PEID        = int;
using HPFloat     = long double;
using LPFloat     = double;

struct KaGenResult {
    // Implicitly convert result of kagen::Generator to KaGenResult
    inline KaGenResult(std::pair<EdgeList, VertexRange> result)
        : edges(std::move(result.first)),
          vertex_range(std::move(result.second)) {}

    EdgeList    edges;
    VertexRange vertex_range;
};

class KaGen {
public:
    KaGen(MPI_Comm comm);
    ~KaGen();

    void SetSeed(int seed);

    void EnableUndirectedGraphVerification();

    void SetNumberOfChunks(SInt k);

    KaGenResult GenerateDirectedGMM(SInt n, SInt m = 0, bool self_loops = false);

    KaGenResult GenerateUndirectedGNM(SInt n, SInt m = 0, bool self_loops = false);

    KaGenResult GenerateDirectedGNP(SInt n, LPFloat p = 0, bool self_loops = false);

    KaGenResult GenerateUndirectedGNP(SInt n, LPFloat p = 0, bool self_loops = false);

    KaGenResult Generate2DRGG(SInt n, LPFloat r = 0);

    KaGenResult Generate3DRGG(SInt n, LPFloat r = 0);

    KaGenResult Generate2DRDG(SInt n = 0);

    KaGenResult Generate3DRDG(SInt n = 0);

    KaGenResult GenerateBA(SInt n, SInt d = 0);

    KaGenResult GenerateRHG(SInt n, LPFloat gamma, SInt d = 0);

    KaGenResult Generate2DGrid(SInt n, LPFloat p, bool periodic = 0);

    KaGenResult Generate2DGrid(SInt grid_x, SInt grid_y, LPFloat p, bool periodic = 0);

    KaGenResult Generate3DGrid(SInt n, LPFloat p, bool periodic = 0);

    KaGenResult Generate3DGrid(SInt grid_x, SInt grid_y, SInt grid_z, LPFloat p, bool periodic = 0);

    KaGenResult GenerateKronecker(SInt n, SInt m = 0);

private:
    void SetDefaults();

    MPI_Comm                          comm_;
    std::unique_ptr<PGeneratorConfig> config_;
};
} // namespace kagen
