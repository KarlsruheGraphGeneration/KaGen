/*******************************************************************************
 interface/kagen_interface.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <algorithm>
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

template <typename IDX = SInt>
void BuildVertexDistribution(
    KaGenResult& graph, IDX** vtxdist, IDX* vtxdist_size, MPI_Datatype idx_mpi_type, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    *vtxdist = new IDX[size + 1];
    if (vtxdist_size != nullptr) {
        *vtxdist_size = size + 1;
    }
    (*vtxdist)[rank + 1] = graph.vertex_range.second;

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, (*vtxdist) + 1, 1, idx_mpi_type, comm);
}

template <typename IDX = SInt>
void BuildCSR(KaGenResult& graph, IDX** xadj, IDX* xadj_size, IDX** adjncy, IDX* adjncy_size) {
    // Edges must be sorted
    if (!std::is_sorted(graph.edges.begin(), graph.edges.end())) {
        std::sort(graph.edges.begin(), graph.edges.end());
    }

    const SInt num_local_nodes = graph.vertex_range.second - graph.vertex_range.first;
    const SInt num_local_edges = graph.edges.size();

    // Alloocate CSR data structures
    *xadj   = new IDX[num_local_nodes + 1];
    *adjncy = new IDX[num_local_edges];
    if (xadj_size != nullptr) {
        *xadj_size = num_local_nodes + 1;
    }
    if (adjncy_size != nullptr) {
        *adjncy_size = num_local_edges;
    }

    // Build CSR graph
    SInt cur_vertex = 0;
    SInt cur_edge   = 0;
    (*xadj)[0]      = 0;
    for (const auto& [from, to]: graph.edges) {
        while (from - graph.vertex_range.first > cur_vertex) {
            (*xadj)[++cur_vertex] = cur_edge;
        }
        (*adjncy)[cur_edge++] = to;
    }
    while (cur_vertex < num_local_nodes) {
        (*xadj)[++cur_vertex] = cur_edge;
    }
}
} // namespace kagen
