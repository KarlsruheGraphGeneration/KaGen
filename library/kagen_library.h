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

using SInt          = unsigned long long;
using SSInt         = long long;
using EdgeList      = std::vector<std::tuple<SInt, SInt>>;
using VertexRange   = std::pair<SInt, SInt>;
using PEID          = int;
using HPFloat       = long double;
using LPFloat       = double;
using Coordinates2D = std::vector<std::tuple<HPFloat, HPFloat>>;
using Coordinates3D = std::vector<std::tuple<HPFloat, HPFloat, HPFloat>>;
using Coordinates   = std::pair<Coordinates2D, Coordinates3D>;

//! Result with 2D coordinates
struct KaGenResult2D {
    inline KaGenResult2D(std::tuple<EdgeList, VertexRange, Coordinates> result)
        : edges(std::move(std::get<0>(result))),
          vertex_range(std::get<1>(result)),
          coordinates(std::move(std::get<2>(result).first)) {}

    EdgeList      edges;
    VertexRange   vertex_range;
    Coordinates2D coordinates;
};

//! Result with 3D coordinates
struct KaGenResult3D {
    inline KaGenResult3D(std::tuple<EdgeList, VertexRange, Coordinates> result)
        : edges(std::move(std::get<0>(result))),
          vertex_range(std::get<1>(result)),
          coordinates(std::move(std::get<2>(result).second)) {}

    EdgeList      edges;
    VertexRange   vertex_range;
    Coordinates3D coordinates;
};

//! Result without coordinates
struct KaGenResult {
    inline KaGenResult(std::tuple<EdgeList, VertexRange, Coordinates> result)
        : edges(std::move(std::get<0>(result))),
          vertex_range(std::move(std::get<1>(result))) {}
    inline KaGenResult(EdgeList edges, VertexRange vertex_range)
        : edges(std::move(edges)),
          vertex_range(vertex_range) {}
    inline KaGenResult(KaGenResult2D&& result) : edges(std::move(result.edges)), vertex_range(result.vertex_range) {}
    inline KaGenResult(KaGenResult3D&& result) : edges(std::move(result.edges)), vertex_range(result.vertex_range) {}

    EdgeList    edges;
    VertexRange vertex_range;
};

class KaGen {
public:
    KaGen(MPI_Comm comm);

    KaGen(const KaGen&) = delete;

    KaGen(KaGen&&) noexcept;

    KaGen& operator=(const KaGen&) = delete;

    KaGen& operator=(KaGen&&) noexcept;

    ~KaGen();

    void SetSeed(int seed);

    void EnableUndirectedGraphVerification();

    void EnableBasicStatistics();

    void EnableAdvancedStatistics();

    void SetNumberOfChunks(SInt k);

    KaGenResult GenerateDirectedGMM(SInt n, SInt m, bool self_loops = false);

    KaGenResult GenerateUndirectedGNM(SInt n, SInt m, bool self_loops = false);

    KaGenResult GenerateDirectedGNP(SInt n, LPFloat p, bool self_loops = false);

    KaGenResult GenerateUndirectedGNP(SInt n, LPFloat p, bool self_loops = false);

    KaGenResult GenerateRGG2D(SInt n, LPFloat r);

    KaGenResult GenerateRGG2D_NM(SInt n, SInt m);

    KaGenResult GenerateRGG2D_MR(SInt m, LPFloat r);

    KaGenResult2D GenerateRGG2D_Coordinates(SInt n, LPFloat r);

    KaGenResult GenerateRGG3D(SInt n, LPFloat r);

    KaGenResult GenerateRGG3D_NM(SInt n, SInt m);

    KaGenResult GenerateRGG3D_MR(SInt m, LPFloat r);

    KaGenResult3D GenerateRGG3D_Coordinates(SInt n, LPFloat r);

    KaGenResult GenerateRDG2D(SInt n, bool periodic);

    KaGenResult GenerateRDG2D_M(SInt m, bool periodic);

    KaGenResult2D GenerateRDG2D_Coordinates(SInt n, bool periodic);

    KaGenResult GenerateRDG3D(SInt n);

    KaGenResult GenerateRDG3D_M(SInt m);

    KaGenResult3D GenerateRDG3D_Coordinates(SInt n);

    KaGenResult GenerateBA(SInt n, SInt d);

    KaGenResult GenerateRHG(LPFloat gamma, SInt n, LPFloat d, bool high_res = false);

    KaGenResult GenerateRHG_NM(LPFloat gamma, SInt n, SInt m, bool high_res = false);

    KaGenResult GenerateRHG_MD(LPFloat gamma, SInt m, LPFloat d, bool high_res = false);

    KaGenResult2D GenerateRHG_Coordinates(LPFloat gamma, SInt n, LPFloat d, bool high_res = false);

    KaGenResult2D GenerateRHG_Coordinates_NM(LPFloat gamma, SInt n, SInt m, bool high_res = false);

    KaGenResult2D GenerateRHG_Coordinates_MD(LPFloat gamma, SInt m, LPFloat d, bool high_res = false);

    KaGenResult GenerateGrid2D(SInt grid_x, SInt grid_y, LPFloat p, bool periodic = false);

    KaGenResult GenerateGrid2D_N(SInt n, LPFloat p, bool periodic = false);

    KaGenResult2D GenerateGrid2D_Coordinates(SInt grid_x, SInt grid_y, LPFloat p, bool periodic = false);

    KaGenResult GenerateGrid3D(SInt grid_x, SInt grid_y, SInt grid_z, LPFloat p, bool periodic = false);

    KaGenResult GenerateGrid3D_N(SInt n, LPFloat p, bool periodic = false);

    KaGenResult3D GenerateGrid3D_Coordinates(SInt grid_x, SInt grid_y, SInt grid_z, LPFloat p, bool periodic = false);

    KaGenResult GenerateKronecker(SInt n, SInt m);

private:
    void SetDefaults();

    MPI_Comm                          comm_;
    std::unique_ptr<PGeneratorConfig> config_;
};

template <typename Graph = KaGenResult, typename IDX = SInt>
void BuildVertexDistribution(Graph& graph, IDX** vtxdist, IDX* vtxdist_size, MPI_Datatype idx_mpi_type, MPI_Comm comm) {
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

template <typename Graph = KaGenResult, typename IDX = SInt>
void BuildCSR(Graph& graph, IDX** xadj, IDX* xadj_size, IDX** adjncy, IDX* adjncy_size) {
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
