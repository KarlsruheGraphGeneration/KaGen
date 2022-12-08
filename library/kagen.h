/*******************************************************************************
 * interface/kagen_interface.h
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
using VertexWeights = std::vector<SSInt>;
using EdgeWeights   = std::vector<SSInt>;

//! Result with 2D coordinates
struct KaGenResult2D {
    inline KaGenResult2D() : edges(), vertex_range(0, 0), coordinates() {}
    inline KaGenResult2D(std::tuple<EdgeList, VertexRange, Coordinates, VertexWeights, EdgeWeights> result)
        : edges(std::move(std::get<0>(result))),
          vertex_range(std::get<1>(result)),
          coordinates(std::move(std::get<2>(result).first)) {}

    EdgeList      edges;
    VertexRange   vertex_range;
    Coordinates2D coordinates;
};

//! Result with 3D coordinates
struct KaGenResult3D {
    inline KaGenResult3D() : edges(), vertex_range(0, 0), coordinates() {}
    inline KaGenResult3D(std::tuple<EdgeList, VertexRange, Coordinates, VertexWeights, EdgeWeights> result)
        : edges(std::move(std::get<0>(result))),
          vertex_range(std::get<1>(result)),
          coordinates(std::move(std::get<2>(result).second)) {}

    EdgeList      edges;
    VertexRange   vertex_range;
    Coordinates3D coordinates;
};

//! Result without coordinates
struct KaGenResult {
    inline KaGenResult() : edges(), vertex_range(0, 0) {}
    inline KaGenResult(std::tuple<EdgeList, VertexRange, Coordinates, VertexWeights, EdgeWeights> result)
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

    /*!
     * Sets the seed for the random number generator (Default: 0).
     * @param seed Seed for random number generation.
     */
    void SetSeed(int seed);

    /*!
     * If enabled, KaGen will verify that the generated graph is simple and undirected.
     * This feature is only for debugging; unless explicitly configured otherwise, KaGen
     * should *never* fail this verification.
     */
    void EnableUndirectedGraphVerification();

    void EnableBasicStatistics();

    void EnableAdvancedStatistics();

    /*!
     * If enabled, KaGen will print information to stdout and stderr (but only on rank 0).
     *
     * @param header If set to true, KaGen will also print a banner and a summary of the
     * configuration parameters.
     */
    void EnableOutput(bool header);

    /*!
     * If set to true, KaGen will use higher-precision floating point numbers (80 bit instead of 64 bit on x86 systems)
     * to generate graphs. Currently, this only affects the random hyperbolic graph generator. Per default, KaGen will
     * decide automatically which precision to use.
     *
     * Note that "higher-precision" refers to "long double". On systems where "long double" has the same precision as
     * "double", this option does nothing.
     *
     * @param state If true, always use higher-precision floating point numbers; if false, never use them.
     */
    void UseHPFloats(bool state);

    /*!
     * Controls the number of chunks KaGen uses for graph generation. If not set explicitly, KaGen will choose the
     * number of chunks automatically. Usually, you do not have to use this option.
     *
     * @param k Number of chunks to be used.
     */
    void SetNumberOfChunks(SInt k);

    /*!
     * Generates a graph with options given by a string of options in `key=value` or `flag` format:
     * `key1=value1;flag1;...`.
     *
     * Use one of the following flags to select a graph model:
     * - gnm_undirected
     * - gnm_directed
     * - gnp_undirected
     * - gnp_directed
     * - rgg2d
     * - rgg3d
     * - grid2d
     * - grid3d
     * - rdg2d
     * - rdg3d
     * - rhg
     * - ba
     * - kronecker
     * - rmat
     *
     * Use the following keys to specify generator properties:
     * - n=<SInt>             -- number of nodes
     * - N=<SInt>             -- number of nodes as a power of 2
     * - m=<SInt>             -- number of edges
     * - M=<SInt>             -- number of edges as a power of 2
     * - k=<SInt>             -- number of chunks
     * - prob=<HPFloat>       -- edge probability (varius generators)
     * - radius=<HPFloat>     -- edge radius (RGG2D/3D)
     * - gamma=<HPFloat>      -- power law exponent (RHG)
     * - avg_degree=<HPFloat> -- average degree (RHG)
     * - min_degree=<SInt>    -- minimum degree (BA)
     * - grid_x=<SInt>        -- grid width (GRID2D/3D)
     * - grid_y=<SInt>        -- grid height (GRID2D/3D)
     * - grid_z=<SInt>        -- grid depth (GRID3D)
     * - rmat_a=<HPFloat>     -- RMat probability for block A (RMAT)
     * - rmat_b=<HPFloat>     -- RMat probability for block B (RMAT)
     * - rmat_c=<HPFloat>     -- RMat probability for block C (RMAT)
     * - periodic[=0|1]       -- periodic boundary condition (various generators)
     *
     * Depending on the selected generator type, some options are mandatory, some are optional and some are ignored.
     * The following example generates a RGG2D graph with 100 nodes and 200 edges: `rgg2d;n=100;m=200`.
     *
     * @type options Options string with key=value pairs.
     * @return The generated graph.
     */
    KaGenResult GenerateFromOptionString(const std::string& options);

    KaGenResult2D GenerateFromOptionString2D(const std::string& options);

    KaGenResult3D GenerateFromOptionString3D(const std::string& options);

    KaGenResult GenerateDirectedGNM(SInt n, SInt m, bool self_loops = false);

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

    KaGenResult GenerateBA(SInt n, SInt d, bool directed = false, bool self_loops = false);

    KaGenResult GenerateBA_NM(SInt n, SInt m, bool directed = false, bool self_loops = false);

    KaGenResult GenerateBA_MD(SInt m, SInt d, bool directed = false, bool self_loops = false);

    KaGenResult GenerateRHG(LPFloat gamma, SInt n, LPFloat d);

    KaGenResult GenerateRHG_NM(LPFloat gamma, SInt n, SInt m);

    KaGenResult GenerateRHG_MD(LPFloat gamma, SInt m, LPFloat d);

    KaGenResult2D GenerateRHG_Coordinates(LPFloat gamma, SInt n, LPFloat d);

    KaGenResult2D GenerateRHG_Coordinates_NM(LPFloat gamma, SInt n, SInt m);

    KaGenResult2D GenerateRHG_Coordinates_MD(LPFloat gamma, SInt m, LPFloat d);

    KaGenResult GenerateGrid2D(SInt grid_x, SInt grid_y, LPFloat p, bool periodic = false);

    KaGenResult GenerateGrid2D_N(SInt n, LPFloat p, bool periodic = false);

    KaGenResult GenerateGrid2D_NM(SInt n, SInt m, bool periodic = false);

    KaGenResult2D GenerateGrid2D_Coordinates(SInt grid_x, SInt grid_y, LPFloat p, bool periodic = false);

    KaGenResult GenerateGrid3D(SInt grid_x, SInt grid_y, SInt grid_z, LPFloat p, bool periodic = false);

    KaGenResult GenerateGrid3D_N(SInt n, LPFloat p, bool periodic = false);

    KaGenResult GenerateGrid3D_NM(SInt n, SInt m, bool periodic = false);

    KaGenResult3D GenerateGrid3D_Coordinates(SInt grid_x, SInt grid_y, SInt grid_z, LPFloat p, bool periodic = false);

    KaGenResult GenerateKronecker(SInt n, SInt m, bool directed = false, bool self_loops = false);

    KaGenResult
    GenerateRMAT(SInt n, SInt m, LPFloat a, LPFloat b, LPFloat c, bool directed = false, bool self_loops = false);

private:
    void SetDefaults();

    MPI_Comm                          comm_;
    std::unique_ptr<PGeneratorConfig> config_;
};

/*!
 * Returns an array A of size |comm| + 1 s.t. PE i contains nodes in the range [A[i], A[i + 1]).
 *
 * @param graph Graph generated by any of the Generate* functions.
 * @param idx_mpi_type MPI type corresponding to the template parameter IDX.
 * @param comm MPI communicator that was used to generate the graph.
 *
 * @tparam IDX Data type to be used for the entries in A. Must be large enough to represent the global number of nodes
 * in the graph.
 * @tparam Graph Return type of the Generate* function that generated the graph.
 *
 * @return Vertex distribution as described above.
 */
template <typename IDX, typename Graph>
std::vector<IDX> BuildVertexDistribution(const Graph& graph, MPI_Datatype idx_mpi_type, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::vector<IDX> distribution(size + 1);
    distribution[rank + 1] = graph.vertex_range.second;
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, distribution.data() + 1, 1, idx_mpi_type, comm);

    return distribution;
}

template <typename IDX>
struct KaGenResultCSR {
    std::vector<IDX> xadj;
    std::vector<IDX> adjncy;
};

/*!
 * Transforms a graph generated by any of the Generate* functions (i.e., in edge list format) into the compressed sparse
 * row (CSR) format.
 *
 * Let n and m be the number of vertices and the number of edges on this PE, respectively. Then, the result of this
 * transformation consists of two arrays: xadj of size n + 1 and adjncy of size m.
 * For a local node u (i.e., 0 <= u < n), its degree is xadj[u + 1] - xadj[u] and its neighbors are adjncy[xadj[u], ...,
 * xadj[u + 1]). Node that adjncy[xadj[u + 1]] is excluded, i.e., is *not* a neighbor of node u.
 * The IDs of neighbors are global node IDs, i.e., not necessarily in the interval [0, n).
 *
 * @param graph Graph generated by any of the Generate* functions. Use `std::move()` if you do not need the edge list
 * after transformation.
 *
 * @tparam IDX Data type to be used for vertex and edge IDs in the transformed graph. Must be large enough to represent
 * the number of vertices and the number of edges in the global graph.
 * @tparam Graph Return type of the Generate* function that was used to generate the graph.
 *
 * @return Graph in CSR format.
 */
template <typename IDX, typename Graph>
KaGenResultCSR<IDX> BuildCSR(Graph graph) {
    // Edges must be sorted
    if (!std::is_sorted(graph.edges.begin(), graph.edges.end())) {
        std::sort(graph.edges.begin(), graph.edges.end());
    }

    const SInt num_local_nodes = graph.vertex_range.second - graph.vertex_range.first;
    const SInt num_local_edges = graph.edges.size();

    KaGenResultCSR<IDX> csr;
    csr.xadj.resize(num_local_nodes + 1);
    csr.adjncy.resize(num_local_edges);

    // Build CSR graph
    SInt cur_vertex = 0;
    SInt cur_edge   = 0;

    for (const auto& [from, to]: graph.edges) {
        while (from - graph.vertex_range.first > cur_vertex) {
            csr.xadj[++cur_vertex] = cur_edge;
        }
        csr.adjncy[cur_edge++] = to;
    }
    while (cur_vertex < num_local_nodes) {
        csr.xadj[++cur_vertex] = cur_edge;
    }

    return csr;
}
} // namespace kagen

