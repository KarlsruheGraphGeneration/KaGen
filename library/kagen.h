/*******************************************************************************
 * interface/kagen_interface.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include <algorithm>
#include <iterator>
#include <memory>
#include <numeric>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <mpi.h>

#include "kagen/interface_definitions.h"

namespace kagen {
struct PGeneratorConfig;

struct KaGenResult {
    inline KaGenResult() : vertex_range(0, 0) {}
    inline KaGenResult(
        std::tuple<VertexRange, EdgeList, XadjArray, AdjncyArray, VertexWeights, EdgeWeights, Coordinates> result)
        : vertex_range(std::move(std::get<0>(result))),
          edges(std::move(std::get<1>(result))),
          xadj(std::move(std::get<2>(result))),
          adjncy(std::move(std::get<3>(result))),
          vertex_weights(std::move(std::get<4>(result))),
          edge_weights(std::move(std::get<5>(result))),
          coordinates_2d(std::move(std::get<6>(result).first)),
          coordinates_3d(std::move(std::get<6>(result).second)) {}

    template <typename T = SInt>
    std::vector<std::pair<T, T>> TakeEdges() {
        return TakeVector<std::pair<T, T>>(edges);
    }

    template <typename T = SInt>
    std::vector<T> TakeXadj() {
        return TakeVector<T>(xadj);
    }

    template <typename T = SInt>
    std::vector<T> TakeAdjncy() {
        return TakeVector<T>(adjncy);
    }

    template <typename T = SSInt>
    std::vector<T> TakeVertexWeights() {
        return TakeVector<T>(vertex_weights);
    }

    template <typename T = SSInt>
    std::vector<T> TakeEdgeWeights() {
        return TakeVector<T>(edge_weights);
    }

    VertexRange vertex_range;

    // Edge list representation
    EdgeList edges;

    // CSR representation
    XadjArray   xadj;
    AdjncyArray adjncy;

    VertexWeights vertex_weights;
    EdgeWeights   edge_weights;

    Coordinates2D coordinates_2d;
    Coordinates3D coordinates_3d;

private:
    template <typename To, typename From, std::enable_if_t<std::is_same_v<typename From::value_type, To>, bool> = true>
    std::vector<To> TakeVector(From& from) {
        return std::move(from);
    }

    template <typename To, typename From, std::enable_if_t<!std::is_same_v<typename From::value_type, To>, bool> = true>
    std::vector<To> TakeVector(From& from) {
        std::vector<To> copy(from.size());
        std::copy(from.begin(), from.end(), copy.begin());
        std::vector<typename From::value_type> free;
        std::swap(from, free);
        return copy;
    }
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
     * Must be the same on all PEs.
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
     * Represents the generated graph as a list of edges (from, to).
     * This representation requires 2 * |E| memory and is the default representation.
     */
    void UseEdgeListRepresentation();

    /*!
     * Represents the generated graph in compressed sparse row format.
     * This representation requires |V| + |E| memory. However, not all generators support this representation;
     * for generators that do not support it, KaGen first generates the graph as a list of edges, sorts it and
     * then builds the CSR data structures. This requires |V|+3|E| memory and |E|*log(|E|) time.
     */
    void UseCSRRepresentation();

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
     * - coordinates          -- generate coordinates (only for geometric graph generators)
     *
     * Depending on the selected generator type, some options are mandatory, some are optional and some are ignored.
     * The following example generates a RGG2D graph with 100 nodes and 200 edges: `rgg2d;n=100;m=200`.
     *
     * @type options Options string with key=value pairs.
     * @return The generated graph.
     */
    KaGenResult GenerateFromOptionString(const std::string& options);

    KaGenResult GenerateDirectedGNM(SInt n, SInt m, bool self_loops = false);

    KaGenResult GenerateUndirectedGNM(SInt n, SInt m, bool self_loops = false);

    KaGenResult GenerateDirectedGNP(SInt n, LPFloat p, bool self_loops = false);

    KaGenResult GenerateUndirectedGNP(SInt n, LPFloat p, bool self_loops = false);

    KaGenResult GenerateRGG2D(SInt n, LPFloat r, bool coordinates = false);

    KaGenResult GenerateRGG2D_NM(SInt n, SInt m, bool coordinates = false);

    KaGenResult GenerateRGG2D_MR(SInt m, LPFloat r, bool coordinates = false);

    KaGenResult GenerateRGG3D(SInt n, LPFloat r, bool coordinates = false);

    KaGenResult GenerateRGG3D_NM(SInt n, SInt m, bool coordinates = false);

    KaGenResult GenerateRGG3D_MR(SInt m, LPFloat r, bool coordinates = false);

    KaGenResult GenerateRDG2D(SInt n, bool periodic, bool coordinates = false);

    KaGenResult GenerateRDG2D_M(SInt m, bool periodic, bool coordinates = false);

    KaGenResult GenerateRDG3D(SInt n, bool coordinates = false);

    KaGenResult GenerateRDG3D_M(SInt m, bool coordinates = false);

    KaGenResult GenerateBA(SInt n, SInt d, bool directed = false, bool self_loops = false);

    KaGenResult GenerateBA_NM(SInt n, SInt m, bool directed = false, bool self_loops = false);

    KaGenResult GenerateBA_MD(SInt m, SInt d, bool directed = false, bool self_loops = false);

    KaGenResult GenerateRHG(LPFloat gamma, SInt n, LPFloat d, bool coordinates = false);

    KaGenResult GenerateRHG_NM(LPFloat gamma, SInt n, SInt m, bool coordinates = false);

    KaGenResult GenerateRHG_MD(LPFloat gamma, SInt m, LPFloat d, bool coordinates = false);

    KaGenResult GenerateGrid2D(SInt grid_x, SInt grid_y, LPFloat p, bool periodic = false, bool coordinates = false);

    KaGenResult GenerateGrid2D_N(SInt n, LPFloat p, bool periodic = false, bool coordinates = false);

    KaGenResult GenerateGrid2D_NM(SInt n, SInt m, bool periodic = false, bool coordinates = false);

    KaGenResult
    GenerateGrid3D(SInt grid_x, SInt grid_y, SInt grid_z, LPFloat p, bool periodic = false, bool coordinates = false);

    KaGenResult GenerateGrid3D_N(SInt n, LPFloat p, bool periodic = false, bool coordinates = false);

    KaGenResult GenerateGrid3D_NM(SInt n, SInt m, bool periodic = false, bool coordinates = false);

    KaGenResult GenerateKronecker(SInt n, SInt m, bool directed = false, bool self_loops = false);

    KaGenResult
    GenerateRMAT(SInt n, SInt m, LPFloat a, LPFloat b, LPFloat c, bool directed = false, bool self_loops = false);

    KaGenResult ReadFromFile(
        std::string const& filename, const StaticGraphFormat format, const StaticGraphDistribution distribution);

private:
    void SetDefaults();

    MPI_Comm                          comm_;
    std::unique_ptr<PGeneratorConfig> config_;
    GraphRepresentation               representation_;
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
} // namespace kagen

