#pragma once

#include "kagen/kagen.h"

#include <mpi.h>

namespace kagen {
/**
 * @brief Adds any missing reverse edges that goes between different PEs.
 * More precisely, for any edge `(u, v)` where `u` is assigned to another PE than `v`, this function adds the reverse
 * edge `(v, u)` if it is not already present.
 * **Note:** in order to do so, the edge list is sorted and duplicate edges are removed.
 *
 * @param edge_list The edge list to add reverse edges to.
 * @param edge_weights The corresponding edge weights if generated (otherwise list is empty).
 * @param vertex_range The vertex range assigned to this PE.
 * @param comm The MPI communicator.
 */
void AddNonlocalReverseEdges(Edgelist& edge_list, EdgeWeights& edge_weights, VertexRange vertex_range, MPI_Comm comm);

/**
 * @brief Redistributes the edges such that each PE owns the edges that it should owned due to the given vertex range.
 * @param edge_list The edge list to redistribute.
 * @param vertex_range The vertex range assigned to this PE.
 * @param comm The MPI communicator.
 * @param use_binary_search If true, use binary search for PE lookup (faster for many PEs).
 */
void RedistributeEdgesByVertexRange(
    Edgelist& edge_list, VertexRange vertex_range, MPI_Comm comm, bool use_binary_search = false);

/**
 * @brief Redistributes the edges such that each PE owns a contiguous range of vertices.
 * The source edge list is sorted, deduplicated, and consumed by this call.
 *
 * @param source The edge list to redistribute (sorted, deduplicated, and consumed by this call).
 * @param destination The edge list to store the redistributed edges in.
 * @param n The number of vertices in the graph.
 * @param remap_round_robin If true, vertices are first remapped round-robin (vertex v is assigned
 *        to PE v % size) before redistribution. This breaks correlation between vertex ID structure and
 *        PE assignment, which is important for generators like R-MAT where low-numbered vertices
 *        have disproportionately high degrees. If false, the consecutive balanced vertex distribution
 *        [0, n/p), [n/p, 2n/p), ... is used directly without remapping vertex IDs.
 * @param comm The MPI communicator.
 * @return The vertex range assigned to this PE.
 */
VertexRange RedistributeEdges(Edgelist& source, Edgelist& destination, SInt n, bool remap_round_robin, MPI_Comm comm);

/**
 * @brief Computes a balanced vertex distribution where each PE gets n/p consecutive vertices.
 * The first n%p PEs get one additional vertex.
 *
 * @param n The number of vertices in the graph.
 * @param comm The MPI communicator.
 * @return The vertex distribution array (size + 1 entries): [0, n0, n0+n1, ..., n].
 */
std::vector<SInt> ComputeBalancedVertexDistribution(SInt n, MPI_Comm comm);

/**
 * @brief Remaps vertex IDs using a round-robin assignment and returns the resulting vertex distribution.
 *
 * @param edges The edge list whose vertex IDs are remapped in-place.
 * @param n The number of vertices in the graph.
 * @param comm The MPI communicator.
 * @return The vertex distribution array (size + 1 entries).
 */
std::vector<SInt> RoundRobinRemapping(Edgelist& edges, SInt n, MPI_Comm comm);

/**
 * @brief Redistributes edges to balance the number of edges per PE.
 * First, vertices are assigned to PEs according to a vertex distribution. Then, the vertex
 * distribution is refined so that each PE owns approximately the same number of edges.
 * The source edge list is sorted, deduplicated, and consumed by this call.
 *
 * @param source The edge list to redistribute (sorted, deduplicated, and consumed by this call).
 * @param destination The edge list to store the redistributed edges in.
 * @param n The number of vertices in the graph.
 * @param remap_round_robin If true, vertices are first remapped round-robin (vertex v is assigned
 *        to PE v % size) before balancing. This breaks correlation between vertex ID structure and
 *        PE assignment, which is important for generators like R-MAT where low-numbered vertices
 *        have disproportionately high degrees. If false, the consecutive balanced vertex distribution
 *        [0, n/p), [n/p, 2n/p), ... is used directly without remapping vertex IDs.
 * @param comm The MPI communicator.
 * @return The vertex range assigned to this PE after redistribution.
 */
VertexRange
RedistributeEdgesBalanced(Edgelist& source, Edgelist& destination, SInt n, bool remap_round_robin, MPI_Comm comm);
} // namespace kagen
