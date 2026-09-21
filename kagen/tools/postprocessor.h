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
 *        to PE v % size) before balancing. If false, the consecutive balanced vertex distribution
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
 *        to PE v % size) before balancing. If false, the consecutive balanced vertex distribution
 *        [0, n/p), [n/p, 2n/p), ... is used directly without remapping vertex IDs.
 * @param comm The MPI communicator.
 * @return The vertex range assigned to this PE after redistribution.
 */
VertexRange
RedistributeEdgesBalanced(Edgelist& source, Edgelist& destination, SInt n, bool remap_round_robin, MPI_Comm comm);

/**
 * @brief Result of RedistributeEdgesTrueBalance().
 */
struct EdgeBalancedDistribution {
    //! A complete, gap-free partition of [0, n): every vertex in the graph, including isolated (degree-0)
    //! vertices, is covered by exactly one PE's range, so summing NumberOfLocalVertices()-equivalent sizes
    //! across all PEs always yields exactly n. This is the right range to use for Graph::vertex_range: headers,
    //! global-vertex-count, and per-vertex-weight-generation code all need exactly this "every vertex counted
    //! once" property, which the graph-generation and I/O guard rails rely on whether or not any split occurred.
    //! It does *not* imply this PE holds every one of a vertex's edges, though: a shared boundary vertex (see
    //! left_partial_vertex/right_partial_vertex) is resolved to exactly one (the lower-rank) of the two PEs
    //! sharing it for counting purposes, but the higher-rank PE still physically holds some of that vertex's
    //! edges too -- see fully_owned_vertex_range for a range that excludes this case.
    VertexRange vertex_range;

    //! The strict subset of vertex_range whose vertices are *exclusively* owned by this PE: every edge of every
    //! vertex in this range is present on this PE and nowhere else (isolated, degree-0 vertices trivially
    //! qualify). A shared boundary vertex -- see left_partial_vertex/right_partial_vertex -- is excluded from
    //! *both* PEs sharing it, since neither alone holds its whole adjacency. Always a subset of vertex_range
    //! (equal to it whenever has_split_vertices is false); use this instead of vertex_range for anything that
    //! assumes single-PE ownership of a vertex's whole adjacency (e.g. treating a local id range as directly
    //! CSR-indexable).
    VertexRange fully_owned_vertex_range;

    //! This PE's first local vertex, if its own edges are shared with (and owned by) the lower-rank neighbor.
    std::optional<SplitVertexInfo> left_partial_vertex;

    //! This PE's last local vertex, if its own edges are shared with the higher-rank neighbor. Owned by this
    //! (lower-rank) PE. A single vertex whose degree spans more than one PE's fair share can appear as a
    //! boundary vertex on more than 2 PEs in total: as right_partial_vertex on the lowest-rank PE holding a
    //! share of it, as left_partial_vertex on the highest-rank one, and as both on every PE strictly in between.
    std::optional<SplitVertexInfo> right_partial_vertex;

    //! True if any vertex anywhere in the graph (on any PE) was split, i.e. some PE's adjacency-grouped view of
    //! the graph is no longer valid. False if this redistribution happened to produce no actual split (e.g. no
    //! vertex has a degree exceeding a PE's fair edge share), in which case the graph is fully compatible with
    //! adjacency-grouped consumers despite BALANCE_EDGES_TRUE having been requested.
    bool has_split_vertices = false;
};

/**
 * @brief Which of this PE's two boundary vertices are shared with an adjacent PE, plus the resulting
 * gap-free vertex ownership. Computed purely from every PE's (first tail, last tail) via one MPI_Allgather
 * plus one MPI_Allreduce(LOR); carries no edge payload, so it is reusable by both the redistribution path
 * and the direct edge-range read path.
 */
struct BoundaryOwnership {
    //! Complete, gap-free partition of [0, n): boundary vertex resolved to the lower-rank PE; absorbs runs
    //! of isolated (degree-0) vertices into the preceding PE. See EdgeBalancedDistribution::vertex_range.
    VertexRange vertex_range;
    //! Strict subset of vertex_range excluding a boundary vertex shared with a neighbor from both sides.
    VertexRange fully_owned_vertex_range;
    //! This PE's first local vertex is shared with (owned by) the left neighbor.
    bool is_left_partial = false;
    //! This PE's last local vertex is shared with the right neighbor (and owned by this PE).
    bool is_right_partial = false;
    //! True if any vertex anywhere in the graph was split across PEs (MPI_Allreduce(LOR) of the above).
    bool has_split_vertices = false;
};

/**
 * @brief Computes gap-free vertex ownership and split-boundary flags from this PE's first/last tail vertex.
 * One MPI_Allgather of (first_tail, last_tail) + one MPI_Allreduce(LOR); no edge movement. When
 * has_local_edges is false, pass any first_tail/last_tail (they are ignored in favor of the sentinel n).
 *
 * @param first_tail The tail vertex of this PE's first local edge (ignored if has_local_edges is false).
 * @param last_tail The tail vertex of this PE's last local edge (ignored if has_local_edges is false).
 * @param has_local_edges Whether this PE holds any local edges.
 * @param n The number of vertices in the graph.
 * @param comm The MPI communicator.
 */
BoundaryOwnership
ComputeBoundaryOwnership(SInt first_tail, SInt last_tail, bool has_local_edges, SInt n, MPI_Comm comm);

/**
 * @brief Derives the full EdgeBalancedDistribution (vertex_range, fully_owned_vertex_range, left_partial_vertex,
 * right_partial_vertex, has_split_vertices) for a PE that already holds its sorted, edge-balanced local edge
 * list. Wraps ComputeBoundaryOwnership and fills left_partial_vertex/right_partial_vertex via binary search over
 * the (tail-sorted) edges. No edge movement occurs.
 *
 * @param sorted_local_edges This PE's local edges, sorted by (tail, head).
 * @param n The number of vertices in the graph.
 * @param comm The MPI communicator.
 */
EdgeBalancedDistribution ComputeEdgeBalancedBoundaries(const Edgelist& sorted_local_edges, SInt n, MPI_Comm comm);

/**
 * @brief Redistributes edges so that each PE owns exactly (+/-1) the same number of edges, splitting a single
 * vertex's own edges across multiple PEs if its degree exceeds one PE's fair share. Unlike
 * RedistributeEdgesBalanced(), this breaks the invariant that a single PE owns a vertex's whole adjacency; only
 * use this when consuming the result as a plain edge list (not with adjacency-grouped output formats,
 * --validate-simple-graph, or edge-weight/statistics code that assumes single ownership).
 * The source edge list is sorted, deduplicated, and consumed by this call.
 *
 * @param source The edge list to redistribute (sorted, deduplicated, and consumed by this call).
 * @param destination The edge list to store the redistributed edges in.
 * @param n The number of vertices in the graph.
 * @param remap_round_robin If true, vertices are first remapped round-robin (vertex v is assigned
 *        to PE v % size) before balancing. If false, the consecutive balanced vertex distribution
 *        [0, n/p), [n/p, 2n/p), ... is used directly without remapping vertex IDs.
 * @param comm The MPI communicator.
 * @return Description of this PE's resulting vertex ownership; see EdgeBalancedDistribution.
 */
EdgeBalancedDistribution
RedistributeEdgesTrueBalance(Edgelist& source, Edgelist& destination, SInt n, bool remap_round_robin, MPI_Comm comm);
} // namespace kagen
