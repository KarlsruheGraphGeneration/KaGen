#include "kagen/tools/postprocessor.h"

#include "kagen/tools/utils.h"

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace kagen {

void AddNonlocalReverseEdges(
    Edgelist& edge_list, EdgeWeights& edge_weights, const VertexRange vertex_range, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const auto ranges           = AllgatherVertexRange(vertex_range, comm);
    const bool has_edge_weights = !edge_weights.empty();

    const auto& [local_from, local_to] = ranges[rank];
    std::unordered_map<PEID, std::vector<SInt>> message_buffers_edges;
    std::unordered_map<PEID, std::vector<SInt>> message_buffers_weights;

    // Each PE gets the edges that we have to that PE
    for (size_t i = 0; i < edge_list.size(); ++i) {
        const auto& [tail, head] = edge_list[i];
        if ((tail >= local_from && tail < local_to) && (head < local_from || head >= local_to)) {
            const SInt pe = static_cast<SInt>(FindPEInRange(head, ranges));
            message_buffers_edges[pe].emplace_back(tail);
            message_buffers_edges[pe].emplace_back(head);
            if (has_edge_weights) {
                message_buffers_weights[pe].emplace_back(edge_weights[i]);
            }
        }
    }

    {
        // exchange edges
        auto recv_buf = ExchangeMessageBuffers(message_buffers_edges, MPI_UINT64_T, comm);

        for (std::size_t i = 0; i < recv_buf.size(); i += 2) {
            edge_list.emplace_back(recv_buf[i + 1], recv_buf[i]);
        }
    }
    {
        // exchange weights
        auto recv_buf = ExchangeMessageBuffers(message_buffers_weights, MPI_INT64_T, comm);

        for (std::size_t i = 0; i < recv_buf.size(); ++i) {
            edge_weights.emplace_back(recv_buf[i]);
        }
    }

    // KaGen sometimes produces duplicate edges
    SortEdgesAndWeights(edge_list, edge_weights);
    RemoveDuplicates(edge_list, edge_weights);
}

void RedistributeEdgesByVertexRange(
    Edgelist& edge_list, const VertexRange vertex_range, MPI_Comm comm, bool use_binary_search) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const auto ranges = AllgatherVertexRange(vertex_range, comm);
    const auto from   = ranges[rank].first;
    const auto to     = ranges[rank].second;

    // Create new edge arrays
    std::vector<std::pair<SInt, SInt>>              local_edges;
    std::vector<std::vector<std::pair<SInt, SInt>>> remote_edges(size);
    for (const auto& [u, v]: edge_list) {
        if (from <= u && u < to) { // Edge starts from local vertex
            local_edges.emplace_back(u, v);
        } else { // Edge starts from remote vertex
            const PEID owner = use_binary_search ? FindPEInRangeWithBinarySearch(u, ranges) : FindPEInRange(u, ranges);
            assert(0 <= owner && owner < size);
            remote_edges[owner].emplace_back(u, v);
        }
    }

    // Exchange edges
    std::vector<SInt> recv_buf;
    std::vector<SInt> send_buf;
    std::vector<int>  send_counts(size);
    std::vector<int>  recv_counts(size);
    std::vector<int>  send_displs(size);
    std::vector<int>  recv_displs(size);
    for (std::size_t i = 0; i < send_counts.size(); ++i) {
        send_counts[i] = remote_edges[i].size() * 2;
    }

    std::exclusive_scan(send_counts.begin(), send_counts.end(), send_displs.begin(), 0);
    const SInt total_send_count = send_displs.back() + send_counts.back();
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm);
    std::exclusive_scan(recv_counts.begin(), recv_counts.end(), recv_displs.begin(), 0);
    const SInt total_recv_count = recv_displs.back() + recv_counts.back();

    send_buf.reserve(total_send_count);
    for (std::size_t i = 0; i < send_counts.size(); ++i) {
        for (const auto& [u, v]: remote_edges[i]) {
            send_buf.push_back(u);
            send_buf.push_back(v);
        }
        {
            [[maybe_unused]] auto _clear = std::move(remote_edges[i]);
        }
    }
    recv_buf.resize(total_recv_count);
    local_edges.reserve(local_edges.size() + total_recv_count / 2);

    MPI_Alltoallv(
        send_buf.data(), send_counts.data(), send_displs.data(), KAGEN_MPI_SINT, recv_buf.data(), recv_counts.data(),
        recv_displs.data(), KAGEN_MPI_SINT, comm);
    {
        [[maybe_unused]] auto _clear = std::move(send_buf);
    }

    for (std::size_t i = 0; i < recv_buf.size(); i += 2) {
        local_edges.emplace_back(recv_buf[i], recv_buf[i + 1]);
    }
    {
        [[maybe_unused]] auto _clear = std::move(recv_buf);
    }

    // Deduplicate edges
    SortAndRemoveDuplicates(local_edges);

    // Set original edge list to new edge list
    std::swap(local_edges, edge_list);
}

std::vector<SInt> ComputeBalancedVertexDistribution(const SInt n, MPI_Comm comm) {
    PEID size;
    MPI_Comm_size(comm, &size);

    const SInt num_vertices_per_pe = n / size;
    const SInt remaining_vertices  = n % size;

    std::vector<SInt> distribution(size + 1);
    for (PEID pe = 0; pe < size; ++pe) {
        distribution[pe] = num_vertices_per_pe + (static_cast<SInt>(pe) < remaining_vertices);
    }
    std::exclusive_scan(distribution.begin(), distribution.end(), distribution.begin(), SInt{0});
    distribution.back() = n;

    return distribution;
}

std::vector<SInt> RoundRobinRemapping(Edgelist& edges, const SInt n, MPI_Comm comm) {
    PEID size;
    MPI_Comm_size(comm, &size);

    std::vector<SInt> distribution = ComputeBalancedVertexDistribution(n, comm);

    auto compute_owner = [&](const SInt id) {
        return id % size;
    };
    auto compute_remap = [&](const SInt id) {
        return distribution[compute_owner(id)] + id / size;
    };
    for (auto& [u, v]: edges) {
        u = compute_remap(u);
        v = compute_remap(v);
    }
    return distribution;
}

class Distribution {
public:
    Distribution(std::vector<SInt> distribution) : distribution_(std::move(distribution)) {}
    PEID compute_owner(SInt v) const {
        auto it = std::upper_bound(distribution_.begin(), distribution_.end(), v);
        assert(it != distribution_.end());
        PEID pe = static_cast<int>(it - distribution_.begin() - 1);
        return pe;
    }

    SInt compute_local_index(SInt v, PEID rank) const {
        assert(distribution_[rank] <= v && v < distribution_[rank + 1]);
        return v - distribution_[rank];
    }

    auto const& get_underlying_dist() const {
        return distribution_;
    }

    SInt local_count(PEID rank) const {
        assert(rank + 1 < static_cast<int>(distribution_.size()));
        return distribution_[rank + 1] - distribution_[rank];
    }

    VertexRange get_vertex_range(PEID rank) const {
        assert(rank + 1 < static_cast<int>(distribution_.size()));
        return VertexRange{distribution_[rank], distribution_[rank + 1]};
    }

private:
    std::vector<SInt> distribution_;
};

std::vector<SInt>
ComputeBalancedEdgeDistribution(Edgelist const& edges, const std::vector<SInt>& vertex_distribution, MPI_Comm comm) {
    PEID size;
    PEID rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    const SInt   n = vertex_distribution.back();
    Distribution input_dist(vertex_distribution);

    // Count out-degree of each source vertex from the local edge list
    std::unordered_map<SInt, SInt> partial_degree;
    for (const auto& [u, v]: edges) {
        ++partial_degree[u];
    }

    // Pack (vertex, degree) pairs into a send buffer, grouped by owning PE
    std::vector<int> send_counts(size);
    for (const auto& [u, degree]: partial_degree) {
        send_counts[input_dist.compute_owner(u)] += 2; // 2 SInts per entry: vertex + degree
    }

    std::vector<int> send_displs(size);
    std::exclusive_scan(send_counts.begin(), send_counts.end(), send_displs.begin(), 0);
    auto              write_offsets = send_displs;
    std::vector<SInt> sendbuf(send_displs.back() + send_counts.back());

    for (const auto& [u, degree]: partial_degree) {
        auto& pos        = write_offsets[input_dist.compute_owner(u)];
        sendbuf[pos]     = u;
        sendbuf[pos + 1] = degree;
        pos += 2;
    }

    // Send partial degrees to the PE that owns each vertex
    std::vector<int> recv_counts(size);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm);
    std::vector<int> recv_displs(size);
    std::exclusive_scan(recv_counts.begin(), recv_counts.end(), recv_displs.begin(), 0);

    std::vector<SInt> recvbuf(recv_counts.back() + recv_displs.back());
    MPI_Alltoallv(
        sendbuf.data(), send_counts.data(), send_displs.data(), KAGEN_MPI_SINT, recvbuf.data(), recv_counts.data(),
        recv_displs.data(), KAGEN_MPI_SINT, comm);
    partial_degree.clear();

    // Aggregate partial degrees from all PEs into the total degree for each owned vertex
    SInt              total_degree = 0;
    const SInt        local_n      = input_dist.local_count(rank);
    std::vector<SInt> local_degree(local_n, 0);
    for (std::size_t i = 0; i < recvbuf.size(); i += 2) {
        const SInt u      = recvbuf[i];
        const SInt degree = recvbuf[i + 1];
        total_degree += degree;
        local_degree[input_dist.compute_local_index(u, rank)] += degree;
    }

    // Compute prefix sum of degrees across PEs and total edge count m
    SInt prefix_sum = 0;
    MPI_Exscan(&total_degree, &prefix_sum, 1, KAGEN_MPI_SINT, MPI_SUM, comm);
    if (rank == 0) {
        prefix_sum = 0;
    }
    SInt m = total_degree;
    MPI_Allreduce(MPI_IN_PLACE, &m, 1, KAGEN_MPI_SINT, MPI_SUM, comm);

    if (m == 0) {
        std::vector<SInt> dist(size + 1, n);
        dist[0] = 0;
        return dist;
    }

    // Find vertex breakpoints that split the edge set into ~equal-sized buckets (one per PE).
    // Walk through owned vertices in global order; whenever the cumulative degree crosses
    // a bucket boundary, record that vertex as a breakpoint in the new distribution.
    //
    // A vertex is only ever recorded as a breakpoint *once*, even if its own degree spans
    // multiple bucket boundaries (i.e. a single hub vertex's degree exceeds bucket_size). This
    // means a hub-owning PE simply absorbs the overflow instead of the naive approach of
    // emitting the same breakpoint repeatedly, which would make Distribution hand several
    // consecutive PEs an empty (v, v) range, starving them entirely. The residual imbalance
    // this leaves behind is bounded by O(max single-vertex degree) on at most one PE; use
    // RedistributeEdgesTrueBalance if that bound is not tight enough.
    std::vector<SInt> breakpoints;
    SInt              cur_sum = prefix_sum;
    assert(size > 0); // size = 0 is not possible in MPI code
    SInt bucket_size_remainder = m % size;
    SInt bucket_size           = std::max(
        static_cast<SInt>(1),
        m / size + static_cast<SInt>(bucket_size_remainder != 0)); // all PEs get one extra edge except for the last one
    if (rank == 0) {
        breakpoints.push_back(0);
    }
    for (std::size_t i = 0; i < local_degree.size(); ++i) {
        SInt degree = local_degree[i];
        SInt low    = (cur_sum + bucket_size - 1) / bucket_size; // ceil(cur_sum / bucket_size)
        SInt start  = std::max(static_cast<SInt>(1), low);
        if (start * bucket_size < cur_sum + degree) {
            // Vertex i itself is the one whose degree crosses the boundary: keep it (and its whole degree) in
            // the *current* bucket -- i.e. actually absorb the overflow here -- and start the next bucket at
            // the following vertex. This breakpoint is therefore always strictly greater than any previously
            // recorded one (i + 1 > i), so it can never collide with -- and spuriously empty out -- an
            // adjacent PE's range, even if vertex i alone already exceeds an entire PE's fair share.
            breakpoints.push_back(i + 1 + vertex_distribution[rank]);
        }
        cur_sum += local_degree[i];
    }
    // Gather all locally computed breakpoints into the global edge-balanced distribution
    std::vector<SInt> edge_balanced_distribution;
    {
        std::vector<int> recvcounts(size, 0);
        int              local_len = static_cast<int>(breakpoints.size());
        MPI_Allgather(&local_len, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, comm);
        std::vector<int> displs(size, 0);
        int              total = 0;
        for (int i = 0; i < size; ++i) {
            displs[i] = total;
            total += recvcounts[i];
        }
        edge_balanced_distribution.resize(total);
        MPI_Allgatherv(
            breakpoints.data(), local_len, KAGEN_MPI_SINT, edge_balanced_distribution.data(), recvcounts.data(),
            displs.data(), KAGEN_MPI_SINT, comm);
    }
    // Pad with n if fewer than size+1 breakpoints were found (low-degree graphs)
    for (std::size_t i = edge_balanced_distribution.size(); i < static_cast<SInt>(size + 1); ++i) {
        edge_balanced_distribution.push_back(n);
    }
    return edge_balanced_distribution;
}

VertexRange RedistributeEdgesBalanced(
    Edgelist& source, Edgelist& destination, const SInt n, bool remap_round_robin, MPI_Comm comm) {
    SortAndRemoveDuplicates(source);
    std::vector<SInt> vertex_distribution;
    if (remap_round_robin) {
        vertex_distribution = RoundRobinRemapping(source, n, comm);
    } else {
        vertex_distribution = ComputeBalancedVertexDistribution(n, comm);
    }
    std::vector<SInt> edge_balanced_distribution = ComputeBalancedEdgeDistribution(source, vertex_distribution, comm);
    Distribution      dist(edge_balanced_distribution);

    PEID rank;
    MPI_Comm_rank(comm, &rank);

    VertexRange vertex_range = dist.get_vertex_range(rank);
    RedistributeEdgesByVertexRange(source, vertex_range, comm, true);
    std::swap(source, destination);

    return vertex_range;
}

VertexRange
RedistributeEdges(Edgelist& source, Edgelist& destination, const SInt n, bool remap_round_robin, MPI_Comm comm) {
    SortAndRemoveDuplicates(source);
    std::vector<SInt> distribution;
    if (remap_round_robin) {
        distribution = RoundRobinRemapping(source, n, comm);
    } else {
        distribution = ComputeBalancedVertexDistribution(n, comm);
    }

    PEID rank;
    MPI_Comm_rank(comm, &rank);
    VertexRange vertex_range = {distribution[rank], distribution[rank + 1]};

    RedistributeEdgesByVertexRange(source, vertex_range, comm, true);
    std::swap(source, destination);

    return vertex_range;
}

BoundaryOwnership ComputeBoundaryOwnership(
    const SInt first_tail, const SInt last_tail, const bool has_local_edges, const SInt n, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Determine which of this PE's boundary vertices are shared with an adjacent PE, and derive a complete
    // (gap-free) partition of [0, n) for this PE's vertex ownership. Both need to see past an empty neighbor to
    // find the next PE that actually holds edges, so both are computed from one cheap MPI_Allgather of every
    // PE's (first tail, last tail) rather than a pairwise exchange.
    const SInt my_first_tail = has_local_edges ? first_tail : n;
    const SInt my_last_tail  = has_local_edges ? last_tail : n;

    std::vector<SInt> all_first_tails(static_cast<std::size_t>(size));
    std::vector<SInt> all_last_tails(static_cast<std::size_t>(size));
    MPI_Allgather(&my_first_tail, 1, KAGEN_MPI_SINT, all_first_tails.data(), 1, KAGEN_MPI_SINT, comm);
    MPI_Allgather(&my_last_tail, 1, KAGEN_MPI_SINT, all_last_tails.data(), 1, KAGEN_MPI_SINT, comm);

    const PEID left_neighbor  = rank > 0 ? rank - 1 : MPI_PROC_NULL;
    const PEID right_neighbor = rank + 1 < size ? rank + 1 : MPI_PROC_NULL;

    BoundaryOwnership result;
    result.is_left_partial  = has_local_edges && left_neighbor != MPI_PROC_NULL
                              && all_last_tails[static_cast<std::size_t>(left_neighbor)] == my_first_tail;
    result.is_right_partial = has_local_edges && right_neighbor != MPI_PROC_NULL
                              && all_first_tails[static_cast<std::size_t>(right_neighbor)] == my_last_tail;

    // Complete, gap-free partition of [0, n): walk all PEs in rank order. A PE with local edges owns through
    // (its own last tail + 1); this resolves a shared boundary vertex to the lower-rank PE (whichever PE's last
    // tail equals it), and lets an empty PE -- or a run of several -- defer the running boundary forward
    // without owning anything new, which is exactly what's needed to absorb a stretch of isolated (degree-0)
    // vertices into whichever PE precedes it. The last PE always closes the partition at n, absorbing any
    // trailing isolated vertices (or, in the degenerate case where no PE has any edges at all, the entire
    // vertex set).
    SInt running_boundary = 0;
    SInt my_range_start   = 0;
    SInt my_range_end     = 0;
    for (PEID pe = 0; pe < size; ++pe) {
        if (pe == rank) {
            my_range_start = running_boundary;
        }
        if (all_last_tails[static_cast<std::size_t>(pe)] != n) {
            running_boundary = all_last_tails[static_cast<std::size_t>(pe)] + 1;
        }
        if (pe == size - 1) {
            running_boundary = n;
        }
        if (pe == rank) {
            my_range_end = running_boundary;
        }
    }
    result.vertex_range = {my_range_start, my_range_end};

    // fully_owned_vertex_range is the strict subset of vertex_range excluding a shared boundary vertex from
    // *both* sides sharing it (neither alone holds its whole adjacency) -- a cheap trim of the already-computed
    // complete range, no new communication needed. Clamp to an empty range at my_range_start if trimming both
    // ends would cross (this PE's only local vertices were entirely boundary-shared ones).
    const SInt owned_start = my_range_start + (result.is_left_partial ? 1 : 0);
    const SInt owned_end   = my_range_end - (result.is_right_partial ? 1 : 0);
    result.fully_owned_vertex_range =
        (owned_start <= owned_end) ? VertexRange{owned_start, owned_end} : VertexRange{owned_start, owned_start};

    bool any_split = result.is_left_partial || result.is_right_partial;
    MPI_Allreduce(MPI_IN_PLACE, &any_split, 1, MPI_C_BOOL, MPI_LOR, comm);
    result.has_split_vertices = any_split;

    return result;
}

EdgeBalancedDistribution ComputeEdgeBalancedBoundaries(const Edgelist& edges, const SInt n, MPI_Comm comm) {
    const bool has_local_edges = !edges.empty();
    const SInt first_tail      = has_local_edges ? edges.front().first : n;
    const SInt last_tail       = has_local_edges ? edges.back().first : n;

    const BoundaryOwnership bo = ComputeBoundaryOwnership(first_tail, last_tail, has_local_edges, n, comm);

    EdgeBalancedDistribution result;
    result.vertex_range             = bo.vertex_range;
    result.fully_owned_vertex_range = bo.fully_owned_vertex_range;
    result.has_split_vertices       = bo.has_split_vertices;

    if (has_local_edges) {
        auto tail_less = [](const std::pair<SInt, SInt>& edge, SInt v) {
            return edge.first < v;
        };
        auto tail_greater = [](SInt v, const std::pair<SInt, SInt>& edge) {
            return v < edge.first;
        };
        if (bo.is_left_partial) {
            const auto first_run_end   = std::upper_bound(edges.begin(), edges.end(), first_tail, tail_greater);
            const auto count           = static_cast<SInt>(std::distance(edges.begin(), first_run_end));
            result.left_partial_vertex = SplitVertexInfo{first_tail, 0, count};
        }
        if (bo.is_right_partial) {
            const auto last_run_begin   = std::lower_bound(edges.begin(), edges.end(), last_tail, tail_less);
            const auto offset           = static_cast<SInt>(std::distance(edges.begin(), last_run_begin));
            result.right_partial_vertex = SplitVertexInfo{last_tail, offset, static_cast<SInt>(edges.size()) - offset};
        }
    }

    return result;
}

EdgeBalancedDistribution RedistributeEdgesTrueBalance(
    Edgelist& source, Edgelist& destination, const SInt n, bool remap_round_robin, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    SortAndRemoveDuplicates(source);

    // A rough vertex distribution used only to route lightweight per-vertex metadata (never edge payloads) to a
    // deterministic "owner" PE for each vertex; this does not influence the final edge-balanced assignment.
    std::vector<SInt> vertex_distribution;
    if (remap_round_robin) {
        vertex_distribution = RoundRobinRemapping(source, n, comm);
        // Remapping changes vertex IDs, so the (tail, head) sort order computed above no longer holds.
        SortAndRemoveDuplicates(source);
    } else {
        vertex_distribution = ComputeBalancedVertexDistribution(n, comm);
    }
    Distribution vertex_owner_dist(vertex_distribution);

    // Group local edges into contiguous per-vertex runs (source is sorted by (tail, head), so all of one
    // vertex's local edges are contiguous).
    struct Run {
        SInt vertex;
        SInt start;
        SInt count;
    };
    std::vector<Run> runs;
    for (std::size_t i = 0; i < source.size();) {
        const SInt  v = source[i].first;
        std::size_t j = i;
        while (j < source.size() && source[j].first == v) {
            ++j;
        }
        runs.push_back({v, static_cast<SInt>(i), static_cast<SInt>(j - i)});
        i = j;
    }

    // Send (vertex, local degree) pairs to each vertex's nominal owner, remembering -- for each local run -- the
    // exact slot its entry ends up in, so a later reply (see below) can be matched back to it in O(1).
    std::vector<int> send_counts(size, 0);
    for (const auto& run: runs) {
        send_counts[vertex_owner_dist.compute_owner(run.vertex)] += 2; // 2 SInts per entry: vertex + count
    }
    std::vector<int> send_displs(size);
    std::exclusive_scan(send_counts.begin(), send_counts.end(), send_displs.begin(), 0);
    std::vector<SInt>        send_buf(send_displs.back() + send_counts.back());
    std::vector<std::size_t> run_reply_slot(runs.size());
    {
        auto write_offsets = send_displs;
        for (std::size_t r = 0; r < runs.size(); ++r) {
            const auto& run   = runs[r];
            const PEID  owner = vertex_owner_dist.compute_owner(run.vertex);
            auto&       pos   = write_offsets[owner];
            send_buf[pos]     = run.vertex;
            send_buf[pos + 1] = run.count;
            run_reply_slot[r] = static_cast<std::size_t>(pos) / 2;
            pos += 2;
        }
    }

    std::vector<int> recv_counts(size);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm);
    std::vector<int> recv_displs(size);
    std::exclusive_scan(recv_counts.begin(), recv_counts.end(), recv_displs.begin(), 0);
    std::vector<SInt> recv_buf(recv_displs.back() + recv_counts.back());
    MPI_Alltoallv(
        send_buf.data(), send_counts.data(), send_displs.data(), KAGEN_MPI_SINT, recv_buf.data(), recv_counts.data(),
        recv_displs.data(), KAGEN_MPI_SINT, comm);

    // Group received (vertex, count) contributions by local vertex index. Processing recv_buf in ascending
    // sender-PE block order means each vertex's contributions end up ordered by ascending sender PE (each sender
    // contributes at most one entry per vertex, since runs groups a sender's own edges by vertex already).
    const SInt                                      local_n = vertex_owner_dist.local_count(rank);
    std::vector<SInt>                               local_degree(local_n, 0);
    std::vector<std::vector<std::pair<PEID, SInt>>> contributions(local_n);
    for (PEID sender = 0; sender < size; ++sender) {
        const int begin = recv_displs[sender];
        const int end   = recv_displs[sender] + recv_counts[sender];
        for (int pos = begin; pos < end; pos += 2) {
            const SInt v      = recv_buf[pos];
            const SInt degree = recv_buf[pos + 1];
            const SInt li     = vertex_owner_dist.compute_local_index(v, rank);
            contributions[li].emplace_back(sender, degree);
            local_degree[li] += degree;
        }
    }

    SInt total_degree = std::accumulate(local_degree.begin(), local_degree.end(), SInt{0});
    SInt prefix_sum   = 0;
    MPI_Exscan(&total_degree, &prefix_sum, 1, KAGEN_MPI_SINT, MPI_SUM, comm);
    if (rank == 0) {
        prefix_sum = 0;
    }
    SInt m = total_degree;
    MPI_Allreduce(MPI_IN_PLACE, &m, 1, KAGEN_MPI_SINT, MPI_SUM, comm);

    if (m == 0) {
        destination.clear();
        EdgeBalancedDistribution result;
        result.vertex_range             = rank == 0 ? VertexRange{0, n} : VertexRange{n, n};
        result.fully_owned_vertex_range = result.vertex_range;
        result.has_split_vertices       = false;
        return result;
    }

    // Closed-form balanced distribution of global edge ranks [0, m) across PEs: PE p owns
    // [edge_distribution[p], edge_distribution[p + 1]).
    const std::vector<SInt> edge_distribution = ComputeBalancedVertexDistribution(m, comm);
    Distribution            edge_owner_dist(edge_distribution);

    // For each local vertex (in ascending order), determine whether its degree crosses more than one PE's fair
    // edge share ("split") and, if not, each contributing sender's base global edge rank -- the global rank of
    // that sender's first local edge of this vertex -- to reply directly to that sender. This is lightweight
    // metadata only (one SInt per (vertex, sender) pair); no edge payload is sent here.
    //
    // A non-split vertex's edges all go to the same single target PE by construction (its whole [cur_sum,
    // cur_sum + degree) range falls in one bucket), so any cross-PE duplicate of one of its edges still ends up
    // co-located and is caught by the final SortAndRemoveDuplicates(destination) below, exactly as for the
    // vertex-atomic modes. That guarantee does *not* hold for a split vertex purely by sender-then-local-offset
    // position (two duplicate copies held by different source PEs could straddle the split and never meet) --
    // see the escalation path further down, which resolves split vertices by actual (deduplicated) edge value
    // instead of by position.
    constexpr SInt    kEscalate = std::numeric_limits<SInt>::max(); // never a valid base global rank ([0, m))
    std::vector<bool> is_split(local_n, false);
    std::vector<SInt> cur_sum_before(local_n);
    std::vector<std::vector<SInt>> reply_send_bufs(size);
    {
        SInt cur_sum = 0;
        for (SInt li = 0; li < local_n; ++li) {
            cur_sum_before[li] = cur_sum;
            const SInt degree  = local_degree[li];
            if (degree > 0) {
                const PEID first_owner = edge_owner_dist.compute_owner(prefix_sum + cur_sum);
                const PEID last_owner  = edge_owner_dist.compute_owner(prefix_sum + cur_sum + degree - 1);
                is_split[li]           = first_owner != last_owner;
            }
            if (is_split[li]) {
                for (const auto& [sender, degree_ignored]: contributions[li]) {
                    reply_send_bufs[sender].push_back(kEscalate);
                }
            } else {
                SInt within_vertex_offset = 0;
                for (const auto& [sender, sender_degree]: contributions[li]) {
                    reply_send_bufs[sender].push_back(prefix_sum + cur_sum + within_vertex_offset);
                    within_vertex_offset += sender_degree;
                }
            }
            cur_sum += degree;
        }
    }

    std::vector<int> reply_send_counts(size);
    for (PEID pe = 0; pe < size; ++pe) {
        reply_send_counts[pe] = static_cast<int>(reply_send_bufs[pe].size());
    }
    std::vector<int> reply_send_displs(size);
    std::exclusive_scan(reply_send_counts.begin(), reply_send_counts.end(), reply_send_displs.begin(), 0);
    std::vector<SInt> reply_send_buf(reply_send_displs.back() + reply_send_counts.back());
    for (PEID pe = 0; pe < size; ++pe) {
        std::copy(
            reply_send_bufs[pe].begin(), reply_send_bufs[pe].end(), reply_send_buf.begin() + reply_send_displs[pe]);
    }

    // The reply traverses the same communication graph as the initial metadata exchange, but in reverse: this
    // PE's original send_counts/send_displs (halved, since the reply carries 1 SInt per original 2-SInt entry)
    // become its receive counts/displs here.
    std::vector<int> reply_recv_counts(size);
    std::vector<int> reply_recv_displs(size);
    for (PEID pe = 0; pe < size; ++pe) {
        reply_recv_counts[pe] = send_counts[pe] / 2;
        reply_recv_displs[pe] = send_displs[pe] / 2;
    }
    std::vector<SInt> reply_recv_buf(reply_recv_displs.back() + reply_recv_counts.back());
    MPI_Alltoallv(
        reply_send_buf.data(), reply_send_counts.data(), reply_send_displs.data(), KAGEN_MPI_SINT,
        reply_recv_buf.data(), reply_recv_counts.data(), reply_recv_displs.data(), KAGEN_MPI_SINT, comm);

    // Route every local edge of a non-split run directly to its final target PE using the base global rank of
    // its run plus its offset within the run. A split run's actual (vertex, head) values are instead escalated
    // to the vertex's owner, which can resolve exact (deduplicated) positions once it has seen every
    // contributing sender's data -- something no single sender can determine on its own.
    std::unordered_map<PEID, std::vector<SInt>> send_buffers;
    std::unordered_map<PEID, std::vector<SInt>> escalate_buffers;
    for (std::size_t r = 0; r < runs.size(); ++r) {
        const auto& run   = runs[r];
        const SInt  reply = reply_recv_buf[run_reply_slot[r]];
        if (reply == kEscalate) {
            const PEID         owner = vertex_owner_dist.compute_owner(run.vertex);
            std::vector<SInt>& buf   = escalate_buffers[owner];
            for (SInt k = 0; k < run.count; ++k) {
                buf.push_back(run.vertex);
                buf.push_back(source[run.start + k].second);
            }
        } else {
            const SInt base_rank = reply;
            for (SInt k = 0; k < run.count; ++k) {
                const SInt         global_rank = base_rank + k;
                const PEID         target      = edge_owner_dist.compute_owner(global_rank);
                const auto&        edge        = source[run.start + k];
                std::vector<SInt>& buf         = send_buffers[target];
                buf.push_back(edge.first);
                buf.push_back(edge.second);
            }
        }
    }

    // Resolve escalated split vertices: gather every contributing sender's actual edges of that vertex (bounded
    // by the split vertices' own total degree, not the whole graph), deduplicate by head value, and forward each
    // remaining edge directly to its true target PE. Because equal head values are adjacent once sorted, this
    // also means any surviving cross-PE duplicate of a split vertex's edge can only ever straddle one boundary
    // between two adjacent PEs -- the same invariant the vertex-atomic modes get for free from single ownership.
    auto escalated = ExchangeMessageBuffers(std::move(escalate_buffers), KAGEN_MPI_SINT, comm);
    {
        std::vector<std::vector<SInt>> split_heads(local_n);
        for (std::size_t i = 0; i < escalated.size(); i += 2) {
            const SInt v  = escalated[i];
            const SInt h  = escalated[i + 1];
            const SInt li = vertex_owner_dist.compute_local_index(v, rank);
            split_heads[li].push_back(h);
        }
        for (SInt li = 0; li < local_n; ++li) {
            if (split_heads[li].empty()) {
                continue;
            }
            auto& heads = split_heads[li];
            std::sort(heads.begin(), heads.end());
            heads.erase(std::unique(heads.begin(), heads.end()), heads.end());

            const SInt vertex = li + vertex_distribution[rank];
            for (std::size_t k = 0; k < heads.size(); ++k) {
                const SInt         global_rank = prefix_sum + cur_sum_before[li] + static_cast<SInt>(k);
                const PEID         target      = edge_owner_dist.compute_owner(global_rank);
                std::vector<SInt>& buf         = send_buffers[target];
                buf.push_back(vertex);
                buf.push_back(heads[k]);
            }
        }
    }

    auto recv_edges = ExchangeMessageBuffers(std::move(send_buffers), KAGEN_MPI_SINT, comm);
    destination.clear();
    destination.reserve(recv_edges.size() / 2);
    for (std::size_t i = 0; i < recv_edges.size(); i += 2) {
        destination.emplace_back(recv_edges[i], recv_edges[i + 1]);
    }
    SortAndRemoveDuplicates(destination);

    // Stage 2: fix residual drift from pre-dedup bucket boundaries. Everything above placed bucket boundaries
    // using *pre-deduplication* degree estimates (m, prefix_sum, edge_owner_dist), so each PE's now-exact local
    // edge count (post dedup) can drift from its +/-1 target by however many duplicate edges were counted in
    // upstream of it -- bounded by total duplicate volume, not by any single vertex's degree (splitting a
    // dominant hub was already handled correctly above, without ever concentrating its whole edge set on one
    // PE). Since every PE now knows its *exact* local count, correcting this is a second, much smaller, purely
    // rank-based rebalance: no vertex grouping, no escalation, no reply protocol needed -- just move each edge
    // to the target PE implied by its exact global rank among the now-deduplicated total.
    {
        SInt exact_local_count = static_cast<SInt>(destination.size());
        SInt exact_prefix      = 0;
        MPI_Exscan(&exact_local_count, &exact_prefix, 1, KAGEN_MPI_SINT, MPI_SUM, comm);
        if (rank == 0) {
            exact_prefix = 0;
        }
        SInt exact_m = exact_local_count;
        MPI_Allreduce(MPI_IN_PLACE, &exact_m, 1, KAGEN_MPI_SINT, MPI_SUM, comm);

        if (exact_m > 0) {
            Distribution exact_edge_owner_dist(ComputeBalancedVertexDistribution(exact_m, comm));

            std::unordered_map<PEID, std::vector<SInt>> rebalance_buffers;
            for (std::size_t i = 0; i < destination.size(); ++i) {
                const SInt global_rank = exact_prefix + static_cast<SInt>(i);
                const PEID target      = exact_edge_owner_dist.compute_owner(global_rank);
                auto&      buf         = rebalance_buffers[target];
                buf.push_back(destination[i].first);
                buf.push_back(destination[i].second);
            }

            auto rebalanced = ExchangeMessageBuffers(std::move(rebalance_buffers), KAGEN_MPI_SINT, comm);
            destination.clear();
            destination.reserve(rebalanced.size() / 2);
            for (std::size_t i = 0; i < rebalanced.size(); i += 2) {
                destination.emplace_back(rebalanced[i], rebalanced[i + 1]);
            }
            std::sort(destination.begin(), destination.end());
        }
    }

    // Derive the gap-free vertex_range, fully_owned_vertex_range, partial_vertices, and has_split_vertices from
    // this PE's now-final sorted local edges. This is the same boundary computation the direct edge-range read
    // path uses, so it lives in a shared helper.
    return ComputeEdgeBalancedBoundaries(destination, n, comm);
}
} // namespace kagen
