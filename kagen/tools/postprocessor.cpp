#include "kagen/tools/postprocessor.h"

#include "kagen/tools/utils.h"

#include <mpi.h>

#include <algorithm>
#include <cassert>
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
        send_buf.data(), send_counts.data(), send_displs.data(), KAGEN_MPI_SINT, recv_buf.data(),
        recv_counts.data(), recv_displs.data(), KAGEN_MPI_SINT, comm);
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
    std::exclusive_scan(distribution.begin(), distribution.end(), distribution.begin(), 0);
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

    std::unordered_map<SInt, SInt> partial_degree;
    for (const auto& [u, v]: edges) {
        ++partial_degree[u];
    }

    // Compute send_counts and send_displs
    std::vector<int> send_counts(size);
    for (const auto& [u, degree]: partial_degree) {
        send_counts[input_dist.compute_owner(u)] += 2;
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

    // Exchange send_counts + send_displs
    std::vector<int> recv_counts(size);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm);
    std::vector<int> recv_displs(size);
    std::exclusive_scan(recv_counts.begin(), recv_counts.end(), recv_displs.begin(), 0);

    // Exchange edges
    std::vector<SInt> recvbuf(recv_counts.back() + recv_displs.back());
    MPI_Alltoallv(
        sendbuf.data(), send_counts.data(), send_displs.data(), KAGEN_MPI_SINT, recvbuf.data(), recv_counts.data(),
        recv_displs.data(), KAGEN_MPI_SINT, comm);
    partial_degree.clear();
    SInt              total_degree = 0;
    const SInt        local_n      = input_dist.local_count(rank);
    std::vector<SInt> local_degree(local_n, 0);
    for (std::size_t i = 0; i < recvbuf.size(); i += 2) {
        const SInt u      = recvbuf[i];
        const SInt degree = recvbuf[i + 1];
        total_degree += degree;
        local_degree[input_dist.compute_local_index(u, rank)] += degree;
    }

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
        for (std::size_t j = start * bucket_size; j < cur_sum + degree; j += bucket_size) {
            breakpoints.push_back(i + vertex_distribution[rank]);
        }
        cur_sum += local_degree[i];
    }
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
} // namespace kagen
