#include "kagen/tools/postprocessor.h"

#include <mpi.h>

#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace kagen {
namespace {
inline PEID FindPEInRange(const SInt node, const std::vector<std::pair<SInt, SInt>>& ranges) {
    for (std::size_t i = 0; i < ranges.size(); ++i) {
        const auto& [local_from, local_to] = ranges[i];

        if (local_from <= node && node < local_to) {
            return i;
        }
    }

    return -1;
}

std::vector<VertexRange> AllgatherVertexRange(const VertexRange vertex_range, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::vector<VertexRange> ranges(static_cast<std::size_t>(size));
    ranges[static_cast<std::size_t>(rank)] = vertex_range;
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, ranges.data(), sizeof(VertexRange), MPI_BYTE, comm);

    return ranges;
}
} // namespace

void AddNonlocalReverseEdges(Edgelist& edge_list, const VertexRange vertex_range, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const auto ranges = AllgatherVertexRange(vertex_range, comm);

    const auto& [local_from, local_to] = ranges[rank];
    std::unordered_map<PEID, std::vector<SInt>> message_buffers;

    // Each PE gets the edges that we have to that PE
    for (const auto& [tail, head]: edge_list) {
        if ((tail >= local_from && tail < local_to) && (head < local_from || head >= local_to)) {
            const SInt pe = static_cast<SInt>(FindPEInRange(head, ranges));
            message_buffers[pe].emplace_back(tail);
            message_buffers[pe].emplace_back(head);
        }
    }

    std::vector<SInt> send_buf;
    std::vector<SInt> recv_buf;
    std::vector<int>  send_counts(size);
    std::vector<int>  recv_counts(size);
    std::vector<int>  send_displs(size);
    std::vector<int>  recv_displs(size);
    for (size_t i = 0; i < send_counts.size(); ++i) {
        send_counts[i] = message_buffers[i].size();
    }

    std::exclusive_scan(send_counts.begin(), send_counts.end(), send_displs.begin(), 0);
    const std::size_t total_send_count = send_displs.back() + send_counts.back();
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm);
    std::exclusive_scan(recv_counts.begin(), recv_counts.end(), recv_displs.begin(), 0);
    const std::size_t total_recv_count = recv_displs.back() + recv_counts.back();
    send_buf.reserve(total_send_count);
    for (size_t i = 0; i < send_counts.size(); ++i) {
        for (const auto& elem: message_buffers[i]) {
            send_buf.push_back(elem);
        }
        message_buffers[i].clear();
        message_buffers[i].resize(0);
    }
    recv_buf.resize(total_recv_count);
    MPI_Alltoallv(
        send_buf.data(), send_counts.data(), send_displs.data(), MPI_UINT64_T, recv_buf.data(), recv_counts.data(),
        recv_displs.data(), MPI_UINT64_T, comm);
    send_buf.clear();
    send_buf.resize(0);
    for (std::size_t i = 0; i < recv_buf.size(); i += 2) {
        edge_list.emplace_back(recv_buf[i + 1], recv_buf[i]);
    }
    recv_buf.clear();
    recv_buf.resize(0);

    // KaGen sometimes produces duplicate edges
    std::sort(edge_list.begin(), edge_list.end());
    auto it = std::unique(edge_list.begin(), edge_list.end());
    edge_list.erase(it, edge_list.end());
}

void RedistributeEdgesByVertexRange(Edgelist& edge_list, const VertexRange vertex_range, MPI_Comm comm) {
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
            const PEID owner = FindPEInRange(u, ranges);
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
        { [[maybe_unused]] auto _clear = std::move(remote_edges[i]); }
    }
    recv_buf.resize(total_recv_count);
    local_edges.reserve(local_edges.size() + total_recv_count / 2);

    MPI_Alltoallv(
        send_buf.data(), send_counts.data(), send_displs.data(), MPI_UNSIGNED_LONG_LONG, recv_buf.data(),
        recv_counts.data(), recv_displs.data(), MPI_UNSIGNED_LONG_LONG, comm);
    { [[maybe_unused]] auto _clear = std::move(send_buf); }

    for (std::size_t i = 0; i < recv_buf.size(); i += 2) {
        local_edges.emplace_back(recv_buf[i], recv_buf[i + 1]);
    }
    { [[maybe_unused]] auto _clear = std::move(recv_buf); }

    // Deduplicate edges
    std::sort(local_edges.begin(), local_edges.end());
    auto it = std::unique(local_edges.begin(), local_edges.end());
    local_edges.erase(it, local_edges.end());

    // Set original edge list to new edge list
    std::swap(local_edges, edge_list);
}

VertexRange RedistributeEdgesRoundRobin(Edgelist32& source, Edgelist& destination, const SInt n, MPI_Comm comm) {
    {
        std::sort(source.begin(), source.end());
        auto it = std::unique(source.begin(), source.end());
        source.erase(it, source.end());
    }

    PEID size;
    PEID rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // Remove vertex distribution (round-robin)
    const int num_vertices_per_pe = n / size;
    const int remaining_vertices  = n % size;

    std::vector<int> distribution(size + 1);
    for (PEID pe = 0; pe < size; ++pe) {
        distribution[pe] = num_vertices_per_pe + (pe < remaining_vertices);
    }
    std::exclusive_scan(distribution.begin(), distribution.end(), distribution.begin(), 0);
    distribution.back() = n;

    // Find number of edges for each PE
    auto compute_owner = [&](const int id) {
        return id % size;
    };
    auto compute_remap = [&](const int id) {
        return distribution[compute_owner(id)] + id / size;
    };

    // Compute send_counts and send_displs
    std::vector<int> send_counts(size);
    for (const auto& [u, v]: source) {
        ++send_counts[compute_owner(u)];
    }
    std::vector<int> send_displs(size);
    std::exclusive_scan(send_counts.begin(), send_counts.end(), send_displs.begin(), 0);

    // Remap edges and build send buffer
    std::vector<long long> sendbuf(source.size());
    std::vector<int>       sendbuf_pos(size);
    for (const auto& [u, v]: source) {
        const PEID u_owner = compute_owner(u);
        const int  u_prime = compute_remap(u);
        const int  v_prime = compute_remap(v);

        const auto index = send_displs[u_owner] + sendbuf_pos[u_owner];
        sendbuf[index]   = (static_cast<long long>(u_prime) << 32) | static_cast<long long>(v_prime);
        ++sendbuf_pos[u_owner];
    }

    // Free the old edge list before allocating the new one
    { [[maybe_unused]] auto tmp = std::move(source); }
    destination.clear();

    // Exchange send_counts + send_displs
    std::vector<int> recv_counts(size);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm);
    std::vector<int> recv_displs(size);
    std::exclusive_scan(recv_counts.begin(), recv_counts.end(), recv_displs.begin(), 0);

    // Exchange edges
    std::vector<long long> recvbuf(recv_counts.back() + recv_displs.back());
    MPI_Alltoallv(
        sendbuf.data(), send_counts.data(), send_displs.data(), MPI_LONG_LONG, recvbuf.data(), recv_counts.data(),
        recv_displs.data(), MPI_LONG_LONG, comm);

    for (const auto& edge: recvbuf) {
        const int u = edge >> 32;
        const int v = edge & 0xFFFFFFFF;
        destination.emplace_back(u, v);
    }

    {
        std::sort(destination.begin(), destination.end());
        auto it = std::unique(destination.begin(), destination.end());
        destination.erase(it, destination.end());
    }

    return {distribution[rank], distribution[rank + 1]};
}
} // namespace kagen
