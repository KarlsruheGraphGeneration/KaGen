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

void SortEdges(EdgeList& edge_list) {
    std::sort(edge_list.begin(), edge_list.end());
}

void AddReverseEdges(EdgeList& edge_list, const VertexRange vertex_range, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const auto ranges = AllgatherVertexRange(vertex_range, comm);

    const SInt edge_list_size_before = edge_list.size();

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

    std::sort(edge_list.begin(), edge_list.end());

    // KaGen sometimes produces duplicate edges
    auto it = std::unique(edge_list.begin(), edge_list.end());
    edge_list.erase(it, edge_list.end());

    const SInt edge_list_size_after = edge_list.size();

    // Generate some statistics
    SInt edge_list_global_size_before = 0;
    MPI_Reduce(&edge_list_size_before, &edge_list_global_size_before, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, comm);
    SInt edge_list_global_size_after = 0;
    MPI_Reduce(&edge_list_size_after, &edge_list_global_size_after, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, comm);
}

void AddReverseEdgesAndRedistribute(EdgeList& edge_list, const VertexRange vertex_range, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const auto ranges = AllgatherVertexRange(vertex_range, comm);
    const auto from   = ranges[rank].first;
    const auto to     = ranges[rank].second;

    // Create new edge arrays
    std::vector<std::tuple<SInt, SInt>>              local_edges;
    std::vector<std::vector<std::tuple<SInt, SInt>>> remote_edges(size);
    for (const auto& [u, v]: edge_list) {
        if (u == v) { // Ignore self loops
            continue;
        }

        if (from <= u && u < to) { // Edge starts from local vertex
            local_edges.emplace_back(u, v);
        } else { // Edge starts from remote vertex
            const PEID owner = FindPEInRange(u, ranges);
            remote_edges[owner].emplace_back(u, v);
        }

        if (from <= v && v < to) {
            local_edges.emplace_back(v, u);
        } else {
            const PEID owner = FindPEInRange(v, ranges);
            remote_edges[owner].emplace_back(v, u);
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

VertexRange RecomputeVertexRanges(const VertexRange vertex_range, MPI_Comm comm) {
    const SInt num_vertices    = vertex_range.second - vertex_range.first;
    SInt       vertices_before = 0;
    MPI_Exscan(&num_vertices, &vertices_before, 1, KAGEN_MPI_SINT, MPI_SUM, comm);
    return std::make_pair(vertices_before, vertices_before + num_vertices);
}
} // namespace kagen
