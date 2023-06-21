#pragma once

#include <mpi.h>

#include "kagen/definitions.h"

namespace kagen {
void AddReverseEdges(EdgeList& edge_list, VertexRange vertex_range, MPI_Comm comm);

void AddReverseEdgesAndRedistribute(
    EdgeList& edge_list, VertexRange vertex_range, bool skip_self_loops, bool add_reverse_edges, MPI_Comm comm);

template <typename FromEdgeList, typename ToEdgeList>
std::pair<SInt, SInt> RedistributeEdges(
    FromEdgeList& from_edge_list, ToEdgeList& to_edge_list, const bool remove_duplicates_before_redistribution,
    const bool remove_duplicates_after_redistribution, const SInt n, MPI_Comm comm) {
    if (remove_duplicates_before_redistribution) {
        std::sort(from_edge_list.begin(), from_edge_list.end());
        auto it = std::unique(from_edge_list.begin(), from_edge_list.end());
        from_edge_list.erase(it, from_edge_list.end());
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
    for (const auto& [u, v]: from_edge_list) {
        ++send_counts[compute_owner(u)];
    }
    std::vector<int> send_displs(size);
    std::exclusive_scan(send_counts.begin(), send_counts.end(), send_displs.begin(), 0);

    // Remap edges and build send buffer
    std::vector<long long> sendbuf(from_edge_list.size());
    std::vector<int>       sendbuf_pos(size);
    for (const auto& [u, v]: from_edge_list) {
        const PEID u_owner = compute_owner(u);
        const int  u_prime = compute_remap(u);
        const int  v_prime = compute_remap(v);

        const auto index = send_displs[u_owner] + sendbuf_pos[u_owner];
        sendbuf[index]   = (static_cast<long long>(u_prime) << 32) | static_cast<long long>(v_prime);
        ++sendbuf_pos[u_owner];
    }

    // Free the old edge list before allocating the new one
    { [[maybe_unused]] auto tmp = std::move(from_edge_list); }
    to_edge_list.clear();

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
        to_edge_list.emplace_back(u, v);
    }

    if (remove_duplicates_after_redistribution) {
        std::sort(to_edge_list.begin(), to_edge_list.end());
        auto it = std::unique(to_edge_list.begin(), to_edge_list.end());
        to_edge_list.erase(it, to_edge_list.end());
    }

    return {distribution[rank], distribution[rank + 1]};
}
} // namespace kagen
