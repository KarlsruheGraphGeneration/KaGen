#include "kagen/generators/graph500_generator.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

namespace kagen {
Graph500Generator::Graph500Generator(const PGeneratorConfig& config) : config_(config) {}

void Graph500Generator::Finalize(MPI_Comm comm) {
    const SInt log_n = std::log2(config_.n);
    const SInt n     = 1ull << log_n;
    { // Remove local duplicates
        std::sort(local_edges_.begin(), local_edges_.end());
        auto it = std::unique(local_edges_.begin(), local_edges_.end());
        local_edges_.erase(it, local_edges_.end());
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
    for (const auto& [u, v]: local_edges_) {
        ++send_counts[compute_owner(u)];
    }
    std::vector<int> send_displs(size);
    std::exclusive_scan(send_counts.begin(), send_counts.end(), send_displs.begin(), 0);

    // Remap edges and build send buffer
    std::vector<long long> sendbuf(local_edges_.size());
    std::vector<int>       sendbuf_pos(size);
    for (const auto& [u, v]: local_edges_) {
        const PEID u_owner = compute_owner(u);
        const int  u_prime = compute_remap(u);
        const int  v_prime = compute_remap(v);

        const auto index = send_displs[u_owner] + sendbuf_pos[u_owner];
        sendbuf[index]   = (static_cast<long long>(u_prime) << 32) | static_cast<long long>(v_prime);
        ++sendbuf_pos[u_owner];
    }

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
        PushEdge(u, v);
    }

    FilterDuplicateEdges();
    SetVertexRange(distribution[rank], distribution[rank + 1]);
}
} // namespace kagen
