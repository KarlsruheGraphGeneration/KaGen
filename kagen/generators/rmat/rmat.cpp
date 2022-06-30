#include "kagen/generators/rmat/rmat.h"

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/generators/rmat/generators/select.hpp"
#include "kagen/generators/rmat/rmat_impl.hpp"

#include <numeric>
#include <mpi.h>

namespace kagen {
PGeneratorConfig RMATFactory::NormalizeParameters(PGeneratorConfig config, bool output) const {
    if (config.rmat_a < 0 || config.rmat_b < 0 || config.rmat_b < 0) {
        throw ConfigurationError("probabilities may not be negative");
    }
    if (config.rmat_a + config.rmat_b + config.rmat_c > 1) {
        throw ConfigurationError("sum of probabilities may not be larger than 1");
    }

    const SInt log_n = std::log2(config.n);
    if (log_n > 31) {
        throw ConfigurationError("number of vertices is too large (cannot be larger than 31 bits)");
    }

    if (output && config.n != 1ull << log_n) {
        std::cout << "Warning: generator requires the number of vertices to be a power of two" << std::endl;
        std::cout << "  Changing the number of vertices to " << (1ull << log_n) << std::endl;
    }

    return config;
}

std::unique_ptr<Generator> RMATFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<RMAT>(config, rank, size);
}

RMAT::RMAT(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size) {}

void RMAT::GenerateImpl() {
    using RNG  = rmat::generators::select_t;
    using RMAT = rmat::rmat<false>;

    const SInt seed  = rank_ + 1;
    const SInt log_n = std::log2(config_.n);
    const SInt n     = 1ull << log_n;
    const SInt depth = std::min<SInt>(9, log_n);

    RNG  gen(seed);
    RNG  gen_scramble(seed + 1000);
    RMAT r(gen_scramble, log_n, config_.rmat_a, config_.rmat_b, config_.rmat_c);
    r.init(depth);

    // Generate local edges
    std::vector<std::tuple<int, int>> local_edges;
    local_edges.reserve(config_.m);
    r.get_edges(
        [&](const auto u, const auto v) {
            if (u != v) {
                local_edges.emplace_back(u, v);
                local_edges.emplace_back(v, u);
            }
        },
        config_.m, gen);

    { // Remove local duplicates
        std::sort(local_edges.begin(), local_edges.end());
        auto it = std::unique(local_edges.begin(), local_edges.end());
        local_edges.erase(it, local_edges.end());
    }

    // Remove vertex distribution (round-robin)
    const int num_vertices_per_pe = n / size_;
    const int remaining_vertices  = n % size_;

    std::vector<int> distribution(size_ + 1);
    for (PEID pe = 0; pe < size_; ++pe) {
        distribution[pe] = num_vertices_per_pe + (pe < remaining_vertices);
    }
    std::partial_sum(distribution.begin(), distribution.end(), distribution.begin() + 1);
    distribution.front() = 0;

    /*
    std::cout << "distribution: ";
    for (const auto &v : distribution) { std::cout << v << " "; }
    std::cout << std::endl;
    std::cout << "n: " << n << std::endl;
    */

    // Find number of edges for each PE
    auto compute_owner = [&](const int id) {
        return id % size_;
    };
    auto compute_remap = [&](const int id) {
        return distribution[compute_owner(id)] + id / size_;
    };

    // Compute send_counts and send_displs
    std::vector<int> send_counts(size_);
    for (const auto& [u, v]: local_edges) {
        ++send_counts[compute_owner(u)];
    }
    std::vector<int> send_displs(size_);
    std::exclusive_scan(send_counts.begin(), send_counts.end(), send_displs.begin(), 0);

    // Remap edges and build send buffer
    std::vector<long long> sendbuf(local_edges.size());
    std::vector<int>       sendbuf_pos(size_);
    for (const auto& [u, v]: local_edges) {
        const PEID u_owner = compute_owner(u);
        const int  u_prime = compute_remap(u);
        const int  v_prime = compute_remap(v);

        // std::cout << "Remap " << u << " --> " << v << " to " << u_prime << " --> " << v_prime << " on PE " << u_owner
        //<< std::endl;
        const auto index = send_displs[u_owner] + sendbuf_pos[u_owner];
        sendbuf[index]   = (static_cast<long long>(u_prime) << 32) | static_cast<long long>(v_prime);
        ++sendbuf_pos[u_owner];
    }

    // Exchange send_counts + send_displs
    std::vector<int> recv_counts(size_);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD); // @todo comm
    std::vector<int> recv_displs(size_);
    std::exclusive_scan(recv_counts.begin(), recv_counts.end(), recv_displs.begin(), 0);

    // Exchange edges
    std::vector<long long> recvbuf(recv_counts.back() + recv_displs.back());
    MPI_Alltoallv(
        sendbuf.data(), send_counts.data(), send_displs.data(), MPI_LONG_LONG, recvbuf.data(), recv_counts.data(),
        recv_displs.data(), MPI_LONG_LONG, MPI_COMM_WORLD);

    for (const auto& edge: recvbuf) {
        const int u = edge >> 32;
        const int v = edge & 0xFFFFFFFF;
        PushEdge(u, v);
    }

    FilterDuplicateEdges();
    SetVertexRange(distribution[rank_], distribution[rank_ + 1]);
}
} // namespace kagen
