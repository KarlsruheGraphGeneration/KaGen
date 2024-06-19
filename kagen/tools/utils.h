#pragma once

#include "kagen/kagen.h"

#include <mpi.h>

#include <limits>
#include <numeric>

namespace kagen {
inline std::pair<SInt, SInt> ComputeRange(const SInt n, const PEID size, const PEID rank) {
    const SInt chunk = n / size;
    const SInt rem   = n % size;
    const SInt from  = rank * chunk + std::min<SInt>(rank, rem);
    const SInt to    = std::min<SInt>(from + ((static_cast<SInt>(rank) < rem) ? chunk + 1 : chunk), n);
    return {from, to};
}

inline SInt FindNumberOfVerticesInEdgelist(const Edgelist& edges, MPI_Comm comm) {
    SInt n = 0;
    for (const auto& [u, v]: edges) {
        n = std::max(n, std::max(u, v));
    }
    MPI_Allreduce(MPI_IN_PLACE, &n, 1, KAGEN_MPI_SINT, MPI_MAX, comm);
    return n + 1;
}

inline PEID GetCommRank(MPI_Comm comm) {
    PEID rank;
    MPI_Comm_rank(comm, &rank);
    return rank;
}

inline PEID GetCommSize(MPI_Comm comm) {
    PEID size;
    MPI_Comm_size(comm, &size);
    return size;
}

template <typename T>
T FloorLog2(const T arg) {
    constexpr std::size_t arg_width = std::numeric_limits<T>::digits;
    static_assert(
        arg_width == std::numeric_limits<unsigned long>::digits
            || arg_width == std::numeric_limits<unsigned int>::digits,
        "unsupported data type width");

    T log2 = static_cast<T>(arg_width);
    if constexpr (arg_width == std::numeric_limits<unsigned int>::digits) { // 32 bit
        log2 -= __builtin_clz(arg);
    } else { // 64 bit
        log2 -= __builtin_clzl(arg);
    }

    return log2 - 1;
}

inline PEID FindPEInRange(const SInt node, const std::vector<std::pair<SInt, SInt>>& ranges) {
    for (std::size_t i = 0; i < ranges.size(); ++i) {
        const auto& [local_from, local_to] = ranges[i];

        if (local_from <= node && node < local_to) {
            return i;
        }
    }

    return -1;
}

inline std::vector<VertexRange> AllgatherVertexRange(const VertexRange vertex_range, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::vector<VertexRange> ranges(static_cast<std::size_t>(size));
    ranges[static_cast<std::size_t>(rank)] = vertex_range;
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, ranges.data(), sizeof(VertexRange), MPI_BYTE, comm);

    return ranges;
}

template <typename T>
std::vector<T> ExchangeMessageBuffers(
    std::unordered_map<PEID, std::vector<T>> message_buffers, MPI_Datatype mpi_datatype, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    std::vector<T>   send_buf;
    std::vector<T>   recv_buf;
    std::vector<int> send_counts(size);
    std::vector<int> recv_counts(size);
    std::vector<int> send_displs(size);
    std::vector<int> recv_displs(size);
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
        send_buf.data(), send_counts.data(), send_displs.data(), mpi_datatype, recv_buf.data(), recv_counts.data(),
        recv_displs.data(), mpi_datatype, comm);
    send_buf.clear();
    send_buf.resize(0);
    return recv_buf;
}

template <typename Comparator = std::less<Edgelist::value_type>>
inline void
SortEdgesAndWeights(Edgelist& edges, EdgeWeights& edge_weights, Comparator cmp = Comparator{}) {
    if (!std::is_sorted(edges.begin(), edges.end(), cmp)) {
        const SInt num_local_edges = edges.size();
        // If we have edge weights, sort them the same way as the edges
        // This not very efficient; ideally, we should probably implement some kind of zip iterator to sort edges
        // and edge weights without extra allocation / expensive permutation step (@todo)
        if (!edge_weights.empty()) {
            std::vector<EdgeWeights::value_type> indices(num_local_edges);
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&](const auto& lhs, const auto& rhs) {
                return cmp(edges[lhs], edges[rhs]);
            });
            for (std::size_t e = 0; e < num_local_edges; ++e) {
                indices[e] = edge_weights[indices[e]];
            }
            std::swap(edge_weights, indices);
        }

        std::sort(edges.begin(), edges.end(), cmp);
    }
}

inline void RemoveDuplicates(Edgelist& edges, EdgeWeights& edge_weights) {
    const SInt num_local_edges = edges.size();
    if (!edge_weights.empty()) {
        // TODO replace with zip view once C++23 is enabled
        using Edge = typename Edgelist::value_type;
        using Weight = typename EdgeWeights::value_type;

        std::vector<std::pair<Edge, Weight>> edge_weights_zip;
        edge_weights_zip.reserve(num_local_edges);
        for (size_t i = 0; i < num_local_edges; ++i) {
            edge_weights_zip.emplace_back(edges[i], edge_weights[i]);
        }
        edges.clear();
        edge_weights.clear();
        auto it = std::unique(edge_weights_zip.begin(), edge_weights_zip.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.first == rhs.first;
        });
        edge_weights_zip.erase(it, edge_weights_zip.end());
        for (const auto& [edge, weight] : edge_weights_zip) {
            edges.push_back(edge);
            edge_weights.push_back(weight);
        }
    } else {
        auto it = std::unique(edges.begin(), edges.end());
        edges.erase(it, edges.end());
    }
}

} // namespace kagen
