#pragma once

#include "kagen/comm/comm.h"
#include "kagen/comm/comm_types.h"
#include "kagen/kagen.h"

#include <cstring>
#include <limits>
#include <numeric>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace kagen {
template <typename T>
inline T DivideOrDefault(T dividend, T divisor, T default_value = static_cast<T>(0)) {
    static_assert(std::is_arithmetic_v<T>, "DivideOrDefault requires an arithmetic type");
    return divisor == static_cast<T>(0) ? default_value : dividend / divisor;
}

inline std::pair<SInt, SInt> ComputeRange(const SInt n, const PEID size, const PEID rank) {
    const SInt chunk = n / size;
    const SInt rem   = n % size;
    const SInt from  = rank * chunk + std::min<SInt>(rank, rem);
    const SInt to    = std::min<SInt>(from + ((static_cast<SInt>(rank) < rem) ? chunk + 1 : chunk), n);
    return {from, to};
}

inline SInt FindNumberOfVerticesInEdgelist(const Edgelist& edges, Comm& comm) {
    SInt n = 0;
    for (const auto& [u, v]: edges) {
        n = std::max(n, std::max(u, v));
    }
    comm.Allreduce(COMM_IN_PLACE, &n, 1, CommDatatype::UNSIGNED_LONG_LONG, CommOp::MAX);
    return n + 1;
}

inline PEID GetCommRank(Comm& comm) {
    return comm.Rank();
}

inline PEID GetCommSize(Comm& comm) {
    return comm.Size();
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

inline PEID FindPEInRange(const SInt node, const std::vector<VertexRange>& ranges) {
    for (std::size_t i = 0; i < ranges.size(); ++i) {
        const auto& [local_from, local_to] = ranges[i];

        if (local_from <= node && node < local_to) {
            return i;
        }
    }

    return -1;
}

inline PEID FindPEInRangeWithBinarySearch(const SInt node, const std::vector<VertexRange>& ranges) {
    // Find first range whose .second > node
    auto it = std::upper_bound(ranges.begin(), ranges.end(), node, [](SInt value, const VertexRange& range) {
        return value < range.second;
    });
    if (it != ranges.end() && it->first <= node && node < it->second) {
        return std::distance(ranges.begin(), it);
    }
    return -1;
}

inline std::vector<VertexRange> AllgatherVertexRange(const VertexRange vertex_range, Comm& comm) {
    PEID rank = comm.Rank();
    PEID size = comm.Size();

    std::vector<VertexRange> ranges(static_cast<std::size_t>(size));
    ranges[static_cast<std::size_t>(rank)] = vertex_range;
    // Use BYTE-based allgather: each element is sizeof(VertexRange) bytes
    comm.Allgather(
        COMM_IN_PLACE, 0, CommDatatype::BYTE, ranges.data(), static_cast<int>(sizeof(VertexRange)),
        CommDatatype::BYTE);

    return ranges;
}

// Exchange message buffers all-to-all using byte-based communication.
// Each element of type T is sent as sizeof(T) bytes.
template <typename T>
std::vector<T> ExchangeMessageBuffers(std::unordered_map<PEID, std::vector<T>> message_buffers, Comm& comm) {
    PEID rank = comm.Rank();
    PEID size = comm.Size();

    std::vector<T>   send_buf;
    std::vector<T>   recv_buf;
    std::vector<int> send_counts(size);
    std::vector<int> recv_counts(size);
    std::vector<int> send_displs(size);
    std::vector<int> recv_displs(size);

    for (size_t i = 0; i < static_cast<std::size_t>(size); ++i) {
        send_counts[i] = message_buffers[i].size();
    }

    std::exclusive_scan(send_counts.begin(), send_counts.end(), send_displs.begin(), 0);
    const std::size_t total_send_count = send_displs.back() + send_counts.back();

    comm.Alltoall(send_counts.data(), 1, CommDatatype::INT, recv_counts.data(), 1, CommDatatype::INT);

    std::exclusive_scan(recv_counts.begin(), recv_counts.end(), recv_displs.begin(), 0);
    const std::size_t total_recv_count = recv_displs.back() + recv_counts.back();

    send_buf.reserve(total_send_count);
    for (size_t i = 0; i < static_cast<std::size_t>(size); ++i) {
        for (const auto& elem: message_buffers[i]) {
            send_buf.push_back(elem);
        }
        message_buffers[i].clear();
        message_buffers[i].resize(0);
    }
    recv_buf.resize(total_recv_count);

    // Convert element counts to byte counts for BYTE-based alltoallv
    const int element_size = static_cast<int>(sizeof(T));
    std::vector<int> send_counts_bytes(size), recv_counts_bytes(size);
    std::vector<int> send_displs_bytes(size), recv_displs_bytes(size);
    for (int i = 0; i < size; ++i) {
        send_counts_bytes[i] = send_counts[i] * element_size;
        recv_counts_bytes[i] = recv_counts[i] * element_size;
        send_displs_bytes[i] = send_displs[i] * element_size;
        recv_displs_bytes[i] = recv_displs[i] * element_size;
    }

    comm.Alltoallv(
        send_buf.data(), send_counts_bytes.data(), send_displs_bytes.data(), CommDatatype::BYTE, recv_buf.data(),
        recv_counts_bytes.data(), recv_displs_bytes.data(), CommDatatype::BYTE);

    send_buf.clear();
    send_buf.resize(0);
    return recv_buf;
}

template <typename Comparator = std::less<Edgelist::value_type>>
inline void SortEdgesAndWeights(Edgelist& edges, EdgeWeights& edge_weights, Comparator cmp = Comparator{}) {
    if (!std::is_sorted(edges.begin(), edges.end(), cmp)) {
        const SInt num_local_edges = edges.size();
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
        using Edge   = typename Edgelist::value_type;
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
        for (const auto& [edge, weight]: edge_weights_zip) {
            edges.push_back(edge);
            edge_weights.push_back(weight);
        }
    } else {
        auto it = std::unique(edges.begin(), edges.end());
        edges.erase(it, edges.end());
    }
}

inline void SortAndRemoveDuplicates(Edgelist& edges) {
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
}

} // namespace kagen
