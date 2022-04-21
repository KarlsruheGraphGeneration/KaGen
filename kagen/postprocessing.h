/*******************************************************************************
 * include/postprocessing.h
 *
 * Copyright (C) 2022 Niklas Uhl <uhl@ira.uka.de>
 * Copyright (C) 2022 Daniel Seemaier <seemaier@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <algorithm>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

#include <mpi.h>

#include "definitions.h"
#include "generator_config.h"

#define KAGEN_MAX_WARNINGS 10

namespace kagen {
namespace internal {
inline PEID FindPEInRange(const SInt node, const std::vector<std::pair<SInt, SInt>>& ranges) {
    for (std::size_t i = 0; i < ranges.size(); ++i) {
        const auto& [local_from, local_to] = ranges[i];

        if (local_from <= node && node < local_to) {
            return i;
        }
    }

    return -1;
}
} // namespace internal

inline bool
ValidateVertexRanges(const EdgeList& edge_list, const std::vector<VertexRange>& ranges, const bool expect_consecutive) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::size_t num_errors = 0;

    if (static_cast<std::size_t>(size) != ranges.size()) {
        std::cerr << "Number of vertex ranges (" << ranges.size() << ") differs from the size of the MPI_COMM_WORLD ("
                  << size << ")\n";
    }

    for (std::size_t i = 0; i < ranges.size(); ++i) {
        const auto& [from, to] = ranges[i];
        if (from > to) {
            std::cerr << "Invalid vertex range on PE " << i << ": " << from << ".." << to << "\n";
            ++num_errors;

            if (num_errors >= KAGEN_MAX_WARNINGS) {
                return false;
            }
        }
    }

    if (expect_consecutive) {
        if (ranges.front().first != 0) {
            std::cerr << "Expected consecutive vertex ranges, but nodes on PE 0 do not start at 0, but "
                      << ranges.front().first << "\n";
        }
        for (std::size_t i = 1; i < ranges.size(); ++i) {
            if (ranges[i].first != ranges[i - 1].second) {
                std::cerr << "Expected consecutive vertex ranges, but end of PE " << i - 1 << " ("
                          << ranges[i - 1].second << ") differs from start of PE " << i << " (" << ranges[i].first
                          << ")\n";
                ++num_errors;

                if (num_errors >= KAGEN_MAX_WARNINGS) {
                    return false;
                }
            }
        }
    }

    const auto& [local_from, local_to] = ranges[rank];
    const SInt global_n                = ranges.back().second;

    for (const auto& [tail, head]: edge_list) {
        if (tail < local_from || tail >= local_to) {
            std::cerr << "Tail of edge (" << tail << " --> " << head << ") is out of range [" << local_from << ", "
                      << local_to << ")\n";
            ++num_errors;
        }

        if (head >= global_n) {
            std::cerr << "Head of edge (" << tail << " --> " << head << ") is outside the global vertex range\n";
            ++num_errors;
        }

        if (num_errors >= KAGEN_MAX_WARNINGS) {
            return false;
        }
    }

    return num_errors == 0;
}

inline void ValidateUndirectedEdgeList(EdgeList& edge_list, const std::vector<VertexRange>& ranges) {
    // Validate vertex ranges
    if (!ValidateVertexRanges(edge_list, ranges, true)) {
        return; // failed, following checks could crash if vertex ranges are broken
    }

    // Sort edges to allow binary search to find reverse edges
    if (!std::is_sorted(edge_list.begin(), edge_list.end())) {
        std::sort(edge_list.begin(), edge_list.end());
    }

    // Check that there are no self-loops
    {
        for (const auto& [from, to]: edge_list) {
            if (from == to) {
                std::cout << "Warning: there are self-loops\n";
                break;
            }
        }
    }

    // Check that there are no duplicate edges
    {
        const std::size_t size_before = edge_list.size();

        auto it = std::unique(edge_list.begin(), edge_list.end());
        edge_list.erase(it, edge_list.end());

        if (size_before != edge_list.size()) {
            std::cerr << "Warning: there are duplicate edges\n";
        }
    }

    // Precompute offset for each node
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const auto [from, to] = ranges[rank];

    std::vector<SInt> node_offset(to - from + 1);
    for (const auto& [u, v]: edge_list) {
        ++node_offset[u - from + 1];
    }
    std::partial_sum(node_offset.begin(), node_offset.end(), node_offset.begin());

    // Check that there are reverse edges for local edges
    for (const auto& [u, v]: edge_list) {
        if (from <= v && v < to) {
            if (!std::binary_search(
                    edge_list.begin() + node_offset[v - from], edge_list.begin() + node_offset[v - from + 1],
                    std::make_tuple(v, u))) {
                std::cerr << "Warning: missing reverse edge " << v << " --> " << u << " (internal)\n";
                break;
            }
        }
    }

    // Check that there are reverse edges for edges across PEs
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<std::vector<SInt>> message_buffers(size);
    for (const auto& [u, v]: edge_list) {
        if (v < from || v >= to) {
            const SInt pe = static_cast<SInt>(internal::FindPEInRange(v, ranges));
            message_buffers[pe].emplace_back(u);
            message_buffers[pe].emplace_back(v);
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
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    std::exclusive_scan(recv_counts.begin(), recv_counts.end(), recv_displs.begin(), 0);
    const std::size_t total_recv_count = recv_displs.back() + recv_counts.back();

    send_buf.reserve(total_send_count);
    for (std::size_t i = 0; i < send_counts.size(); ++i) {
        for (const auto& elem: message_buffers[i]) {
            send_buf.push_back(elem);
        }
        { [[maybe_unused]] auto clear = std::move(message_buffers[i]); }
    }

    recv_buf.resize(total_recv_count);
    MPI_Alltoallv(
        send_buf.data(), send_counts.data(), send_displs.data(), MPI_UINT64_T, recv_buf.data(), recv_counts.data(),
        recv_displs.data(), MPI_UINT64_T, MPI_COMM_WORLD);

    for (std::size_t i = 0; i < recv_buf.size(); i += 2) {
        const SInt u = recv_buf[i];
        const SInt v = recv_buf[i + 1];

        // Check that v --> u exists
        if (!std::binary_search(
                edge_list.begin() + node_offset[v - from], edge_list.begin() + node_offset[v - from + 1],
                std::make_tuple(v, u))) {
            std::cerr << "Warning: missing reverse edge " << v << " --> " << u << " (external)\n";
            break;
        }
    }
}

inline void FixEdgeList(EdgeList& edge_list, const std::vector<VertexRange>& ranges) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const SInt edge_list_size_before = edge_list.size();

    const auto& [local_from, local_to] = ranges[rank];
    std::unordered_map<PEID, std::vector<SInt>> message_buffers;

    // Each PE gets the edges that we have to that PE
    for (const auto& [tail, head]: edge_list) {
        if ((tail >= local_from && tail < local_to) && (head < local_from || head >= local_to)) {
            const SInt pe = static_cast<SInt>(internal::FindPEInRange(head, ranges));
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
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
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
        recv_displs.data(), MPI_UINT64_T, MPI_COMM_WORLD);
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
    MPI_Reduce(
        &edge_list_size_before, &edge_list_global_size_before, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT,
        MPI_COMM_WORLD);
    SInt edge_list_global_size_after = 0;
    MPI_Reduce(
        &edge_list_size_after, &edge_list_global_size_after, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, MPI_COMM_WORLD);

    if (rank == ROOT) {
        std::cout << "- Number of edges (before): " << edge_list_global_size_before << "\n";
        std::cout << "- Number of edges  (after): " << edge_list_global_size_after << "\n";
        std::cout << "- Changed by: ............. "
                  << static_cast<SSInt>(edge_list_global_size_after) - static_cast<SSInt>(edge_list_global_size_before)
                  << " edges\n";
    }
}

inline void
RedistributeGraph(EdgeList& edge_list, const std::vector<VertexRange>& ranges, const bool make_undirected = true) {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const auto from = ranges[rank].first;
    const auto to   = ranges[rank].second;

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
            const PEID owner = internal::FindPEInRange(u, ranges);
            remote_edges[owner].emplace_back(u, v);
        }

        if (make_undirected) {
            if (from <= v && v < to) {
                local_edges.emplace_back(v, u);
            } else {
                const PEID owner = internal::FindPEInRange(v, ranges);
                remote_edges[owner].emplace_back(v, u);
            }
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
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
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
        recv_counts.data(), recv_displs.data(), MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
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

inline void Postprocess(Postprocessing option, EdgeList& edge_list, const std::vector<VertexRange>& ranges) {
    switch (option) {
        case Postprocessing::VALIDATE_RANGES:
        case Postprocessing::VALIDATE_RANGES_CONSECUTIVE:
            ValidateVertexRanges(edge_list, ranges, option == Postprocessing::VALIDATE_RANGES_CONSECUTIVE);
            return;

        case Postprocessing::VALIDATE_UNDIRECTED:
            ValidateUndirectedEdgeList(edge_list, ranges);
            return;

        case Postprocessing::FIX_UNDIRECTED_EDGE_LIST:
            FixEdgeList(edge_list, ranges);
            return;

        case Postprocessing::REDISTRIBUTE_GRAPH:
            RedistributeGraph(edge_list, ranges);
            return;

        case Postprocessing::SKIP:
            return;
    }

    __builtin_unreachable();
}

inline void Postprocess(Postprocessing option, EdgeList& edges, VertexRange vertex_range) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<VertexRange> ranges(static_cast<std::size_t>(size));
    ranges[static_cast<std::size_t>(rank)] = vertex_range;

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, ranges.data(), sizeof(VertexRange), MPI_BYTE, MPI_COMM_WORLD);

    Postprocess(option, edges, ranges);
}
} // namespace kagen
