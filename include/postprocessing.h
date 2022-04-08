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
#include "io/generator_io.h"

#define KAGEN_MAX_WARNINGS 10

namespace kagen {
namespace internal {
inline PEID FindPEInRange(const SInt node, const std::vector<std::pair<SInt, SInt>>& ranges) {
    for (std::size_t i = 0; i < ranges.size(); ++i) {
        const auto& [local_from, local_to] = ranges[i];

        if (local_from <= node && node <= local_to) {
            return i;
        }
    }

    return -1;
}
} // namespace internal

inline void
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
                return;
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
                    return;
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
            return;
        }
    }
}

inline void ValidateUndirectedEdgeList(const EdgeList& edge_list, const std::vector<VertexRange>& ranges) {
    ((void)edge_list);
    ((void)ranges);
    // @todo
}

inline void FixEdgeList(EdgeList& edge_list, const std::vector<VertexRange>& ranges) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const auto& [local_from, local_to] = ranges[rank];
    std::unordered_map<PEID, std::vector<SInt>> message_buffers;

    for (const auto& [tail, head]: edge_list) {
        if (tail >= local_from && tail < local_to) {
            if (head < local_from || head >= local_to) {
                message_buffers[internal::FindPEInRange(head, ranges)].emplace_back(tail);
                message_buffers[internal::FindPEInRange(head, ranges)].emplace_back(head);
            }
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
    const std::size_t total_send_count = send_displs[size - 1] + send_counts[size - 1];
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    std::exclusive_scan(recv_counts.begin(), recv_counts.end(), recv_displs.begin(), 0);
    const std::size_t total_recv_count = recv_displs[size - 1] + recv_counts[size - 1];
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

    // kagen sometimes produces duplicate edges
    auto it = std::unique(edge_list.begin(), edge_list.end());
    edge_list.erase(it, edge_list.end());
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

        case Postprocessing::SKIP:
            return;
    }

    __builtin_unreachable();
}

template <typename Generator>
void Postprocess(Postprocessing option, Generator& generator) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<VertexRange> ranges(static_cast<std::size_t>(size));
    ranges[static_cast<std::size_t>(rank)] = generator.GetVertexRange();

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, ranges.data(), sizeof(VertexRange), MPI_BYTE, MPI_COMM_WORLD);

    Postprocess(option, generator.IO().GetEdges(), ranges);
}
} // namespace kagen
