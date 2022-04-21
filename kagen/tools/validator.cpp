#include "kagen/tools/validator.h"

#include <algorithm>
#include <iostream>
#include <numeric>

#include <mpi.h>

namespace kagen {
namespace {
PEID FindPEInRange(const SInt node, const std::vector<std::pair<SInt, SInt>>& ranges) {
    for (std::size_t i = 0; i < ranges.size(); ++i) {
        const auto& [local_from, local_to] = ranges[i];

        if (local_from <= node && node < local_to) {
            return i;
        }
    }

    return -1;
}

std::vector<VertexRange> AllgatherVertexRange(VertexRange vertex_range) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<VertexRange> ranges(static_cast<std::size_t>(size));
    ranges[static_cast<std::size_t>(rank)] = vertex_range;
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, ranges.data(), sizeof(VertexRange), MPI_BYTE, MPI_COMM_WORLD);

    return ranges;
}
} // namespace

bool ValidateVertexRanges(const EdgeList& edge_list, VertexRange vertex_range, const bool expect_consecutive) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const auto ranges = AllgatherVertexRange(vertex_range);

    if (static_cast<std::size_t>(size) != ranges.size()) {
        std::cerr << "Number of vertex ranges (" << ranges.size() << ") differs from the size of the MPI_COMM_WORLD ("
                  << size << ")\n";
    }

    for (std::size_t i = 0; i < ranges.size(); ++i) {
        const auto& [from, to] = ranges[i];
        if (from > to) {
            std::cerr << "Invalid vertex range on PE " << i << ": " << from << ".." << to << "\n";
            return false;
        }
    }

    if (expect_consecutive) {
        if (ranges.front().first != 0) {
            std::cerr << "Expected consecutive vertex ranges, but nodes on PE 0 do not start at 0, but "
                      << ranges.front().first << "\n";
            return false;
        }

        for (std::size_t i = 1; i < ranges.size(); ++i) {
            if (ranges[i].first != ranges[i - 1].second) {
                std::cerr << "Expected consecutive vertex ranges, but end of PE " << i - 1 << " ("
                          << ranges[i - 1].second << ") differs from start of PE " << i << " (" << ranges[i].first
                          << ")\n";
                return false;
            }
        }
    }

    const auto& [local_from, local_to] = ranges[rank];
    const SInt global_n                = ranges.back().second;

    for (const auto& [tail, head]: edge_list) {
        if (tail < local_from || tail >= local_to) {
            std::cerr << "Tail of edge (" << tail << " --> " << head << ") is out of range [" << local_from << ", "
                      << local_to << ")\n";
            return false;
        }

        if (head >= global_n) {
            std::cerr << "Head of edge (" << tail << " --> " << head << ") is outside the global vertex range\n";
            return false;
        }
    }

    return true;
}

bool ValidateSimpleGraph(EdgeList& edge_list, VertexRange vertex_range) {
    // Validate vertex ranges first
    if (!ValidateVertexRanges(edge_list, vertex_range, true)) {
        return false; // failed, following checks could crash if vertex ranges are broken
    }

    const auto ranges = AllgatherVertexRange(vertex_range);

    // Sort edges to allow binary search to find reverse edges
    if (!std::is_sorted(edge_list.begin(), edge_list.end())) {
        std::sort(edge_list.begin(), edge_list.end());
    }

    // Check that there are no self-loops
    for (const auto& [from, to]: edge_list) {
        if (from == to) {
            std::cerr << "Graph contains self-loops\n";
            return false;
        }
    }

    // Check that there are no duplicate edges
    {
        const std::size_t size_before = edge_list.size();

        auto it = std::unique(edge_list.begin(), edge_list.end());
        edge_list.erase(it, edge_list.end());

        if (size_before != edge_list.size()) {
            std::cerr << "Graph contains duplicated edges\n";
            return false;
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
                std::cerr << "Missing reverse edge " << v << " --> " << u << " (internal)\n";
                return false;
            }
        }
    }

    // Check that there are reverse edges for edges across PEs
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<std::vector<SInt>> message_buffers(size);
    for (const auto& [u, v]: edge_list) {
        if (v < from || v >= to) {
            const SInt pe = static_cast<SInt>(FindPEInRange(v, ranges));
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
            std::cerr << "Missing reverse edge " << v << " --> " << u << " (external)\n";
            return false;
        }
    }

    return true;
}
} // namespace kagen
