#include "kagen/tools/validator.h"

#include "kagen/tools/converter.h"

#include <mpi.h>

#include <algorithm>
#include <iostream>
#include <numeric>

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

bool ValidateVertexRanges(const Edgelist& edge_list, const VertexRange vertex_range, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const auto ranges = AllgatherVertexRange(vertex_range, comm);

    if (static_cast<std::size_t>(size) != ranges.size()) {
        std::cerr << "Number of vertex ranges (" << ranges.size() << ") differs from the size of the communicator ("
                  << size << ")\n";
    }

    for (std::size_t i = 0; i < ranges.size(); ++i) {
        const auto& [from, to] = ranges[i];
        if (from > to) {
            std::cerr << "Invalid vertex range on PE " << i << ": " << from << ".." << to << "\n";
            return false;
        }
    }

    if (ranges.front().first != 0) {
        std::cerr << "Expected consecutive vertex ranges, but nodes on PE 0 do not start at 0, but "
                  << ranges.front().first << "\n";
        return false;
    }

    for (std::size_t i = 1; i < ranges.size(); ++i) {
        if (ranges[i].first != ranges[i - 1].second) {
            std::cerr << "Expected consecutive vertex ranges, but end of PE " << i - 1 << " (" << ranges[i - 1].second
                      << ") differs from start of PE " << i << " (" << ranges[i].first << ")\n";
            return false;
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

bool ValidateGraph(
    Graph& graph, const bool allow_self_loops, const bool allow_directed_graphs, const bool allow_multi_edges,
    MPI_Comm comm) {
    // Validation in CSR representation is currently not implemented, so we convert to it to edge list
    if (graph.representation == GraphRepresentation::CSR) {
        graph.edges       = BuildEdgeListFromCSR(graph.vertex_range, graph.xadj, graph.adjncy);
        const bool result = ValidateGraph(graph, allow_self_loops, allow_directed_graphs, allow_multi_edges, comm);
        { [[maybe_unused]] auto tmp = std::move(graph.edges); } // Free edge list
        return result;
    }

    auto& edge_list      = graph.edges;
    auto& vertex_range   = graph.vertex_range;
    auto& vertex_weights = graph.vertex_weights;
    auto& edge_weights   = graph.edge_weights;

    // Validate vertex ranges first
    if (!ValidateVertexRanges(edge_list, vertex_range, comm)) {
        return false; // failed, following checks could crash if vertex ranges are broken
    }

    const auto ranges = AllgatherVertexRange(vertex_range, comm);

    if (!vertex_weights.empty() && vertex_weights.size() != vertex_range.second - vertex_range.first) {
        std::cerr << "There are " << vertex_weights.size() << " vertex weights for "
                  << vertex_range.second - vertex_range.first << " vertices\n";
        return false;
    }
    if (!edge_weights.empty() && edge_list.size() != edge_weights.size()) {
        std::cerr << "There are " << edge_weights.size() << " edge weights for " << edge_list.size() << " edges\n";
        return false;
    }

    // Sort edges to allow binary search to find reverse edges
    using Edge = std::tuple<SInt, SInt, SSInt>; // from, to, weight
    std::vector<Edge> sorted_edges(edge_list.size());
    std::transform(edge_list.begin(), edge_list.end(), sorted_edges.begin(), [&](const auto& edge) {
        return std::make_tuple(std::get<0>(edge), std::get<1>(edge), 1);
    });
    if (!edge_weights.empty()) {
        for (std::size_t e = 0; e < edge_weights.size(); ++e) {
            std::get<2>(sorted_edges[e]) = edge_weights[e];
        }
    }

    std::sort(sorted_edges.begin(), sorted_edges.end());

    // Check that there are no self-loops
    if (!allow_self_loops) {
        for (const auto& [from, to, weight]: sorted_edges) {
            if (from == to) {
                std::cerr << "Graph contains a self-loops: " << from << " --> " << to << "; aborting\n";
                return false;
            }
        }
    }

    // Check that there are no duplicate edges
    if (!allow_multi_edges) {
        for (std::size_t i = 1; i < sorted_edges.size(); ++i) {
            if (std::get<0>(sorted_edges[i - 1]) == std::get<0>(sorted_edges[i])
                && std::get<1>(sorted_edges[i - 1]) == std::get<1>(sorted_edges[i])) {
                const auto& [from, to, weight] = sorted_edges[i];
                std::cerr << "Graph contains a duplicated edge: " << from << " --> " << to << "; aborting\n";
                return false;
            }
        }
    }

    // Precompute offset for each node
    if (!allow_directed_graphs) {
        int rank;
        MPI_Comm_rank(comm, &rank);
        const auto [from, to] = ranges[rank];

        std::vector<SInt> node_offset(to - from + 1);
        for (const auto& [u, v, weight]: sorted_edges) {
            ++node_offset[u - from + 1];
        }
        std::partial_sum(node_offset.begin(), node_offset.end(), node_offset.begin());

        // Check that there are reverse edges for local edges
        for (const auto& [u, v, weight]: sorted_edges) {
            if (from <= v && v < to) {
                if (!std::binary_search(
                        sorted_edges.begin() + node_offset[v - from], sorted_edges.begin() + node_offset[v - from + 1],
                        std::make_tuple(v, u, weight))) {
                    std::cerr << "Missing reverse edge " << v << " --> " << u << " with weight " << weight
                              << " (internal); the reverse edge might exist with a different edge weight\n";
                    return false;
                }
            }
        }

        // Check that there are reverse edges for edges across PEs
        int size;
        MPI_Comm_size(comm, &size);

        std::vector<std::vector<SInt>> message_buffers(size);
        for (const auto& [u, v, weight]: sorted_edges) {
            if (v < from || v >= to) {
                const SInt pe = static_cast<SInt>(FindPEInRange(v, ranges));
                message_buffers[pe].emplace_back(u);
                message_buffers[pe].emplace_back(v);
                message_buffers[pe].emplace_back(static_cast<SInt>(weight));
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
        for (std::size_t i = 0; i < send_counts.size(); ++i) {
            for (const auto& elem: message_buffers[i]) {
                send_buf.push_back(elem);
            }
            { [[maybe_unused]] auto clear = std::move(message_buffers[i]); }
        }

        recv_buf.resize(total_recv_count);
        MPI_Alltoallv(
            send_buf.data(), send_counts.data(), send_displs.data(), MPI_UINT64_T, recv_buf.data(), recv_counts.data(),
            recv_displs.data(), MPI_UINT64_T, comm);

        for (std::size_t i = 0; i < recv_buf.size(); i += 3) {
            const SInt  u      = recv_buf[i];
            const SInt  v      = recv_buf[i + 1];
            const SSInt weight = static_cast<SSInt>(recv_buf[i + 2]);

            // Check that v --> u exists
            if (!std::binary_search(
                    sorted_edges.begin() + node_offset[v - from], sorted_edges.begin() + node_offset[v - from + 1],
                    std::make_tuple(v, u, weight))) {
                std::cerr << "Missing reverse edge " << v << " --> " << u << " with weight " << weight
                          << " (external); the reverse edge might exist with a different edge weight\n";
                return false;
            }
        }
    }

    return true;
}
} // namespace kagen
