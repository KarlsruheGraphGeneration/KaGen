#include "kagen/tools/statistics.h"

#include "kagen/definitions.h"

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <unordered_set>

namespace kagen {
// First invalid node on the last PE is the number of nodes in the graph
SInt FindNumberOfGlobalNodes(const VertexRange vertex_range, MPI_Comm comm) {
    PEID size;
    MPI_Comm_size(comm, &size);

    SInt first_invalid_node = vertex_range.second;
    MPI_Bcast(&first_invalid_node, 1, MPI_UNSIGNED_LONG_LONG, size - 1, comm);

    return first_invalid_node;
}

// Length of all edge lists is the number of edges in the graph
SInt FindNumberOfGlobalEdges(const Edgelist& edges, MPI_Comm comm) {
    SInt local_num_edges = edges.size();
    SInt global_num_edges;
    MPI_Allreduce(&local_num_edges, &global_num_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
    return global_num_edges;
}

namespace {
std::vector<SInt> GatherValue(const SInt value, MPI_Comm comm) {
    PEID size;
    MPI_Comm_size(comm, &size);
    std::vector<SInt> values(size);
    MPI_Allgather(&value, 1, MPI_UNSIGNED_LONG_LONG, values.data(), 1, MPI_UNSIGNED_LONG_LONG, comm);
    return values;
}
} // namespace

std::vector<SInt> GatherNumberOfEdges(const Edgelist& edges, MPI_Comm comm) {
    return GatherValue(edges.size(), comm);
}

SInt ReduceSum(const SInt value, MPI_Comm comm) {
    SInt sum = 0;
    MPI_Reduce(&value, &sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, comm);
    return sum;
}

SInt ReduceMin(const SInt value, MPI_Comm comm) {
    SInt min = 0;
    MPI_Reduce(&value, &min, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, ROOT, comm);
    return min;
}

LPFloat ReduceMean(const SInt value, MPI_Comm comm) {
    SInt sum = 0;
    MPI_Reduce(&value, &sum, 1, KAGEN_MPI_SINT, MPI_SUM, ROOT, comm);

    PEID size;
    MPI_Comm_size(comm, &size);

    return 1.0 * sum / size;
}

SInt ReduceMax(const SInt value, MPI_Comm comm) {
    SInt max = 0;
    MPI_Reduce(&value, &max, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, ROOT, comm);
    return max;
}

LPFloat ReduceSD(const SInt value, MPI_Comm comm) {
    const auto values = GatherValue(value, comm);
    const auto mean   = ReduceMean(value, comm);

    LPFloat sd_sum = 0.0;
    for (const auto& e: values) {
        sd_sum += (1.0 * e - mean) * (1.0 * e - mean);
    }

    if (sd_sum != 0) { // root
        return std::sqrt(1.0 / (1.0 * values.size()) * sd_sum);
    }
    return 0.0; // non-root
}

DegreeStatistics ReduceDegreeStatistics(const Edgelist& edges, const SInt global_num_nodes, MPI_Comm comm) {
    assert(std::is_sorted(edges.begin(), edges.end()));

    SInt min = std::numeric_limits<SInt>::max();
    SInt sum = 0;
    SInt max = std::numeric_limits<SInt>::lowest();

    SInt cur_from   = edges.empty() ? 0 : std::get<0>(edges.front());
    SInt cur_degree = 0;

    auto update = [&](const SInt deg) {
        min = std::min(min, deg);
        max = std::max(max, deg);
        sum += deg;
    };

    for (const auto& [from, to]: edges) {
        if (from == cur_from) {
            ++cur_degree;
        } else {
            update(cur_degree);
            cur_degree = 1;
            if (cur_from + 1 != from) {
                min = 0; // Skipped dgeree 0 vertex
            }
            cur_from = from;
        }
    }
    update(cur_degree);

    SInt global_min = 0, global_sum = 0, global_max = 0;
    MPI_Reduce(&min, &global_min, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, ROOT, comm);
    MPI_Reduce(&sum, &global_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, comm);
    MPI_Reduce(&max, &global_max, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, ROOT, comm);

    PEID size = 0;
    MPI_Comm_size(comm, &size);

    return {global_min, 1.0 * global_sum / global_num_nodes, global_max};
}

std::vector<SInt> ComputeDegreeBins(const Edgelist& edges, const VertexRange vertex_range, MPI_Comm comm) {
    assert(std::is_sorted(edges.begin(), edges.end()));

    std::vector<SInt> bins(std::numeric_limits<SInt>::digits);
    SInt              cur_from   = edges.empty() ? 0 : std::get<0>(edges.front());
    SInt              cur_degree = 0;

    auto yield = [&](const SInt deg) {
        const SInt bin = (deg == 0) ? 0 : (std::log2(deg) + 1);
        ++bins[bin];
    };

    for (const auto& [from, to]: edges) {
        if (from == cur_from) {
            ++cur_degree;
        } else {
            yield(cur_degree);
            cur_degree = 1;
            while (++cur_from < from) {
                ++bins[0];
            }
        }
    }
    if (!edges.empty()) {
        yield(cur_degree);
        while (++cur_from < vertex_range.second) {
            ++bins[0];
        }
    }

    std::vector<SInt> global_bins(bins.size());
    MPI_Reduce(bins.data(), global_bins.data(), bins.size(), KAGEN_MPI_SINT, MPI_SUM, ROOT, comm);

    return global_bins;
}

double ComputeEdgeLocalicty(const Edgelist& edges, const VertexRange vertex_range, MPI_Comm comm) {
    const SInt num_local_cut_edges = std::count_if(edges.begin(), edges.end(), [&vertex_range](const auto& edge) {
        return std::get<0>(edge) < vertex_range.first || std::get<1>(edge) >= vertex_range.second;
    });
    const SInt num_local_edges     = edges.size();

    SInt num_global_cut_edges = 0;
    SInt num_global_edges     = 0;

    MPI_Reduce(&num_local_cut_edges, &num_global_cut_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, comm);
    MPI_Reduce(&num_local_edges, &num_global_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, comm);

    return 1.0 - 1.0 * num_global_cut_edges / num_global_edges;
}

SInt ComputeNumberOfGhostNodes(const Edgelist& edges, const VertexRange vertex_range, MPI_Comm comm) {
    std::unordered_set<SInt> ghost_nodes;

    for (const auto& [from, to]: edges) {
        if (to < vertex_range.first || to >= vertex_range.second) {
            ghost_nodes.insert(to);
        }
    }

    const SInt num_local_ghost_nodes  = ghost_nodes.size();
    SInt       num_global_ghost_nodes = 0;
    MPI_Reduce(&num_local_ghost_nodes, &num_global_ghost_nodes, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, comm);
    return num_global_ghost_nodes;
}

namespace {
struct DistributedElements {
    SInt   global;
    SInt   min;
    double mean;
    SInt   max;
    double sd;
};

void PrintBasicStatistics(const DistributedElements& vertices, const DistributedElements& edges) {
    const int global_space = std::max<int>(std::log10(vertices.global), std::log10(edges.global)) + 1;
    const int local_space  = std::max<int>(std::log10(vertices.max), std::log10(edges.max)) + 1;

    const double vertex_imbalance = 1.0 * vertices.max / vertices.mean;
    const double edge_imbalance   = 1.0 * edges.max / edges.mean;

    std::cout << "Number of vertices: " << std::setw(global_space) << vertices.global << " ["
              << "Min=" << std::setw(local_space) << vertices.min << " | "
              << "Mean=" << std::setw(local_space + 2) << std::fixed << std::setprecision(1) << vertices.mean << " | "
              << "Max=" << std::setw(local_space) << vertices.max << " | "
              << "SD=" << std::setw(local_space + 3) << std::fixed << std::setprecision(2) << vertices.sd << "]\n";
    std::cout << "  Vertex imbalance: " << std::fixed << std::setprecision(3) << vertex_imbalance << std::endl;
    std::cout << "Number of edges:    " << std::setw(global_space) << edges.global << " ["
              << "Min=" << std::setw(local_space) << edges.min << " | "
              << "Mean=" << std::setw(local_space + 2) << std::fixed << std::setprecision(1) << edges.mean << " | "
              << "Max=" << std::setw(local_space) << edges.max << " | "
              << "SD=" << std::setw(local_space + 3) << std::fixed << std::setprecision(2) << edges.sd << "]\n";
    std::cout << "  Edge imbalance:   " << std::fixed << std::setprecision(3) << edge_imbalance << std::endl;
}

void PrintBasicStatistics(const SInt local_num_vertices, const SInt local_num_edges, const bool root, MPI_Comm comm) {
    const auto global_num_vertices = ReduceSum(local_num_vertices, comm);
    const auto local_min_vertices  = ReduceMin(local_num_vertices, comm);
    const auto local_mean_vertices = ReduceMean(local_num_vertices, comm);
    const auto local_max_vertices  = ReduceMax(local_num_vertices, comm);
    const auto local_sd_vertices   = ReduceSD(local_num_vertices, comm);

    const auto global_num_edges = ReduceSum(local_num_edges, comm);
    const auto local_min_edges  = ReduceMin(local_num_edges, comm);
    const auto local_mean_edges = ReduceMean(local_num_edges, comm);
    const auto local_max_edges  = ReduceMax(local_num_edges, comm);
    const auto local_sd_edges   = ReduceSD(local_num_edges, comm);

    // Print statistics on root
    if (root) {
        PrintBasicStatistics(
            {
                global_num_vertices,
                local_min_vertices,
                local_mean_vertices,
                local_max_vertices,
                local_sd_vertices,
            },
            {
                global_num_edges,
                local_min_edges,
                local_mean_edges,
                local_max_edges,
                local_sd_edges,
            });
    }
}
} // namespace

void PrintBasicStatistics(
    const XadjArray& xadj, const AdjncyArray& adjncy, VertexRange, const bool root, MPI_Comm comm) {
    PrintBasicStatistics(xadj.size() - 1, adjncy.size(), root, comm);
}

void PrintBasicStatistics(const Edgelist& edges, const VertexRange vertex_range, const bool root, MPI_Comm comm) {
    PrintBasicStatistics(vertex_range.second - vertex_range.first, edges.size(), root, comm);
}

void PrintAdvancedStatistics(Edgelist& edges, const VertexRange vertex_range, const bool root, MPI_Comm comm) {
    // Sort edges for degree computation
    if (!std::is_sorted(edges.begin(), edges.end())) {
        std::sort(edges.begin(), edges.end());
    }

    // Compute degree statistics
    const auto local_num_nodes  = vertex_range.second - vertex_range.first;
    const auto global_num_nodes = ReduceSum(local_num_nodes, comm);
    const auto local_num_edges  = edges.size();
    const auto global_num_edges = ReduceSum(local_num_edges, comm);

    const double density = 1.0 * global_num_edges / global_num_nodes / (global_num_nodes - 1);
    const auto [min_degree, mean_degree, max_degree] = ReduceDegreeStatistics(edges, global_num_nodes, comm);
    const auto degree_bins                           = ComputeDegreeBins(edges, vertex_range, comm);

    // Compute locality statistics
    const double edge_locality          = ComputeEdgeLocalicty(edges, vertex_range, comm);
    const SInt   global_num_ghost_nodes = ComputeNumberOfGhostNodes(edges, vertex_range, comm);
    const double ghost_node_fraction    = 1.0 * global_num_ghost_nodes / (global_num_nodes + global_num_ghost_nodes);

    // Print on root
    if (root) {
        std::cout << "Density: " << std::fixed << std::setprecision(4) << density << "\n";
        std::cout << "Degrees: [Min=" << min_degree << " | Mean=" << std::fixed << std::setprecision(1) << mean_degree
                  << " | Max=" << max_degree << "]\n";

        // Find last non-empty degree bin
        SInt last_nonempty_degree_bin = 0;
        for (SInt i = 0; i < degree_bins.size(); ++i) {
            if (degree_bins[i] > 0) {
                last_nonempty_degree_bin = i;
            }
        }

        // Print degree bins
        const SInt digits10 = std::log10(1 << last_nonempty_degree_bin) + 1;

        std::cout << "Degree bins:\n";
        for (SInt i = 0; i <= last_nonempty_degree_bin; ++i) {
            const SInt from = (i == 0) ? 0 : 1 << (i - 1);
            const SInt to   = 2 * from;
            std::cout << "  Degree in [" << std::setw(digits10) << from << ", " << std::setw(digits10) << to
                      << "): " << degree_bins[i] << "\n";
        }

        // Print locality statistics
        std::cout << "Edge locality: " << std::fixed << std::setprecision(4) << edge_locality << std::endl;
        std::cout << "Fraction of ghost nodes: " << std::fixed << std::setprecision(4) << ghost_node_fraction
                  << std::endl;
        std::cout << "  There are " << global_num_nodes << " real vertices and " << global_num_ghost_nodes
                  << " ghost vertices" << std::endl;
    }
}
} // namespace kagen
