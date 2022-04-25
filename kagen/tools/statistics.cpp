#include "kagen/tools/statistics.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <unordered_set>

#include <mpi.h>

namespace kagen {
// First invalid node on the last PE is the number of nodes in the graph
SInt FindNumberOfGlobalNodes(const VertexRange vertex_range) {
    PEID size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    SInt first_invalid_node = vertex_range.second;
    MPI_Bcast(&first_invalid_node, 1, MPI_UNSIGNED_LONG_LONG, size - 1, MPI_COMM_WORLD);

    return first_invalid_node;
}

// Length of all edge lists is the number of edges in the graph
SInt FindNumberOfGlobalEdges(const EdgeList& edges) {
    SInt local_num_edges = edges.size();
    SInt global_num_edges;
    MPI_Allreduce(&local_num_edges, &global_num_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    return global_num_edges;
}

namespace {
std::vector<SInt> GatherValue(const SInt value) {
    PEID size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<SInt> values(size);
    MPI_Allgather(&value, 1, MPI_UNSIGNED_LONG_LONG, values.data(), 1, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
    return values;
}
} // namespace

std::vector<SInt> GatherNumberOfEdges(const EdgeList& edges) {
    return GatherValue(edges.size());
}

SInt ReduceSum(const SInt value) {
    SInt sum;
    MPI_Reduce(&value, &sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, MPI_COMM_WORLD);
    return sum;
}

SInt ReduceMin(const SInt value) {
    SInt min;
    MPI_Reduce(&value, &min, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, ROOT, MPI_COMM_WORLD);
    return min;
}

LPFloat ReduceMean(const SInt value) {
    SInt sum;
    MPI_Reduce(&value, &sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, MPI_COMM_WORLD);

    PEID size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    return 1.0 * sum / size;
}

SInt ReduceMax(const SInt value) {
    SInt max;
    MPI_Reduce(&value, &max, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, ROOT, MPI_COMM_WORLD);
    return max;
}

LPFloat ReduceSD(const SInt value) {
    const auto values = GatherValue(value);
    const auto mean   = ReduceMean(value);

    LPFloat sd_sum = 0.0;
    for (const auto& e: values) {
        sd_sum += (1.0 * e - mean) * (1.0 * e - mean);
    }

    if (sd_sum != 0) { // root
        return std::sqrt(1.0 / (1.0 * values.size()) * sd_sum);
    }
    return 0.0; // non-root
}

DegreeStatistics ReduceDegreeStatistics(const EdgeList& edges, const SInt global_num_nodes) {
    assert(std::is_sorted(edges.begin(), edges.end()));

    SInt min = std::numeric_limits<SInt>::max();
    SInt sum = 0;
    SInt max = std::numeric_limits<SInt>::lowest();

    SInt cur_from   = std::get<0>(edges.front());
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

    SInt global_min, global_sum, global_max;
    MPI_Reduce(&min, &global_min, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, ROOT, MPI_COMM_WORLD);
    MPI_Reduce(&sum, &global_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, MPI_COMM_WORLD);
    MPI_Reduce(&max, &global_max, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, ROOT, MPI_COMM_WORLD);

    PEID size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    return {global_min, 1.0 * global_sum / global_num_nodes, global_max};
}

std::vector<SInt> ComputeDegreeBins(const EdgeList& edges, const VertexRange vertex_range) {
    assert(std::is_sorted(edges.begin(), edges.end()));

    std::vector<SInt> bins(std::numeric_limits<SInt>::digits);
    SInt              cur_from   = std::get<0>(edges.front());
    SInt              cur_degree = 0;

    auto yield = [&](const SInt deg) {
        const SInt bin = std::log2(deg) + 1;
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
    yield(cur_degree);
    while (++cur_from < vertex_range.second) {
        ++bins[0];
    }

    std::vector<SInt> global_bins(bins.size());
    MPI_Reduce(bins.data(), global_bins.data(), bins.size(), MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, MPI_COMM_WORLD);

    return global_bins;
}

double ComputeEdgeLocalicty(const EdgeList& edges, const VertexRange vertex_range) {
    const SInt num_local_cut_edges = std::count_if(edges.begin(), edges.end(), [&vertex_range](const auto& edge) {
        return std::get<0>(edge) < vertex_range.first || std::get<1>(edge) >= vertex_range.second;
    });
    const SInt num_local_edges     = edges.size();

    SInt num_global_cut_edges;
    SInt num_global_edges;

    MPI_Reduce(&num_local_cut_edges, &num_global_cut_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, MPI_COMM_WORLD);
    MPI_Reduce(&num_local_edges, &num_global_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, MPI_COMM_WORLD);

    return 1.0 * num_global_cut_edges / num_global_edges;
}

SInt ComputeNumberOfGhostNodes(const EdgeList& edges, const VertexRange vertex_range) {
    std::unordered_set<SInt> ghost_nodes;

    for (const auto& [from, to]: edges) {
        if (to < vertex_range.first || to >= vertex_range.second) {
            ghost_nodes.insert(to);
        }
    }

    const SInt num_local_ghost_nodes = ghost_nodes.size();
    SInt       num_global_ghost_nodes;
    MPI_Reduce(
        &num_local_ghost_nodes, &num_global_ghost_nodes, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, MPI_COMM_WORLD);
    return num_global_ghost_nodes;
}

void PrintBasicStatistics(const EdgeList& edges, const VertexRange vertex_range, const bool root) {
    // Compute statistics
    const auto local_num_nodes  = vertex_range.second - vertex_range.first;
    const auto global_num_nodes = ReduceSum(local_num_nodes);
    const auto local_min_nodes  = ReduceMin(local_num_nodes);
    const auto local_mean_nodes = ReduceMean(local_num_nodes);
    const auto local_max_nodes  = ReduceMax(local_num_nodes);
    const auto local_sd_nodes   = ReduceSD(local_num_nodes);

    const auto local_num_edges  = edges.size();
    const auto global_num_edges = ReduceSum(local_num_edges);
    const auto local_min_edges  = ReduceMin(local_num_edges);
    const auto local_mean_edges = ReduceMean(local_num_edges);
    const auto local_max_edges  = ReduceMax(local_num_edges);
    const auto local_sd_edges   = ReduceSD(local_num_edges);

    const double edge_imbalance = 1.0 * local_max_edges / local_mean_edges;

    // Print statistics on root
    if (root) {
        const int global_space = std::max<int>(std::log10(global_num_nodes), std::log10(global_num_edges)) + 1;
        const int local_space  = std::max<int>(std::log10(local_max_nodes), std::log10(local_max_edges)) + 1;

        std::cout << "Number of vertices: " << std::setw(global_space) << global_num_nodes << " ["
                  << "Min=" << std::setw(local_space) << local_min_nodes << " | "
                  << "Mean=" << std::setw(local_space + 2) << std::fixed << std::setprecision(1) << local_mean_nodes
                  << " | "
                  << "Max=" << std::setw(local_space) << local_max_nodes << " | "
                  << "SD=" << std::setw(local_space + 3) << std::fixed << std::setprecision(2) << local_sd_nodes
                  << "]\n";
        std::cout << "Number of edges:    " << std::setw(global_space) << global_num_edges << " ["
                  << "Min=" << std::setw(local_space) << local_min_edges << " | "
                  << "Mean=" << std::setw(local_space + 2) << std::fixed << std::setprecision(1) << local_mean_edges
                  << " | "
                  << "Max=" << std::setw(local_space) << local_max_edges << " | "
                  << "SD=" << std::setw(local_space + 3) << std::fixed << std::setprecision(2) << local_sd_edges
                  << "]\n";
        std::cout << "  Edge imbalance: " << std::fixed << std::setprecision(3) << edge_imbalance << std::endl;
    }
}

void PrintAdvancedStatistics(EdgeList& edges, const VertexRange vertex_range, const bool root) {
    // Sort edges for degree computation
    if (!std::is_sorted(edges.begin(), edges.end())) {
        std::sort(edges.begin(), edges.end());
    }

    // Compute degree statistics
    const auto local_num_nodes  = vertex_range.second - vertex_range.first;
    const auto global_num_nodes = ReduceSum(local_num_nodes);
    const auto local_num_edges  = edges.size();
    const auto global_num_edges = ReduceSum(local_num_edges);

    const double density = 1.0 * global_num_edges / global_num_nodes / (global_num_nodes - 1);
    const auto [min_degree, mean_degree, max_degree] = ReduceDegreeStatistics(edges, global_num_nodes);
    const auto degree_bins                           = ComputeDegreeBins(edges, vertex_range);

    // Compute locality statistics
    const double edge_locality          = ComputeEdgeLocalicty(edges, vertex_range);
    const SInt   global_num_ghost_nodes = ComputeNumberOfGhostNodes(edges, vertex_range);
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
	std::cout << "  There are " << global_num_nodes << " real vertices and " << global_num_ghost_nodes << " ghost vertices" << std::endl;
    }
}
} // namespace kagen
