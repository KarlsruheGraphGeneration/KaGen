#include "kagen/tools/statistics.h"

#include "kagen/definitions.h"
#include "kagen/kagen.h"
#include "kagen/tools/utils.h"

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

namespace {
// Reduces already-computed local min/sum/max degrees into global DegreeStatistics.
DegreeStatistics ReduceLocalDegreeStatistics(
    const SInt local_min, const SInt local_sum, const SInt local_max, const SInt global_num_nodes, MPI_Comm comm) {
    return {ReduceMin(local_min, comm), 1.0 * ReduceSum(local_sum, comm) / global_num_nodes, ReduceMax(local_max, comm)};
}
} // namespace

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

    return ReduceLocalDegreeStatistics(min, sum, max, global_num_nodes, comm);
}

DegreeStatistics ReduceDegreeStatistics(const XadjArray& xadj, const SInt global_num_nodes, MPI_Comm comm) {
    SInt min = std::numeric_limits<SInt>::max();
    SInt sum = 0;
    SInt max = std::numeric_limits<SInt>::lowest();

    for (std::size_t i = 0; i + 1 < xadj.size(); ++i) {
        const SInt deg = xadj[i + 1] - xadj[i];
        min            = std::min(min, deg);
        max            = std::max(max, deg);
        sum += deg;
    }

    return ReduceLocalDegreeStatistics(min, sum, max, global_num_nodes, comm);
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

std::vector<SInt> ComputeDegreeBins(const XadjArray& xadj, MPI_Comm comm) {
    std::vector<SInt> bins(std::numeric_limits<SInt>::digits);

    for (std::size_t i = 0; i + 1 < xadj.size(); ++i) {
        const SInt deg = xadj[i + 1] - xadj[i];
        const SInt bin = (deg == 0) ? 0 : (std::log2(deg) + 1);
        ++bins[bin];
    }

    std::vector<SInt> global_bins(bins.size());
    MPI_Reduce(bins.data(), global_bins.data(), bins.size(), KAGEN_MPI_SINT, MPI_SUM, ROOT, comm);

    return global_bins;
}

namespace {
// Shared tail for ComputeEdgeLocality: given how many of the local edge endpoints are cut edges (point outside
// vertex_range) and the local edge count, reduces both and turns them into a locality fraction.
double ReduceEdgeLocality(const SInt num_local_cut_edges, const SInt num_local_edges, MPI_Comm comm) {
    const SInt num_global_cut_edges = ReduceSum(num_local_cut_edges, comm);
    const SInt num_global_edges     = ReduceSum(num_local_edges, comm);
    return 1.0 - DivideOrDefault(static_cast<double>(num_global_cut_edges), static_cast<double>(num_global_edges), 0.0);
}

// Shared tail for ComputeNumberOfGhostNodes: reduces the count of distinct out-of-range endpoints found locally.
SInt ReduceNumberOfGhostNodes(const std::unordered_set<SInt>& ghost_nodes, MPI_Comm comm) {
    return ReduceSum(static_cast<SInt>(ghost_nodes.size()), comm);
}
} // namespace

double ComputeEdgeLocality(const Edgelist& edges, const VertexRange vertex_range, MPI_Comm comm) {
    const SInt num_local_cut_edges = std::count_if(edges.begin(), edges.end(), [&vertex_range](const auto& edge) {
        return std::get<1>(edge) < vertex_range.first || std::get<1>(edge) >= vertex_range.second;
    });
    return ReduceEdgeLocality(num_local_cut_edges, edges.size(), comm);
}

double ComputeEdgeLocality(const AdjncyArray& adjncy, const VertexRange vertex_range, MPI_Comm comm) {
    const SInt num_local_cut_edges = std::count_if(adjncy.begin(), adjncy.end(), [&vertex_range](const SInt to) {
        return to < vertex_range.first || to >= vertex_range.second;
    });
    return ReduceEdgeLocality(num_local_cut_edges, adjncy.size(), comm);
}

SInt ComputeNumberOfGhostNodes(const Edgelist& edges, const VertexRange vertex_range, MPI_Comm comm) {
    std::unordered_set<SInt> ghost_nodes;

    for (const auto& [from, to]: edges) {
        if (to < vertex_range.first || to >= vertex_range.second) {
            ghost_nodes.insert(to);
        }
    }

    return ReduceNumberOfGhostNodes(ghost_nodes, comm);
}

SInt ComputeNumberOfGhostNodes(const AdjncyArray& adjncy, const VertexRange vertex_range, MPI_Comm comm) {
    std::unordered_set<SInt> ghost_nodes;

    for (const SInt to: adjncy) {
        if (to < vertex_range.first || to >= vertex_range.second) {
            ghost_nodes.insert(to);
        }
    }

    return ReduceNumberOfGhostNodes(ghost_nodes, comm);
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
    const XadjArray&, const AdjncyArray& adjncy, const VertexRange vertex_range, const bool root, MPI_Comm comm) {
    // vertex_range.second - vertex_range.first, not xadj.size() - 1: the latter is the physically-present row
    // count, which double-counts a shared boundary vertex on a split graph (see Graph::PhysicalVertexRange()).
    PrintBasicStatistics(vertex_range.second - vertex_range.first, adjncy.size(), root, comm);
}

void PrintBasicStatistics(const Edgelist& edges, const VertexRange vertex_range, const bool root, MPI_Comm comm) {
    PrintBasicStatistics(vertex_range.second - vertex_range.first, edges.size(), root, comm);
}

namespace {
struct AdvancedStatistics {
    double             density;
    DegreeStatistics   degree;
    std::vector<SInt>  degree_bins;
    double             edge_locality;
    SInt               global_num_nodes;
    SInt               global_num_ghost_nodes;
};

void PrintAdvancedStatistics(const AdvancedStatistics& stats) {
    std::cout << "Density: " << std::fixed << std::setprecision(4) << stats.density << "\n";
    std::cout << "Degrees: [Min=" << stats.degree.min << " | Mean=" << std::fixed << std::setprecision(1)
               << stats.degree.mean << " | Max=" << stats.degree.max << "]\n";

    // Find last non-empty degree bin
    SInt last_nonempty_degree_bin = 0;
    for (SInt i = 0; i < stats.degree_bins.size(); ++i) {
        if (stats.degree_bins[i] > 0) {
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
                  << "): " << stats.degree_bins[i] << "\n";
    }

    // Print locality statistics
    const double ghost_node_fraction =
        1.0 * stats.global_num_ghost_nodes / (stats.global_num_nodes + stats.global_num_ghost_nodes);
    std::cout << "Edge locality: " << std::fixed << std::setprecision(4) << stats.edge_locality << std::endl;
    std::cout << "Fraction of ghost nodes: " << std::fixed << std::setprecision(4) << ghost_node_fraction
              << std::endl;
    std::cout << "  There are " << stats.global_num_nodes << " real vertices and " << stats.global_num_ghost_nodes
              << " ghost vertices" << std::endl;
}
} // namespace

void PrintAdvancedStatistics(
    const XadjArray& xadj, const AdjncyArray& adjncy, const VertexRange vertex_range, const bool root, MPI_Comm comm) {
    const auto local_num_nodes  = xadj.size() - 1;
    const auto global_num_nodes = ReduceSum(local_num_nodes, comm);
    const auto global_num_edges = ReduceSum(adjncy.size(), comm);

    const double density  = 1.0 * global_num_edges / global_num_nodes / (global_num_nodes - 1);
    const auto   degree   = ReduceDegreeStatistics(xadj, global_num_nodes, comm);
    const auto   degree_bins = ComputeDegreeBins(xadj, comm);

    const double edge_locality          = ComputeEdgeLocality(adjncy, vertex_range, comm);
    const SInt   global_num_ghost_nodes = ComputeNumberOfGhostNodes(adjncy, vertex_range, comm);

    if (root) {
        PrintAdvancedStatistics(
            AdvancedStatistics{density, degree, degree_bins, edge_locality, global_num_nodes, global_num_ghost_nodes});
    }
}

void PrintAdvancedStatistics(Edgelist& edges, const VertexRange vertex_range, const bool root, MPI_Comm comm) {
    // Sort edges for degree computation
    if (!std::is_sorted(edges.begin(), edges.end())) {
        std::sort(edges.begin(), edges.end());
    }

    // Compute degree statistics
    const auto local_num_nodes  = vertex_range.second - vertex_range.first;
    const auto global_num_nodes = ReduceSum(local_num_nodes, comm);
    const auto global_num_edges = ReduceSum(edges.size(), comm);

    const double density     = 1.0 * global_num_edges / global_num_nodes / (global_num_nodes - 1);
    const auto   degree      = ReduceDegreeStatistics(edges, global_num_nodes, comm);
    const auto   degree_bins = ComputeDegreeBins(edges, vertex_range, comm);

    // Compute locality statistics
    const double edge_locality          = ComputeEdgeLocality(edges, vertex_range, comm);
    const SInt   global_num_ghost_nodes = ComputeNumberOfGhostNodes(edges, vertex_range, comm);

    if (root) {
        PrintAdvancedStatistics(
            AdvancedStatistics{density, degree, degree_bins, edge_locality, global_num_nodes, global_num_ghost_nodes});
    }
}
} // namespace kagen
