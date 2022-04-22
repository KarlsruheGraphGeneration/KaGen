#include "kagen/tools/statistics.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

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
} // namespace kagen
