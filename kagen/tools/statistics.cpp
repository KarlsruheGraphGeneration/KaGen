#include "kagen/tools/statistics.h"

#include <cmath>

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
} // namespace kagen
