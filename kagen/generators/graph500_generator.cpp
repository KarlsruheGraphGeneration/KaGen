#include "kagen/generators/graph500_generator.h"

#include "kagen/tools/postprocessor.h"
#include "kagen/tools/utils.h"

#include <mpi.h>

#include <cmath>

namespace kagen {
Graph500Generator::Graph500Generator(const PGeneratorConfig& config) : config_(config) {}

void Graph500Generator::FinalizeEdgeList(MPI_Comm comm) {
    const SInt log_n = std::log2(config_.n);
    const SInt n     = 1ull << log_n;

    const bool remap_round_robin = true;
    switch (config_.redistribution) {
        case kagen::GraphRedistribution::BALANCE_EDGES:
            graph_.vertex_range = RedistributeEdgesBalanced(local_edges_, graph_.edges, n, remap_round_robin, comm);
            break;
        case kagen::GraphRedistribution::BALANCE_VERTICES:
            graph_.vertex_range = RedistributeEdges(local_edges_, graph_.edges, n, remap_round_robin, comm);
            break;
        case kagen::GraphRedistribution::BALANCE_EDGES_TRUE: {
            const EdgeBalancedDistribution distribution =
                RedistributeEdgesTrueBalance(local_edges_, graph_.edges, n, remap_round_robin, comm);
            graph_.vertex_range = distribution.vertex_range;
            SetHasSplitVertices(distribution.has_split_vertices);
            SetPartialVertices(distribution.left_partial_vertex, distribution.right_partial_vertex);
            break;
        }
    }
}
} // namespace kagen
