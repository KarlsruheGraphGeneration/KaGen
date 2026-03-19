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
    switch (config_.distribution) {
        case kagen::GraphDistribution::BALANCE_EDGES:
            graph_.vertex_range = RedistributeEdgesBalanced(local_edges_, graph_.edges, n, remap_round_robin, comm);
            break;
        case kagen::GraphDistribution::BALANCE_VERTICES:
            graph_.vertex_range = RedistributeEdges(local_edges_, graph_.edges, n, remap_round_robin, comm);
            break;
    }
}
} // namespace kagen
