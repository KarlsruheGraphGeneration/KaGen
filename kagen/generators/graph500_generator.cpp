#include "kagen/generators/graph500_generator.h"

#include "kagen/tools/postprocessor.h"

#include <mpi.h>

#include <cmath>

namespace kagen {
Graph500Generator::Graph500Generator(const PGeneratorConfig& config) : config_(config) {}

void Graph500Generator::FinalizeEdgeList(MPI_Comm comm) {
    const SInt log_n    = std::log2(config_.n);
    const SInt n        = 1ull << log_n;
    graph_.vertex_range = RedistributeEdgesRoundRobin(local_edges_, graph_.edges, n, comm);
}
} // namespace kagen
