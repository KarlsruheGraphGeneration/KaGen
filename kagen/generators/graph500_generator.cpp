#include "kagen/generators/graph500_generator.h"

#include "kagen/tools/postprocessor.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

namespace kagen {
Graph500Generator::Graph500Generator(const PGeneratorConfig& config) : config_(config) {}

void Graph500Generator::FinalizeEdgeList(MPI_Comm comm) {
    const SInt log_n = std::log2(config_.n);
    const SInt n     = 1ull << log_n;
    vertex_range_    = RedistributeEdgesRoundRobin(local_edges_, edges_, n, comm);
}
} // namespace kagen
