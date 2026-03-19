#include "kagen/generators/graph500_generator.h"

#include "kagen/tools/postprocessor.h"

#include <mpi.h>

#include <cmath>

namespace kagen {
Graph500Generator::Graph500Generator(const PGeneratorConfig& config) : config_(config) {}

void Graph500Generator::FinalizeEdgeList(MPI_Comm comm) {
    const SInt log_n = std::log2(config_.n);
    const SInt n     = 1ull << log_n;
    switch (config_.distribution) {
        case kagen::GraphDistribution::BALANCE_EDGES: {
            const bool remap_round_robin = false;
            graph_.vertex_range = RedistributeEdgesBalanced(local_edges_, graph_.edges, n, remap_round_robin, comm);
            std::cout << "edge balanced" << std::endl;
            break;
        }
        case kagen::GraphDistribution::BALANCE_VERTICES:
            graph_.vertex_range = RedistributeEdgesRoundRobin(local_edges_, graph_.edges, n, comm);
            std::cout << "vertex balanced" << std::endl;
            break;
        case kagen::GraphDistribution::EXPLICIT:
            throw std::runtime_error("not supported");
            break;
        case kagen::GraphDistribution::ROOT:
            throw std::runtime_error("not supported");
            break;
        default:
            graph_.vertex_range = RedistributeEdgesRoundRobin(local_edges_, graph_.edges, n, comm);
            break;
    }
}
} // namespace kagen
