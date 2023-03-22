#include "kagen/generators/path/path_directed.h"
#include "kagen/tools/random_permutation.h"

namespace kagen {

std::unique_ptr<Generator>
PathDirectedFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<PathDirected>(config, rank, size);
}

PathDirected::PathDirected(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size) {}

std::pair<SInt, bool> PathDirected::GetTargetVertex(SInt i) const {
    if (!config_.permute) {
        const SInt j = (i + 1) % config_.n;
        return std::make_pair(j, j != 0 || config_.periodic);
    } else {
        static LCGPseudoRandomPermutation permutator(config_.n - 1);
        SInt                              j = (permutator.finv(i) + 1) % config_.n;
        std::cout << "i: " << i << " permutator.finv(): " << permutator.finv(i) + 1 << std::endl;
        return std::make_pair(permutator.f(j), j != 0 || config_.periodic);
    }
}

void PathDirected::GenerateEdgeList() {
    if (config_.n <= 1) {
        return;
    }

    // all PE get n / size nodes
    // the first n modulo size PEs obtain one additional node.
    const SInt nodes_per_pe                = config_.n / size_;
    const PEID num_pe_with_additional_node = config_.n % size_;
    const bool has_pe_additional_node      = rank_ < num_pe_with_additional_node;
    const SInt begin_nodes                 = std::min(num_pe_with_additional_node, rank_) + rank_ * nodes_per_pe;
    const SInt end_nodes                   = begin_nodes + nodes_per_pe + has_pe_additional_node;


    for (SInt i = begin_nodes; i < end_nodes; ++i) {
        const auto [j, is_valid] = GetTargetVertex(i);
        if (is_valid) {
            PushEdge(i, j);
        }
    }
}

} // namespace kagen
