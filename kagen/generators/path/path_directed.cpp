#include "kagen/generators/path/path_directed.h"

#ifdef KAGEN_XXHASH_FOUND
    #include "kagen/tools/random_permutation.h"
#endif // KAGEN_XXHASH_FOUND

namespace kagen {

std::unique_ptr<Generator>
PathDirectedFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<PathDirected>(config, rank, size);
}

PGeneratorConfig PathDirectedFactory::NormalizeParameters(PGeneratorConfig config, PEID, PEID, bool) const {
#ifndef KAGEN_XXHASH_FOUND
    if (config.permute) {
        throw ConfigurationError("path permutation requires xxHash, but build was configured without xxHash");
    }
#endif // KAGEN_XXHASH_FOUND

    return config;
}

PathDirected::PathDirected(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size) {}

void PathDirected::GenerateEdgeList() {
    if (config_.n <= 1) {
        return;
    }

#ifdef KAGEN_XXHASH_FOUND
    auto permutator = random_permutation::FeistelPseudoRandomPermutation::buildPermutation(config_.n - 1, 0);
#endif // KAGEN_XXHASH_FOUND

    // all PE get n / size nodes
    // the first n modulo size PEs obtain one additional node.
    const SInt nodes_per_pe                = config_.n / size_;
    const PEID num_pe_with_additional_node = config_.n % size_;
    const bool has_pe_additional_node      = rank_ < num_pe_with_additional_node;
    const SInt begin_nodes                 = std::min(num_pe_with_additional_node, rank_) + rank_ * nodes_per_pe;
    const SInt end_nodes                   = begin_nodes + nodes_per_pe + has_pe_additional_node;
    // this is an exclusive range
    SetVertexRange(begin_nodes, end_nodes);

    for (SInt i = begin_nodes; i < end_nodes; ++i) {
        const auto [j, is_valid] = [&]() -> std::pair<SInt, bool> {
#ifdef KAGEN_XXHASH_FOUND
            if (config_.permute) {
                SInt j = (permutator.finv(i) + 1) % config_.n;
                return std::make_pair(permutator.f(j), j != 0 || config_.periodic);
            }
#endif // KAGEN_XXHASH_FOUND

            const SInt j = (i + 1) % config_.n;
            return std::make_pair(j, j != 0 || config_.periodic);
        }();

        if (is_valid) {
            PushEdge(i, j);
        }
    }
}

} // namespace kagen
