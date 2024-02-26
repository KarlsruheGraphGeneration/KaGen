#include "kagen/generators/gnm/gnm_directed.h"

#include "kagen/sampling/hash.hpp"

namespace kagen {
std::unique_ptr<Generator>
GNMDirectedFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<GNMDirected>(config, rank, size);
}

GNMDirected::GNMDirected(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size),
      rng_(config),
      edges_per_node_(config_.self_loops ? config_.n : config_.n - 1) {}

void GNMDirected::GenerateEdgeList() {
    // Chunk distribution
    SInt leftover_chunks = config_.k % size_;
    SInt num_chunks      = (config_.k / size_) + ((SInt)rank_ < leftover_chunks);
    SInt start_chunk     = rank_ * num_chunks + ((SInt)rank_ >= leftover_chunks ? leftover_chunks : 0);
    SInt end_chunk       = start_chunk + num_chunks;

    SInt nodes_per_chunk = config_.n / config_.k;
    SInt remaining_nodes = config_.n % config_.k;

    start_node_ = start_chunk * nodes_per_chunk + std::min(remaining_nodes, start_chunk);
    end_node_   = end_chunk * nodes_per_chunk + std::min(remaining_nodes, end_chunk);
    num_nodes_  = end_node_ - start_node_ - 1;

    // Generate chunks
    for (SInt i = 0; i < num_chunks; i++)
        GenerateChunk(start_chunk++);

    SetVertexRange(start_node_, start_node_ + num_nodes_);
}

void GNMDirected::GenerateChunk(const SInt chunk_id) {
    GenerateChunk(config_.n, config_.m, config_.k, chunk_id, 0, 0, 1);
}

void GNMDirected::GenerateChunk(
    const SInt n, const SInt m, const SInt k, const SInt chunk_id, const SInt chunk_start, const SInt node_start,
    const SInt level) {
    // Stop if there are no edges left
    if (m <= 0)
        return;

    // Base Case if only one chunk is left
    if (k == 1) {
        GenerateEdges(n, m, chunk_id, node_start);
        return;
    }

    // Determine splitter and hypergeometric random variate
    SInt k_split = (k + 1) / 2;
    SInt n_split = k_split * (n / k) + std::min(n % k, k_split);

    // Generate variate
    SInt h       = sampling::Spooky::hash(config_.seed + level * config_.n + chunk_start);
    SInt variate = rng_.GenerateHypergeometric(h, n_split * edges_per_node_, m, n * edges_per_node_);

    // Distributed splitting of chunks
    if (chunk_id < chunk_start + k_split) {
        GenerateChunk(n_split, variate, (k + 1) / 2, chunk_id, chunk_start, node_start, level + 1);
    } else {
        GenerateChunk(
            n - n_split, m - variate, k / 2, chunk_id, chunk_start + k_split, node_start + n_split, level + 1);
    }
}

void GNMDirected::GenerateEdges(const SInt n, const SInt m, const SInt chunk_id, const SInt offset) {
    // Sample from [1, num_edges]
    SInt h = sampling::Spooky::hash(config_.seed + chunk_id);
    rng_.GenerateSample(h, n * edges_per_node_, m, [&](SInt sample) {
        SInt source = (sample - 1) / edges_per_node_ + offset;
        SInt target = (sample - 1) % edges_per_node_;
        if (!config_.self_loops)
            target += (((sample - 1) % edges_per_node_) >= source);
        PushEdge(source, target);
    });
}
} // namespace kagen
