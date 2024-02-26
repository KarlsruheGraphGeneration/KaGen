#include "kagen/generators/gnp/gnp_directed.h"

#include "kagen/sampling/hash.hpp"

namespace kagen {
std::unique_ptr<Generator>
GNPDirectedFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<GNPDirected>(config, rank, size);
}

GNPDirected::GNPDirected(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size),
      rng_(config) {
    // Init variables
    if (!config_.self_loops)
        edges_per_node = config_.n - 1;
    else
        edges_per_node = config_.n;
}

void GNPDirected::GenerateEdgeList() {
    // Chunk distribution
    SInt nodes_per_chunk = config_.n / config_.k;
    SInt leftover_chunks = config_.k % size_;
    SInt remaining_nodes = config_.n % config_.k;
    SInt num_chunks      = (config_.k / size_) + ((SInt)rank_ < leftover_chunks);
    SInt start_chunk     = rank_ * num_chunks + ((SInt)rank_ >= leftover_chunks ? leftover_chunks : 0);
    SInt end_chunk       = start_chunk + num_chunks;

    start_node_ = start_chunk * nodes_per_chunk + std::min(start_chunk, remaining_nodes);
    end_node_   = end_chunk * nodes_per_chunk + std::min(end_chunk, remaining_nodes);
    num_nodes_  = end_node_ - start_node_ - 1;

    // Generate chunks
    SInt current_node = start_node_;
    for (SInt i = 0; i < num_chunks; i++) {
        SInt nodes_for_chunk = nodes_per_chunk + (start_chunk < remaining_nodes);
        GenerateChunk(start_chunk++, current_node, nodes_for_chunk);
        current_node += nodes_for_chunk;
    }

    SetVertexRange(start_node_, start_node_ + num_nodes_);
}

void GNPDirected::GenerateChunk(const SInt chunk_id, const SInt node_id, const SInt n) {
    GenerateEdges(n, config_.p, chunk_id, node_id);
}

void GNPDirected::GenerateEdges(const SInt n, const double p, const SInt chunk_id, const SInt offset) {
    // Generate variate
    SInt h         = sampling::Spooky::hash(config_.seed + chunk_id);
    SInt num_edges = rng_.GenerateBinomial(h, n * edges_per_node, p);

    // Sample from [1, num_edges]
    rng_.GenerateSample(h, n * edges_per_node, num_edges, [&](SInt sample) {
        SInt source = (sample - 1) / edges_per_node + offset;
        SInt target = (sample - 1) % edges_per_node;
        if (!config_.self_loops)
            target += ((sample - 1) % edges_per_node >= source);
        PushEdge(source, target);
    });
}
} // namespace kagen
