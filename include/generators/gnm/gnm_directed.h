/*******************************************************************************
 * include/generators/gnm/gnm_directed.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _GNM_DIRECTED_H_
#define _GNM_DIRECTED_H_

#include <iostream>
#include <vector>

#include "definitions.h"
#include "generator_config.h"
#include "generator_io.h"
#include "rng_wrapper.h"
#include "hash.hpp"

namespace kagen {

template <typename EdgeCallback> 
class GNMDirected {
 public:
  GNMDirected(const PGeneratorConfig &config, const PEID /* rank */,
              const EdgeCallback &cb)
      : config_(config), rng_(config), io_(config), cb_(cb) {
    // Init variables
    if (!config_.self_loops)
      edges_per_node_ = config_.n - 1;
    else
      edges_per_node_ = config_.n;
  }

  void Generate() {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Chunk distribution
    SInt leftover_chunks = config_.k % size;
    SInt num_chunks = (config_.k / size) + ((SInt)rank < leftover_chunks);
    SInt start_chunk = rank * num_chunks +
                         ((SInt)rank >= leftover_chunks ? leftover_chunks : 0);
    SInt end_chunk = start_chunk + num_chunks;

    SInt nodes_per_chunk = config_.n / config_.k;
    SInt remaining_nodes = config_.n % config_.k;

    start_node_ = start_chunk * nodes_per_chunk + std::min(remaining_nodes, start_chunk);
    end_node_ = end_chunk * nodes_per_chunk + std::min(remaining_nodes, end_chunk);
    num_nodes_ = end_node_ - start_node_ - 1;

    // Generate chunks
    for (SInt i = 0; i < num_chunks; i++) GenerateChunk(start_chunk++);
  }

  void Output() const { 
#ifdef OUTPUT_EDGES
    io_.OutputEdges(); 
#else
    io_.OutputDist(); 
#endif
  }

  std::pair<SInt, SInt> GetVertexRange() {
    return std::make_pair(start_node_, start_node_ + num_nodes_);
  }

  SInt NumberOfEdges() const { return io_.NumEdges(); }

 private:
  // Config
  PGeneratorConfig config_;

  // Variates
  RNGWrapper rng_;

  // I/O
  GeneratorIO<> io_;
  EdgeCallback cb_; 

  // Constants and variables
  SInt edges_per_node_;
  SInt start_node_, end_node_, num_nodes_;

  void GenerateChunk(const SInt chunk_id) {
    GenerateChunk(config_.n, config_.m, config_.k, chunk_id, 0, 0, 1);
  }

  void GenerateChunk(const SInt n, const SInt m, const SInt k,
                     const SInt chunk_id, const SInt chunk_start,
                     const SInt node_start, const SInt level) {
    // Stop if there are no edges left
    if (m <= 0) return;

    // Base Case if only one chunk is left
    if (k == 1) {
      GenerateEdges(n, m, chunk_id, node_start);
      return;
    }

    // Determine splitter and hypergeometric random variate
    SInt k_split = (k + 1) / 2;
    SInt n_split = k_split * (n / k) + std::min(n % k, k_split);

    // Generate variate
    SInt h = sampling::Spooky::hash(config_.seed + level * config_.n + chunk_start);
    SInt variate = rng_.GenerateHypergeometric(h, n_split * edges_per_node_, m, n * edges_per_node_);

    // Distributed splitting of chunks
    if (chunk_id < chunk_start + k_split) {
      GenerateChunk(n_split, variate, (k + 1) / 2, chunk_id, chunk_start,
                    node_start, level + 1);
    } else {
      GenerateChunk(n - n_split, m - variate, k / 2, chunk_id,
                    chunk_start + k_split, node_start + n_split, level + 1);
    }
  }

  void GenerateEdges(const SInt n, const SInt m, const SInt chunk_id,
                     const SInt offset) {
    // Sample from [1, num_edges]
    SInt h = sampling::Spooky::hash(config_.seed + chunk_id);
    rng_.GenerateSample(h, n * edges_per_node_, m, [&](SInt sample) {
      SInt source = (sample - 1) / edges_per_node_ + offset;
      SInt target = (sample - 1) % edges_per_node_;
      if (!config_.self_loops)
        target += (((sample - 1) % edges_per_node_) >= source);
      cb_(source, target);
#ifdef OUTPUT_EDGES
      io_.PushEdge(source, target);
#else
      io_.UpdateDist(source);
      io_.UpdateDist(target);
#endif
    });
  }
};

}
#endif
