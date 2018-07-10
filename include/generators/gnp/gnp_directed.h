/*******************************************************************************
 * include/generators/gnp/gnp_directed.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _GNP_DIRECTED_H_
#define _GNP_DIRECTED_H_

#include <iostream>
#include <vector>

#include "definitions.h"
#include "generator_config.h"
#include "generator_io.h"
#include "rng_wrapper.h"
#include "hash.hpp"

namespace kagen {

template <typename EdgeCallback> 
class GNPDirected {
 public:
  GNPDirected(const PGeneratorConfig &config, const PEID /* rank */,
              const EdgeCallback &cb)
      : config_(config), rng_(config), io_(config), cb_(cb) {
    // Init variables
    if (!config_.self_loops)
      edges_per_node = config_.n - 1;
    else
      edges_per_node = config_.n;
  }

  void Generate() {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Chunk distribution
    SInt nodes_per_chunk = config_.n / config_.k;
    SInt leftover_chunks = config_.k % size;
    SInt remaining_nodes = config_.n % config_.k;
    SInt num_chunks = (config_.k / size) + ((SInt)rank < leftover_chunks);
    SInt start_chunk = rank * num_chunks +
                    ((SInt)rank >= leftover_chunks ? leftover_chunks : 0);
    SInt end_chunk = start_chunk + num_chunks;

    start_node_ = start_chunk * nodes_per_chunk + std::min(start_chunk, remaining_nodes);
    end_node_ = end_chunk * nodes_per_chunk + std::min(end_chunk, remaining_nodes);
    num_nodes_ = end_node_ - start_node_ - 1;

    // Generate chunks
    SInt current_node = start_node_;
    for (SInt i = 0; i < num_chunks; i++) {
      SInt nodes_for_chunk = nodes_per_chunk + (start_chunk < remaining_nodes);
      GenerateChunk(start_chunk++, current_node, nodes_for_chunk);
      current_node += nodes_for_chunk;
    }
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
  RNGWrapper<> rng_;

  // I/O
  GeneratorIO<> io_;
  EdgeCallback cb_;

  // Constants and variables
  SInt edges_per_node;
  SInt start_node_, end_node_, num_nodes_;

  void GenerateChunk(const SInt chunk_id, const SInt node_id, const SInt n) {
    GenerateEdges(n, config_.p, chunk_id, node_id);
  }

  void GenerateEdges(const SInt n, const double p, const SInt chunk_id,
                     const SInt offset) {
    // Generate variate
    SInt h = sampling::Spooky::hash(config_.seed + chunk_id);
    SInt num_edges = (SInt)rng_.GenerateBinomial(h, n * edges_per_node, p);

    // Sample from [1, num_edges]
    rng_.GenerateSample(h, (HPFloat)n * edges_per_node, num_edges,
                        [&](SInt sample) {
                          SInt source = (sample - 1) / edges_per_node + offset;
                          SInt target = (sample - 1) % edges_per_node;
                          if (!config_.self_loops)
                            target += ((sample - 1) % edges_per_node >= source);
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
