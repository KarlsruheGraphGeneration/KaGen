/******************************************************************************
 * gnm_directed.h
 *
 * Source of the graph generator
 ******************************************************************************
 * Copyright (C) 2016 Sebastian Lamm <lamm@ira.uka.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef _GNM_DIRECTED_H_
#define _GNM_DIRECTED_H_

#include <iostream>
#include <vector>

#include "definitions.h"
#include "generator_config.h"
#include "generator_io.h"
#include "rng_wrapper.h"
#include "hash.hpp"

class GNMDirected {
 public:
  GNMDirected(const PGeneratorConfig &config, const PEID rank)
      : config_(config), rng_(config), io_(config) {
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
    SInt current_chunk = rank * num_chunks +
                         ((SInt)rank >= leftover_chunks ? leftover_chunks : 0);

    // Generate chunks
    for (SInt i = 0; i < num_chunks; i++) GenerateChunk(current_chunk++);
  }

  void Output() const { 
#ifdef OUTPUT_EDGES
    io_.OutputEdges(); 
#else
    io_.OutputDist(); 
#endif
  }

  SInt NumberOfEdges() const { return io_.NumEdges(); }

 private:
  // Config
  PGeneratorConfig config_;

  // Variates
  RNGWrapper<> rng_;

  // I/O
  GeneratorIO<> io_;

  // Constants and variables
  SInt edges_per_node_;

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
    SInt variate = (SInt)rng_.GenerateHypergeometric(
        h, (HPFloat)n_split * edges_per_node_, m, (HPFloat)n * edges_per_node_);

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
    rng_.GenerateSample(h, (HPFloat)n * edges_per_node_, m, [&](SInt sample) {
      SInt source = (sample - 1) / edges_per_node_ + offset;
      SInt target = (sample - 1) % edges_per_node_;
      if (!config_.self_loops)
        target += (((sample - 1) % edges_per_node_) >= source);
#ifdef OUTPUT_EDGES
      io_.PushEdge(source, target);
#else
      io_.UpdateDist(source);
      io_.UpdateDist(target);
#endif
    });
  }
};

#endif
