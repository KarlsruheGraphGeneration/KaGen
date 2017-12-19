/******************************************************************************
 * gnp_directed.h
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

#ifndef _GNP_DIRECTED_H_
#define _GNP_DIRECTED_H_

#include <iostream>
#include <vector>

#include "definitions.h"
#include "generator_config.h"
#include "generator_io.h"
#include "rng_wrapper.h"
#include "hash.hpp"

class GNPDirected {
 public:
  GNPDirected(const PGeneratorConfig &config, const PEID rank)
      : config_(config), rng_(config), io_(config) {
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
    SInt leftover_nodes = config_.n % config_.k;
    SInt num_chunks = (config_.k / size) + ((SInt)rank < leftover_chunks);
    SInt chunk_id = rank * num_chunks +
                    ((SInt)rank >= leftover_chunks ? leftover_chunks : 0);
    SInt node_id =
        chunk_id * nodes_per_chunk + std::min(chunk_id, leftover_nodes);

    // Generate chunks
    for (SInt i = 0; i < num_chunks; i++) {
      SInt nodes_for_chunk = nodes_per_chunk + (chunk_id < leftover_nodes);
      GenerateChunk(chunk_id++, node_id, nodes_for_chunk);
      node_id += nodes_for_chunk;
    }
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
  SInt edges_per_node;

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
