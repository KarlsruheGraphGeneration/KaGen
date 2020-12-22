/*******************************************************************************
 * include/generators/gnp/gnp_undirected.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _GNP_UNDIRECTED_H_
#define _GNP_UNDIRECTED_H_

#include <iostream>
#include <vector>

#include "definitions.h"
#include "generator_config.h"
#include "generator_io.h"
#include "rng_wrapper.h"
#include "hash.hpp"

namespace kagen {

template <typename EdgeCallback> 
class GNPUndirected {
 public:
  GNPUndirected(const PGeneratorConfig &config, const PEID /* rank */,
                const EdgeCallback &cb)
      : config_(config), rng_(config), io_(config), cb_(cb) { }

  void Generate() {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Chunk distribution
    nodes_per_chunk = config_.n / config_.k;
    SInt leftover_chunks = config_.k % size;
    SInt remaining_nodes = config_.n % config_.k;
    SInt nodes_per_chunk = config_.n / config_.k;
    SInt num_chunks = config_.k / size + ((SInt)rank < leftover_chunks);
    SInt start_chunk = rank * num_chunks +
                    ((SInt)rank >= leftover_chunks ? leftover_chunks : 0);
    SInt end_chunk = start_chunk + num_chunks;

    start_node_ = start_chunk * nodes_per_chunk + std::min(remaining_nodes, start_chunk);
    end_node_ = end_chunk * nodes_per_chunk + std::min(remaining_nodes, end_chunk);
    num_nodes_ = end_node_ - start_node_ - 1;

    // Generate chunks
    for (SInt i = 0; i < num_chunks; i++) {
      SInt row_n = 0;
      SInt column_n = 0;
      SInt row_node_id =
          start_chunk * nodes_per_chunk + std::min(start_chunk, remaining_nodes);
      SInt column_node_id = 0;
      SInt current_row = start_chunk++;
      SInt current_column = 0;
      // Iterate current_row
      while (current_column < current_row) {
        row_n = nodes_per_chunk + (current_row < remaining_nodes);
        column_n = nodes_per_chunk + (current_column < remaining_nodes);
        GenerateRectangleChunk(current_row, current_column++, row_node_id,
                               column_node_id, row_n, column_n);
        column_node_id += column_n;
      }
      // Handle triangular section
      if (current_row < config_.k) {
        row_n = nodes_per_chunk + (current_row < remaining_nodes);
        column_n = nodes_per_chunk + (current_column < remaining_nodes);
        // TODO: Triangle chunk
        GenerateTriangleChunk(current_row++, current_column,
                              row_node_id + (!config_.self_loops), column_node_id,
                              row_n, column_n);
        row_node_id += row_n;
      }
      // Iterate current_column
      while (current_row < config_.k) {
        row_n = nodes_per_chunk + (current_row < remaining_nodes);
        column_n = nodes_per_chunk + (current_column < remaining_nodes);
        GenerateRectangleChunk(current_row++, current_column, row_node_id,
                               column_node_id, row_n, column_n);
        row_node_id += row_n;
      }
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
  RNGWrapper rng_;

  // I/O
  GeneratorIO<> io_;
  EdgeCallback cb_;

  // Constants and variables
  SInt nodes_per_chunk;
  SInt start_node_, end_node_, num_nodes_;

  void GenerateTriangleChunk(const SInt row_id, const SInt column_id,
                             const SInt row_node_id, const SInt column_node_id,
                             const SInt row_n, const SInt column_n) {
    GenerateTriangularEdges(row_n, column_n, config_.p, row_id, column_id,
                            row_node_id, column_node_id);
  }

  void GenerateRectangleChunk(const SInt row_id, const SInt column_id,
                              const SInt row_node_id, const SInt column_node_id,
                              const SInt row_n, const SInt column_n) {
    GenerateRectangleEdges(row_n, column_n, config_.p, row_id, column_id,
                           row_node_id, column_node_id);
  }

  void GenerateTriangularEdges(const SInt row_n, const SInt column_n,
                               const double p, const SInt row_id,
                               const SInt column_id, const SInt offset_row,
                               const SInt offset_column) {
    // Number of edges
    HPFloat total_edges = 0;
    if (!config_.self_loops)
      total_edges = row_n * (column_n - 1) / 2;
    else
      total_edges = row_n * (column_n + 1) / 2;

    // Generate variate
    SInt h =
        sampling::Spooky::hash(config_.seed + (((row_id + 1) * row_id) / 2) + column_id);
    SInt num_edges = rng_.GenerateBinomial(h, total_edges, p);

    // Sample from [1, num_edges]
    rng_.GenerateSample(h, total_edges, num_edges, [&](SInt sample) {
      // Absolute triangular point
      // if (loops) sqr = (sqrt(8*((double)sample-1)+1) - 1)/2 + 1;
      SInt sqr = sqrt(8 * (sample - 1) + 1);
      // TODO: Nasty hack
      while (sqr * sqr > 8 * (sample - 1) + 1) sqr--;
      SInt i = (sqr - 1) / 2;
      SInt j = (sample - 1) - i * (i + 1) / 2;
      cb_(i + offset_row, j + offset_column);
#ifdef OUTPUT_EDGES
      io_.PushEdge(i + offset_row, j + offset_column);
#else
      io_.UpdateDist(i + offset_row);
      io_.UpdateDist(j + offset_column);
#endif
    });
  }

  void GenerateRectangleEdges(const SInt row_n, const SInt column_n,
                              const double p, const SInt row_id,
                              const SInt column_id, const SInt offset_row,
                              const SInt offset_column) {
    // Generate variate
    SInt h =
        sampling::Spooky::hash(config_.seed + (((row_id + 1) * row_id) / 2) + column_id);
    SInt num_edges = rng_.GenerateBinomial(h, row_n * column_n, p);

    // Sample from [1, num_edges]
    rng_.GenerateSample(h, row_n * column_n, num_edges, [&](SInt sample) {
                          SInt i = (sample - 1) / column_n;
                          SInt j = (sample - 1) % column_n;
                          cb_(i + offset_row, j + offset_column);
#ifdef OUTPUT_EDGES
                          io_.PushEdge(i + offset_row, j + offset_column);
#else
                          io_.UpdateDist(i + offset_row);
                          io_.UpdateDist(j + offset_column);
#endif
                        });
  }
};

}
#endif
