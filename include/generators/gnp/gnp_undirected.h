/******************************************************************************
 * gnp_undirected.h
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

#ifndef _GNP_UNDIRECTED_H_
#define _GNP_UNDIRECTED_H_

#include <iostream>
#include <vector>

#include "definitions.h"
#include "generator_config.h"
#include "generator_io.h"
#include "rng_wrapper.h"
#include "hash.hpp"

class GNPUndirected {
 public:
  GNPUndirected(const PGeneratorConfig &config, const PEID rank)
      : config_(config), rng_(config), io_(config) { }

  void Generate() {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Chunk distribution
    nodes_per_chunk = config_.n / config_.k;
    SInt leftover_chunks = config_.k % size;
    SInt leftover_nodes = config_.n % config_.k;
    SInt nodes_per_chunk = config_.n / config_.k;
    SInt num_chunks = config_.k / size + ((SInt)rank < leftover_chunks);
    SInt init_row = rank * num_chunks +
                    ((SInt)rank >= leftover_chunks ? leftover_chunks : 0);

    // Generate chunks
    for (SInt i = 0; i < num_chunks; i++) {
      SInt row_n = 0;
      SInt column_n = 0;
      SInt row_node_id =
          init_row * nodes_per_chunk + std::min(init_row, leftover_nodes);
      SInt column_node_id = 0;
      SInt current_row = init_row++;
      SInt current_column = 0;
      // Iterate current_row
      while (current_column < current_row) {
        row_n = nodes_per_chunk + (current_row < leftover_nodes);
        column_n = nodes_per_chunk + (current_column < leftover_nodes);
        GenerateRectangleChunk(current_row, current_column++, row_node_id,
                               column_node_id, row_n, column_n);
        column_node_id += column_n;
      }
      // Handle triangular section
      if (current_row < config_.k) {
        row_n = nodes_per_chunk + (current_row < leftover_nodes);
        column_n = nodes_per_chunk + (current_column < leftover_nodes);
        GenerateTriangleChunk(current_row++, current_column,
                              row_node_id + config_.self_loops, column_node_id,
                              row_n, column_n);
        row_node_id += row_n;
      }
      // Iterate current_column
      while (current_row < config_.k) {
        row_n = nodes_per_chunk + (current_row < leftover_nodes);
        column_n = nodes_per_chunk + (current_column < leftover_nodes);
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

  SInt NumberOfEdges() const { return io_.NumEdges(); }

 private:
  // Config
  PGeneratorConfig config_;

  // Variates
  RNGWrapper<> rng_;

  // I/O
  GeneratorIO<> io_;

  // Constants and variables
  SInt nodes_per_chunk;

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
    SInt num_edges = (SInt)rng_.GenerateBinomial(h, total_edges, p);

    // Sample from [1, num_edges]
    rng_.GenerateSample(h, total_edges, num_edges, [&](SInt sample) {
      // Absolute triangular point
      // if (loops) sqr = (sqrt(8*((double)sample-1)+1) - 1)/2 + 1;
      SInt sqr = sqrt(8 * (sample - 1) + 1);
      // TODO: Nasty hack
      while (sqr * sqr > 8 * (sample - 1) + 1) sqr--;
      SInt i = (sqr - 1) / 2;
      SInt j = (sample - 1) - i * (i + 1) / 2;
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
    SInt num_edges = (SInt)rng_.GenerateBinomial(h, row_n * column_n, p);

    // Sample from [1, num_edges]
    rng_.GenerateSample(h, (HPFloat)row_n * column_n, num_edges,
                        [&](SInt sample) {
                          SInt i = (sample - 1) / column_n;
                          SInt j = (sample - 1) % column_n;
#ifdef OUTPUT_EDGES
                          io_.PushEdge(i + offset_row, j + offset_column);
#else
                          io_.UpdateDist(i + offset_row);
                          io_.UpdateDist(j + offset_column);
#endif
                        });
  }
};

#endif
