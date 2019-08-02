/*******************************************************************************
 * include/generators/grid/grid_2d.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _GRID_2D_H_
#define _GRID_2D_H_

#include <iostream>
#include <vector>

#include "morton2D.h"
#include "definitions.h"
#include "generator_config.h"
#include "generator_io.h"
#include "rng_wrapper.h"
#include "hash.hpp"

namespace kagen {

template <typename EdgeCallback> 
class Grid2D {
 public:
  Grid2D(const PGeneratorConfig &config, const PEID /* rank */,
       const EdgeCallback &cb)
      : config_(config), rng_(config), io_(config), cb_(cb) { }

  void Generate() {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Init dimensions
    // TODO: Only tested for powers of two
    total_rows_ = log2(config_.n);
    total_cols_ = log2(config_.m);
    edge_probability_ = config_.p;

    // Init chunks
    total_chunks_ = config_.k;
    chunks_per_dim_ = sqrt(total_chunks_);

    SInt leftover_chunks = total_chunks_ % size;
    SInt num_chunks = (total_chunks_ / size) + ((SInt)rank < leftover_chunks);
    SInt start_chunk = rank * num_chunks +
                         ((SInt)rank >= leftover_chunks ? leftover_chunks : 0);
    SInt end_chunk = start_chunk + num_chunks;

    // Chunk distribution
    rows_per_chunk_ = total_rows_ / chunks_per_dim_;
    remaining_rows_ = total_rows_ % chunks_per_dim_;

    cols_per_chunk_ = total_cols_ / chunks_per_dim_;
    remaining_cols_ = total_cols_ % chunks_per_dim_;

    // Determine node count
    SInt start_chunk_row, start_chunk_col;
    Decode(start_chunk, start_chunk_col, start_chunk_row);
    SInt start_vertex_row = start_chunk_row * rows_per_chunk_ + std::min(remaining_rows_, start_chunk_row);
    SInt start_vertex_col = start_chunk_col * cols_per_chunk_ + std::min(remaining_cols_, start_chunk_col);

    SInt end_chunk_row, end_chunk_col;
    Decode(end_chunk, end_chunk_col, end_chunk_row);

    SInt end_vertex_row = end_chunk_row * rows_per_chunk_ + std::min(remaining_rows_, end_chunk_row);
    SInt end_vertex_col = end_chunk_col * cols_per_chunk_ + std::min(remaining_cols_, end_chunk_col);

    start_node_ = Encode(start_vertex_col, start_vertex_row);
    end_node_ = Encode(end_vertex_col, end_vertex_row);
    num_nodes_ = end_node_ - start_node_;

    // std::cout << "R" << rank << " chunks [" << start_chunk << "," << end_chunk << ")" << std::endl;
    // std::cout << "R" << rank << " nodes [" << start_node_ << "," << end_node_ << ") num " << num_nodes_ << std::endl;
    for (SInt i = 0; i < num_chunks; i++) {
      GenerateChunk(start_chunk + i);
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
    return std::make_pair(start_node_, start_node_ + num_nodes_ - 1);
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
  SInt start_node_, end_node_, num_nodes_;
  LPFloat edge_probability_;
  SInt total_rows_, total_cols_;
  SInt total_chunks_, chunks_per_dim_;
  SInt rows_per_chunk_, cols_per_chunk_;
  SInt remaining_rows_, remaining_cols_;

  void GenerateChunk(const SInt chunk) {
    // Compute position of chunk
    SInt chunk_row, chunk_col;
    Decode(chunk, chunk_col, chunk_row);

    // Compute start vertex
    SInt vertex_row = chunk_row * rows_per_chunk_ + std::min(remaining_rows_, chunk_row);
    SInt vertex_col = chunk_col * cols_per_chunk_ + std::min(remaining_cols_, chunk_col);

    SInt num_rows = rows_per_chunk_ + (chunk_row < remaining_rows_);
    SInt num_cols = cols_per_chunk_ + (chunk_col < remaining_cols_);

    for (SInt i = 0; i < num_rows; i++) {
      for (SInt j = 0; j < num_cols; j++) {
        SInt row = vertex_row + i;
        SInt col = vertex_col + j;
        GenerateEdges(row, col);
      }
    }
  }

  void GenerateEdges(const SInt row, const SInt col) {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (SSInt r = -1; r <= 1; r++) {
      for (SSInt c = -1; c <= 1; c++) {
        if (r != 0 && c != 0) continue;
        if (r == 0 && c == 0) continue;
        SSInt neighbor_row = (SSInt)row + r;
        SSInt neighbor_col = (SSInt)col + c;

        // Check if neighbor is valid 
        if (config_.periodic) {
          neighbor_row = (neighbor_row + total_rows_) % total_rows_;
          neighbor_col = (neighbor_col + total_cols_) % total_cols_;
        } 
        if (neighbor_row < 0 || neighbor_row >= total_rows_) continue;
        if (neighbor_col < 0 || neighbor_col >= total_cols_) continue;

        SInt source = Encode(col, row);
        SInt target = Encode((SInt)neighbor_col, (SInt)neighbor_row);

        SInt edge_seed = std::min(source, target) * total_rows_ * total_cols_ + std::max(source, target);
        SInt h = sampling::Spooky::hash(config_.seed + edge_seed);
        if (rng_.GenerateBinomial(h, 1, edge_probability_)) {
          // std::cout << "R" << rank << " e " << source << " -> " << target << std::endl;
          cb_(source, target);
          cb_(target, source);
#ifdef OUTPUT_EDGES
          io_.PushEdge(source, target);
#else
          io_.UpdateDist(source);
          io_.UpdateDist(target);
#endif
        }
      }
    }
  }

  inline void Decode(const SInt id, SInt &x, SInt &y) const {
    m2D_d_sLUT(id, y, x);
  }

  // Chunk/vertex coding
  inline SInt Encode(const SInt x, const SInt y) const {
    return m2D_e_sLUT<SInt>(y, x);
  }
};

}
#endif
