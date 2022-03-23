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

// #include "morton2D.h"
#include "definitions.h"
#include "generator_config.h"
#include "generator_io.h"
#include "rng_wrapper.h"
#include "hash.hpp"

namespace kagen {

template <typename EdgeCallback> 
class Grid2D {
 public:
  Grid2D(PGeneratorConfig &config, const PEID /* rank */,
       const EdgeCallback &cb)
      : config_(config), rng_(config), io_(config), cb_(cb) { }

  void Generate() {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Init dimensions
    // TODO: Only tested for cube PEs and one chunk per PE
    total_rows_ = config_.grid_x;
    total_cols_ = config_.grid_y;
    config_.n = total_rows_ * total_cols_;
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

    vertices_per_chunk_ = rows_per_chunk_ * cols_per_chunk_;

    start_node_ = OffsetForChunk(start_chunk);
    end_node_ = OffsetForChunk(end_chunk);
    num_nodes_ = end_node_ - start_node_;

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
  PGeneratorConfig &config_;

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
  SInt vertices_per_chunk_;

  void GenerateChunk(const SInt chunk) {
    SInt start_vertex = OffsetForChunk(chunk);
    SInt end_vertex = OffsetForChunk(chunk + 1);
    for (SInt i = start_vertex; i < end_vertex; i++) {
      GenerateEdges(chunk, i);
    }
  }

  void GenerateEdges(const SInt chunk, const SInt vertex) {
    QueryInDirection(chunk, vertex, Direction::Right);
    QueryInDirection(chunk, vertex, Direction::Left);
    QueryInDirection(chunk, vertex, Direction::Up);
    QueryInDirection(chunk, vertex, Direction::Down);
  }

  void QueryInDirection(const SInt chunk, const SInt vertex, Direction direction) {
    SInt offset = OffsetForChunk(chunk);
    SInt local_vertex = vertex - offset;

    SInt chunk_row, chunk_col;
    Decode(chunk, chunk_row, chunk_col);

    SInt rows = rows_per_chunk_ + (chunk_row < remaining_rows_);
    SInt cols = cols_per_chunk_ + (chunk_col < remaining_cols_);

    SInt local_row = local_vertex / cols;
    SInt local_col = local_vertex % cols;

    SInt local_neighbor_row = local_row + DirectionRow(direction);
    SInt local_neighbor_col = local_col + DirectionColumn(direction);

    if (IsLocalVertex(local_neighbor_row, local_neighbor_col, rows, cols)) {
      SInt neighbor_vertex = offset + (local_neighbor_row * cols + local_neighbor_col);
      GenerateEdge(vertex, neighbor_vertex);
    } else {
      // Determine neighboring chunk
      SSInt neighbor_chunk_row = (SSInt)chunk_row + DirectionRow(direction);
      SSInt neighbor_chunk_col = (SSInt)chunk_col + DirectionColumn(direction);
      if (config_.periodic) {
        neighbor_chunk_row = (neighbor_chunk_row + chunks_per_dim_) % chunks_per_dim_;
        neighbor_chunk_col = (neighbor_chunk_col + chunks_per_dim_) % chunks_per_dim_;
      }
      if (!IsValidChunk(neighbor_chunk_row, neighbor_chunk_col)) return;

      SInt neighbor_chunk = Encode(neighbor_chunk_row, neighbor_chunk_col);
      SInt neighbor_vertex = LocateVertexInChunk(neighbor_chunk, local_row, local_col, direction);
      GenerateEdge(vertex, neighbor_vertex);
    }
  }

  bool IsLocalVertex(const SInt local_row, const SInt local_col, 
                     const SInt rows, const SInt cols) {
    if (local_row >= rows) return false;
    if (local_col >= cols) return false;
    return true;
  }

  bool IsValidChunk(const SInt chunk_row, const SInt chunk_col) {
    if (chunk_row >= chunks_per_dim_) return false;
    if (chunk_col >= chunks_per_dim_) return false;
    return true;
  }

  SInt LocateVertexInChunk(const SInt chunk, const SInt local_row, const SInt local_col, Direction direction) {
    SInt offset = OffsetForChunk(chunk);

    SInt chunk_row, chunk_col;
    Decode(chunk, chunk_row, chunk_col);
      
    SInt rows = rows_per_chunk_ + (chunk_row < remaining_rows_);
    SInt cols = cols_per_chunk_ + (chunk_col < remaining_cols_);

    SInt local_neighbor_row, local_neighbor_col;
    switch(direction) {
      case Right:
        local_neighbor_row = local_row;
        local_neighbor_col = 0;
        break;
      case Left:
        local_neighbor_row = local_row;
        local_neighbor_col = cols - 1;
        break;
      case Up:
        local_neighbor_row = rows - 1;
        local_neighbor_col = local_col;
        break;
      case Down:
        local_neighbor_row = 0;
        local_neighbor_col = local_col;
        break;
      default:
        break;
    }
    return offset + (local_neighbor_row * cols + local_neighbor_col);
  }

  void GenerateEdge(const SInt source, const SInt target) {
    SInt edge_seed = std::min(source, target) * total_rows_ * total_cols_ + std::max(source, target);
    SInt h = sampling::Spooky::hash(config_.seed + edge_seed);
    if (rng_.GenerateBinomial(h, 1, edge_probability_)) {
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

  inline SSInt DirectionRow(Direction direction) {
    switch(direction) {
      case Up:
        return -1;
      case Down:
        return 1;
      default:
        return 0;
    }
  }

  inline SSInt DirectionColumn(Direction direction) {
    switch(direction) {
      case Left:
        return -1;
      case Right:
        return 1;
      default:
        return 0;
    }
  }

  SInt OffsetForChunk(const SInt chunk) {
    SInt chunk_row, chunk_col;
    Decode(chunk, chunk_row, chunk_col);
      
    // Compute start vertex coordinates from chunk
    SInt vertex_row = chunk_row * rows_per_chunk_ + std::min(chunk_row, remaining_rows_);
    SInt vertex_col = chunk_col * cols_per_chunk_ + std::min(chunk_col, remaining_cols_);

    // Compute offset of start vertex
    SInt upper_rectangle = vertex_row * total_cols_;
    SInt upperleft_rectangle = vertex_col * vertex_row; 

    SInt next_vertex_row = (chunk_row + 1) * rows_per_chunk_ + std::min(chunk_row + 1, remaining_rows_);
    SInt left_rectangle = vertex_col * next_vertex_row;

    return upper_rectangle + left_rectangle - upperleft_rectangle;
  }

  inline void Decode(const SInt id, SInt &x, SInt &y) const {
    x = id / chunks_per_dim_;
    y = id % chunks_per_dim_;
    // m2D_d_sLUT(id, y, x);
  }

  // Chunk/vertex coding
  inline SInt Encode(const SInt x, const SInt y) const {
    return x * chunks_per_dim_ + y;
    // return m2D_e_sLUT<SInt>(y, x);
  }
};

}
#endif
