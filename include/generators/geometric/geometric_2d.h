/*******************************************************************************
 * include/generators/geometric/geometric_2d.h
 *
 * Copyright (C) 2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Daniel Funke <funke@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _GEOMETRIC_2D_H_
#define _GEOMETRIC_2D_H_

#include <google/dense_hash_map>
#include <iostream>
#include <tuple>
#include <vector>

#include "definitions.h"
#include "generator_config.h"
#include "generator_io.h"
#include "geometry.h"
#include "libmorton/morton2D.h"
#include "rng_wrapper.h"
#include "mersenne.h"
#include "hash.hpp"

namespace kagen {

class Geometric2D {
 public:
  // n, x_off, y_off, generated, offset
  using Chunk = std::tuple<SInt, LPFloat, LPFloat, bool, SInt>;
  // n, x_off, y_off, generated, offset
  using Cell = std::tuple<SInt, LPFloat, LPFloat, bool, SInt>;
  // x, y, id
  using Vertex = std::tuple<LPFloat, LPFloat, SInt>;

  Geometric2D(PGeneratorConfig &config, const PEID /* rank */)
      : config_(config), rng_(config) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);

    // Vertex range
    start_node_ = std::numeric_limits<SInt>::max();
    num_nodes_ = 0;
  }

  void Generate() {
    // Generate per PE point distribution
    for (SInt i = local_chunk_start_; i < local_chunk_end_; ++i)
      ComputeChunk(i);

    // Generate local chunks and edges
    for (SInt i = local_chunk_start_; i < local_chunk_end_; ++i)
      GenerateChunk(i);
  }

  std::pair<SInt, SInt> GetVertexRange() {
    return std::make_pair(start_node_, start_node_ + num_nodes_ - 1);
  }

  virtual void Output() const = 0;

  virtual SInt NumberOfEdges() const = 0;

 protected:
  // Config
  PGeneratorConfig &config_;
  PEID rank_, size_;

  // Variates
  RNGWrapper rng_;
  Mersenne mersenne;

  // Constants and variables
  LPFloat chunk_size_;
  SInt total_chunks_, chunks_per_dim_, local_chunk_start_, local_chunk_end_;
  LPFloat cell_size_;
  SInt cells_per_chunk_, cells_per_dim_;
  SInt start_node_, num_nodes_;

  // Data structures
  // std::vector<Chunk> chunks_;
  google::dense_hash_map<SInt, Chunk> chunks_;
  // std::vector<Cell> cells_;
  google::dense_hash_map<SInt, Cell> cells_;
  // std::vector<std::vector<Vertex>> vertices_;
  google::dense_hash_map<SInt, std::vector<Vertex>> vertices_;

  void InitDatastructures() {
    // Chunk distribution
    SInt leftover_chunks = total_chunks_ % size_;
    SInt local_chunks = total_chunks_ / size_ + ((SInt)rank_ < leftover_chunks);

    // Compute local chunk range
    local_chunk_start_ = rank_ * local_chunks +
                         ((SInt)rank_ >= leftover_chunks ? leftover_chunks : 0);
    local_chunk_end_ = local_chunk_start_ + local_chunks;

    // Init data structures
    chunks_.set_empty_key(total_chunks_);
    cells_.set_empty_key(total_chunks_ * cells_per_chunk_);
    vertices_.set_empty_key(total_chunks_ * cells_per_chunk_);
    // chunks_.resize(total_chunks_);
    // cells_.resize(total_chunks_ * cells_per_chunk_);
    // vertices_.resize(total_chunks_ * cells_per_chunk_);
    // edge_file = fopen((config_.debug_output + std::to_string(rank_)).c_str(),
    // "w");
  }

  void ComputeChunk(const SInt chunk_id) {
    ComputeChunk(chunk_id, config_.n, chunks_per_dim_, chunks_per_dim_, 0, 0, 1, 0);
  }

  void ComputeChunk(const SInt chunk_id, const SInt n, const SInt row_k,
                    const SInt column_k, const SInt chunk_start_row,
                    const SInt chunk_start_column, const SInt level, const SInt offset) {
    // Stop if chunk exists
    if (chunks_.find(chunk_id) != end(chunks_)) return;

    // Stop if no nodes remain
    if (n <= 0 || n > config_.n) return;

    SInt chunk_row, chunk_column;
    Decode(chunk_id, chunk_column, chunk_row);

    SInt chunk_start = Encode(chunk_start_column, chunk_start_row);

    // Base case
    if (row_k == 1 && column_k == 1) {
      chunks_[chunk_start] =
          std::make_tuple(n, chunk_start_row * chunk_size_,
                          chunk_start_column * chunk_size_, false, offset);
      if (IsLocalChunk(chunk_id)) {
        if (start_node_ > offset) start_node_ = offset;
        num_nodes_ += n;
      }
      return;
    }

    // Find splitter
    SInt row_splitter = (row_k + 1) / 2;
    SInt column_splitter = (column_k + 1) / 2;

    // Generate variate for upper/lower half
    SInt h = sampling::Spooky::hash(config_.seed + chunk_start + level * total_chunks_);
    SInt v_variate = rng_.GenerateBinomial(h, n, (LPFloat)row_splitter / row_k);

    // Upper half
    if (chunk_row < row_splitter + chunk_start_row) {
      // Generate variate for left/right half
      SInt h_variate = rng_.GenerateBinomial(h, v_variate, (LPFloat)column_splitter / column_k);

      // Upper left/right quadrant
      if (chunk_column < column_splitter + chunk_start_column)
        ComputeChunk(chunk_id, h_variate, row_splitter, column_splitter,
                     chunk_start_row, chunk_start_column, level + 1, offset);
      else
        ComputeChunk(chunk_id, v_variate - h_variate, row_splitter,
                     column_k - column_splitter, chunk_start_row,
                     chunk_start_column + column_splitter, level + 1, offset + h_variate);
    } else {
      // Lower half
      // Generate variate for left/right half
      SInt h_variate = rng_.GenerateBinomial(h, n - v_variate, (LPFloat)column_splitter / column_k);

      // Lower left/right quadrant
      if (chunk_column < column_splitter + chunk_start_column)
        ComputeChunk(chunk_id, h_variate, row_k - row_splitter, column_splitter,
                     chunk_start_row + row_splitter, chunk_start_column,
                     level + 1, offset + v_variate);
      else
        ComputeChunk(chunk_id, n - v_variate - h_variate, row_k - row_splitter,
                     column_k - column_splitter, chunk_start_row + row_splitter,
                     chunk_start_column + column_splitter, level + 1, offset + v_variate + h_variate);
    }
  }

  void GenerateChunk(const SInt chunk_id) {
    SInt chunk_row, chunk_column;
    Decode(chunk_id, chunk_column, chunk_row);
    // Generate local cells and fill
    GenerateCells(chunk_id);
    for (SInt i = 0; i < cells_per_chunk_; ++i) GenerateVertices(chunk_id, i);
    // Generate edges and vertices on demand
    GenerateEdges(chunk_row, chunk_column);
  }

  virtual void GenerateCells(const SInt chunk_id) {
    // Lazily compute chunk
    if (chunks_.find(chunk_id) == end(chunks_)) {
      ComputeChunk(chunk_id);
    }

    auto &chunk = chunks_[chunk_id];

    // Stop if cell distribution already generated
    if (std::get<3>(chunk)) return;

    SInt seed = 0;
    SInt n = std::get<0>(chunk);
    SInt offset = std::get<4>(chunk);
    LPFloat total_area = chunk_size_ * chunk_size_;
    LPFloat cell_area = cell_size_ * cell_size_;

    for (SInt i = 0; i < cells_per_chunk_; ++i) {
      seed = config_.seed + chunk_id * cells_per_chunk_ + i +
             total_chunks_ * cells_per_chunk_;
      SInt h = sampling::Spooky::hash(seed);
      SInt cell_vertices = rng_.GenerateBinomial(h, n, cell_area / total_area);
      LPFloat cell_start_x =
          std::get<1>(chunk) + (i / cells_per_dim_) * cell_size_;
      LPFloat cell_start_y =
          std::get<2>(chunk) + (i % cells_per_dim_) * cell_size_;
      if (cell_vertices != 0) {
        cells_[ComputeGlobalCellId(chunk_id, i)] =
            std::make_tuple(cell_vertices, cell_start_x, cell_start_y, false, offset);
      }

      // Update for multinomial
      n -= cell_vertices;
      offset += cell_vertices;
      total_area -= cell_area;
    }
    std::get<3>(chunk) = true;
  }

  void GenerateVertices(const SInt chunk_id, const SInt cell_id) {
    // Lazily compute chunk
    if (chunks_.find(chunk_id) == end(chunks_)) ComputeChunk(chunk_id);
    auto &chunk = chunks_[chunk_id];

    // Lazily compute cell distribution
    if (!std::get<3>(chunk)) GenerateCells(chunk_id);

    // Stop if cell already generated
    SInt global_cell_id = ComputeGlobalCellId(chunk_id, cell_id);
    if (cells_.find(global_cell_id) == end(cells_)) return;
    auto &cell = cells_[global_cell_id];
    if (std::get<3>(cell)) return;

    // Compute vertex distribution
    SInt n = std::get<0>(cell);
    SInt offset = std::get<4>(cell);
    LPFloat start_x = std::get<1>(cell);
    LPFloat start_y = std::get<2>(cell);

    SInt seed = config_.seed + chunk_id * cells_per_chunk_ + cell_id;
    SInt h = sampling::Spooky::hash(seed);
    mersenne.RandomInit(h);
    std::vector<Vertex> &cell_vertices = vertices_[global_cell_id];
    cell_vertices.reserve(n);
    for (SInt i = 0; i < n; ++i) {
      // Compute coordinates
      LPFloat x = mersenne.Random() * cell_size_ + start_x;
      LPFloat y = mersenne.Random() * cell_size_ + start_y;

      cell_vertices.emplace_back(x, y, offset + i);
      // fprintf(edge_file, "v %f %f\n", x, y);
    }
    std::get<3>(cell) = true;
  }

  void GenerateVertices(const SInt chunk_id, const SInt cell_id,
                        std::vector<Vertex> &vertex_buffer) {
    // Lazily compute chunk
    if (chunks_.find(chunk_id) == end(chunks_)) ComputeChunk(chunk_id);
    auto &chunk = chunks_[chunk_id];

    // Lazily compute cell distribution
    if (!std::get<3>(chunk)) GenerateCells(chunk_id);

    // Compute vertex distribution
    SInt global_cell_id = ComputeGlobalCellId(chunk_id, cell_id);
    if (cells_.find(global_cell_id) == end(cells_)) return;
    auto &cell = cells_[global_cell_id];

    SInt n = std::get<0>(cell);
    SInt offset = std::get<4>(cell);
    LPFloat start_x = std::get<1>(cell);
    LPFloat start_y = std::get<2>(cell);

    SInt seed = config_.seed + chunk_id * cells_per_chunk_ + cell_id;
    SInt h = sampling::Spooky::hash(seed);
    mersenne.RandomInit(h);
    vertex_buffer.clear();
    vertex_buffer.reserve(n);
    for (SInt i = 0; i < n; ++i) {
      // Compute coordinates
      LPFloat x = mersenne.Random() * cell_size_ + start_x;
      LPFloat y = mersenne.Random() * cell_size_ + start_y;

      vertex_buffer.emplace_back(x, y, offset + i);
      // fprintf(edge_file, "v %f %f\n", x, y);
    }
  }

  virtual void GenerateEdges(const SInt chunk_row, const SInt chunk_column) = 0;

  inline SInt ComputeGlobalCellId(const SInt chunk_id, const SInt cell_id) const {
    return chunk_id * cells_per_chunk_ + cell_id;
  }

  inline bool IsLocalChunk(const SInt chunk_id) const {
    return (chunk_id >= local_chunk_start_ && chunk_id < local_chunk_end_);
  }

  // Chunk coding
  inline SInt Encode(const SInt x, const SInt y) const {
      return libmorton::m2D_e_sLUT<SInt>(x, y);
    // return x + y * chunks_per_dim_;
  }

  inline void Decode(const SInt id, SInt &x, SInt &y) const {
      libmorton::m2D_d_sLUT(id, x, y);
    // x = id % chunks_per_dim_;
    // y = id / chunks_per_dim_;
  }
};

}
#endif
