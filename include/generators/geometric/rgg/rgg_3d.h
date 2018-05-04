/*******************************************************************************
 * include/generators/geometric/rgg/rgg_3d.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Daniel Funke <funke@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _RGG_3D_H_
#define _RGG_3D_H_

#include "geometric/geometric_3d.h"

template <typename EdgeCallback>
class RGG3D : public Geometric3D {
 public:
  RGG3D(const PGeneratorConfig &config, const PEID rank,
        const EdgeCallback &cb)
      : Geometric3D(config, rank), io_(config), cb_(cb) {
    // Chunk variables
    total_chunks_ = config_.k;
    chunks_per_dim_ = cbrt(config_.k);
    chunk_size_ = 1.0 / chunks_per_dim_;

    // Cell variables
    cells_per_dim_ = floor(chunk_size_ / config_.r);
    cells_per_chunk_ = cells_per_dim_ * cells_per_dim_ * cells_per_dim_;
    cell_size_ =
        config_.r + (chunk_size_ - cells_per_dim_ * config_.r) / cells_per_dim_;
    target_r_ = config_.r * config_.r;

    InitDatastructures();
  }

  void Output() const override { 
#ifdef OUTPUT_EDGES
    io_.OutputEdges(); 
#else
    io_.OutputDist(); 
#endif
  }

  SInt NumberOfEdges() const override { return io_.NumEdges(); }

 public:
  // I/O
  GeneratorIO<> io_;
  EdgeCallback cb_;
  // FILE* edge_file;

  LPFloat target_r_;

  void GenerateEdges(const SInt chunk_row, const SInt chunk_column,
                     const SInt chunk_depth) override {
    // Iterate grid cells
    for (SInt cell_row = 0; cell_row < cells_per_dim_; ++cell_row) {
      for (SInt cell_column = 0; cell_column < cells_per_dim_; ++cell_column) {
        for (SInt cell_depth = 0; cell_depth < cells_per_dim_; ++cell_depth) {
          // Iterate neighboring cells
          for (SSInt i = -1; i <= 1; i++) {
            SSInt neighbor_row = cell_row + i;
            for (SSInt j = -1; j <= 1; j++) {
              SSInt neighbor_column = cell_column + j;
              for (SSInt k = -1; k <= 1; k++) {
                SSInt neighbor_depth = cell_depth + k;

                // Compute diffs
                int horizontal_diff = 0;
                int vertical_diff = 0;
                int depth_diff = 0;
                if (neighbor_depth < 0)
                  depth_diff = -1;
                else if (neighbor_depth >= (SSInt)cells_per_dim_)
                  depth_diff = 1;
                if (neighbor_column < 0)
                  horizontal_diff = -1;
                else if (neighbor_column >= (SSInt)cells_per_dim_)
                  horizontal_diff = 1;
                if (neighbor_row < 0)
                  vertical_diff = -1;
                else if (neighbor_row >= (SSInt)cells_per_dim_)
                  vertical_diff = 1;

                // Get correct grid cells
                SInt neighbor_cell_row =
                    (neighbor_row % (SSInt)cells_per_dim_ + cells_per_dim_) %
                    cells_per_dim_;
                SInt neighbor_cell_column =
                    (neighbor_column % (SSInt)cells_per_dim_ + cells_per_dim_) %
                    cells_per_dim_;
                SInt neighbor_cell_depth =
                    (neighbor_depth % (SSInt)cells_per_dim_ + cells_per_dim_) %
                    cells_per_dim_;

                // Skip invalid cells
                if ((SSInt)chunk_row + vertical_diff < 0 ||
                    chunk_row + vertical_diff >= chunks_per_dim_ ||
                    (SSInt)chunk_column + horizontal_diff < 0 ||
                    chunk_column + horizontal_diff >= chunks_per_dim_ ||
                    (SSInt)chunk_depth + depth_diff < 0 ||
                    chunk_depth + depth_diff >= chunks_per_dim_)
                  continue;

                // Get grid buckets for each cell
                SInt chunk_id = Encode(chunk_column, chunk_row, chunk_depth);
                SInt cell_id = cell_row * cells_per_dim_ + cell_column +
                               (cells_per_dim_ * cells_per_dim_) * cell_depth;
                SInt neighbor_id =
                    Encode(chunk_column + horizontal_diff,
                           chunk_row + vertical_diff, chunk_depth + depth_diff);
                SInt neighbor_cell_id =
                    neighbor_cell_row * cells_per_dim_ + neighbor_cell_column +
                    (cells_per_dim_ * cells_per_dim_) * neighbor_cell_depth;

                // If neighbor is local chunk skip
                if (chunk_id > neighbor_id && IsLocalChunk(neighbor_id))
                  continue;
                // Skip grid cells with lower id
                if (chunk_id == neighbor_id && cell_id > neighbor_cell_id)
                  continue;

                GenerateGridEdges(chunk_id, cell_id, neighbor_id,
                                  neighbor_cell_id);
              }
            }
          }
        }
      }
    }
  }

  void GenerateGridEdges(const SInt first_chunk_id, const SInt first_cell_id,
                         const SInt second_chunk_id,
                         const SInt second_cell_id) {
    // Check if vertices not generated
    SInt first_global_cell_id = ComputeGlobalCellId(first_chunk_id, first_cell_id);
    SInt second_global_cell_id =
        ComputeGlobalCellId(second_chunk_id, second_cell_id);
    GenerateVertices(first_chunk_id, first_cell_id);
    GenerateVertices(second_chunk_id, second_cell_id);

    // Gather vertices
    const std::vector<Vertex> &vertices_first = vertices_[first_global_cell_id];
    const std::vector<Vertex> &vertices_second = vertices_[second_global_cell_id];
    // GenerateVertices(first_chunk_id, first_cell_id, first_vertex_buffer_);
    // GenerateVertices(second_chunk_id, second_cell_id, second_vertex_buffer_);
    // const std::vector<Vertex> &vertices_first = first_vertex_buffer_;
    // const std::vector<Vertex> &vertices_second = second_vertex_buffer_;
    // Generate edges
    // Same cell
    if (first_chunk_id == second_chunk_id && first_cell_id == second_cell_id) {
      for (SInt i = 0; i < vertices_first.size(); ++i) {
        const Vertex &v1 = vertices_first[i];
        for (SInt j = i + 1; ++j < vertices_second.size();) {
          const Vertex &v2 = vertices_second[j];
          // Euclidean distance
          LPFloat x = std::get<0>(v1) - std::get<0>(v2);
          LPFloat y = std::get<1>(v1) - std::get<1>(v2);
          LPFloat z = std::get<2>(v1) - std::get<2>(v2);
          if (x * x + y * y + z * z <= target_r_) {
            cb_(std::get<3>(v1), std::get<3>(v2));
#ifdef OUTPUT_EDGES
            io_.PushEdge(std::get<3>(v1), std::get<3>(v2));
#else
            io_.UpdateDist(std::get<3>(v1));
            io_.UpdateDist(std::get<3>(v2));
#endif
            // fprintf(edge_file, "e %f %f %f %f %f %f\n", std::get<0>(v1),
            // std::get<1>(v1), std::get<2>(v1), std::get<0>(v2),
            // std::get<1>(v2), std::get<2>(v2));
          }
        }
      }
    } else {
      for (SInt i = 0; i < vertices_first.size(); ++i) {
        const Vertex &v1 = vertices_first[i];
        for (SInt j = 0; j < vertices_second.size(); ++j) {
          const Vertex &v2 = vertices_second[j];
          LPFloat x = std::get<0>(v1) - std::get<0>(v2);
          LPFloat y = std::get<1>(v1) - std::get<1>(v2);
          LPFloat z = std::get<2>(v1) - std::get<2>(v2);
          if (x * x + y * y + z * z <= target_r_) {
            cb_(std::get<3>(v1), std::get<3>(v2));
#ifdef OUTPUT_EDGES
            io_.PushEdge(std::get<3>(v1), std::get<3>(v2));
#else
            io_.UpdateDist(std::get<3>(v1));
            io_.UpdateDist(std::get<3>(v2));
#endif
            // fprintf(edge_file, "e %f %f %f %f %f %f\n", std::get<0>(v1),
            // std::get<1>(v1), std::get<2>(v1), std::get<0>(v2),
            // std::get<1>(v2), std::get<2>(v2));
          }
        }
      }
    }
  }

  void GenerateCells(const SInt chunk_id) override {
    // Lazily compute chunk
    if (chunks_.find(chunk_id) == end(chunks_)) {
      ComputeChunk(chunk_id);
    }

    auto &chunk = chunks_[chunk_id];

    // Stop if cell distribution already generated
    if (std::get<4>(chunk)) return;

    SInt seed = 0;
    SInt n = std::get<0>(chunk);
    SInt offset = std::get<5>(chunk);
    LPFloat total_area = chunk_size_ * chunk_size_ * chunk_size_;
    LPFloat cell_area = cell_size_ * cell_size_ * cell_size_;

    for (SInt i = 0; i < cells_per_chunk_; ++i) {
      seed = config_.seed + chunk_id * cells_per_chunk_ + i +
             total_chunks_ * cells_per_chunk_;
      SInt h = sampling::Spooky::hash(seed);
      SInt cell_vertices =
          (SInt)rng_.GenerateBinomial(h, n, (LPFloat)cell_area / total_area);
      LPFloat cell_start_x =
          std::get<1>(chunk) +
          ((i / cells_per_dim_) % cells_per_dim_) * cell_size_;
      LPFloat cell_start_y =
          std::get<2>(chunk) + (i % cells_per_dim_) * cell_size_;
      LPFloat cell_start_z =
          std::get<3>(chunk) +
          (i / (cells_per_dim_ * cells_per_dim_)) * cell_size_;

      // Only store cells that are adjacent to local ones
      if (IsLocalChunk(chunk_id) || IsAdjacentCell(chunk_id, i))
        cells_[ComputeGlobalCellId(chunk_id, i)] = std::make_tuple(
            cell_vertices, cell_start_x, cell_start_y, cell_start_z, false, offset);

      // Update for multinomial
      n -= cell_vertices;
      offset += cell_vertices;
      total_area -= cell_area;
    }
    std::get<4>(chunk) = true;
  }

  bool IsAdjacentCell(const SInt chunk_id, const SInt cell_id) {
    SInt chunk_row, chunk_column, chunk_depth;
    Decode(chunk_id, chunk_column, chunk_row, chunk_depth);
    SInt cell_row, cell_column, cell_depth;
    DecodeCell(cell_id, cell_column, cell_row, cell_depth);
    // Iterate neighboring cells
    for (SSInt i = -1; i <= 1; i++) {
      SSInt neighbor_row = cell_row + i;
      for (SSInt j = -1; j <= 1; j++) {
        SSInt neighbor_column = cell_column + j;
        for (SSInt k = -1; k <= 1; k++) {
          SSInt neighbor_depth = cell_depth + k;

          // Compute diffs
          int horizontal_diff = 0;
          int vertical_diff = 0;
          int depth_diff = 0;
          if (neighbor_depth < 0)
            depth_diff = -1;
          else if (neighbor_depth >= (SSInt)cells_per_dim_)
            depth_diff = 1;
          if (neighbor_column < 0)
            horizontal_diff = -1;
          else if (neighbor_column >= (SSInt)cells_per_dim_)
            horizontal_diff = 1;
          if (neighbor_row < 0)
            vertical_diff = -1;
          else if (neighbor_row >= (SSInt)cells_per_dim_)
            vertical_diff = 1;

          // Skip invalid cells
          if ((SSInt)chunk_row + vertical_diff < 0 ||
              chunk_row + vertical_diff >= chunks_per_dim_ ||
              (SSInt)chunk_column + horizontal_diff < 0 ||
              chunk_column + horizontal_diff >= chunks_per_dim_ ||
              (SSInt)chunk_depth + depth_diff < 0 ||
              chunk_depth + depth_diff >= chunks_per_dim_)
            continue;

          SInt neighbor_id =
              Encode(chunk_column + horizontal_diff, chunk_row + vertical_diff,
                     chunk_depth + depth_diff);
          if (IsLocalChunk(neighbor_id)) return true;
        }
      }
    }
    return false;
  }

  inline SInt EncodeCell(const SInt x, const SInt y, const SInt z) const {
    return x + y * cells_per_dim_ + z * (cells_per_dim_ * cells_per_dim_);
  }

  inline void DecodeCell(const SInt id, SInt &x, SInt &y, SInt &z) const {
    x = id % cells_per_dim_;
    y = (id / cells_per_dim_) % cells_per_dim_;
    z = id / (cells_per_dim_ * cells_per_dim_);
  }
};

#endif
