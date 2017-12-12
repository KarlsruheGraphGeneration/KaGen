/******************************************************************************
 * gnm_undirected.h
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

#ifndef _GNM_UNDIRECTED_H_
#define _GNM_UNDIRECTED_H_

#include <iostream>
#include <vector>

#include "definitions.h"
#include "generator_config.h"
#include "generator_io.h"
#include "rng_wrapper.h"
#include "tools/spooky_hash.h"

class GNMUndirected {
 public:
  GNMUndirected(const PGeneratorConfig &config, const PEID rank)
      : config_(config), rng_(config), io_(config) { }

  void Generate() {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    SInt leftover_chunks = config_.k % size;
    SInt num_chunks = config_.k / size + ((SInt)rank < leftover_chunks);
    SInt row = rank * num_chunks +
               ((SInt)rank >= leftover_chunks ? leftover_chunks : 0);

    nodes_per_chunk_ = config_.n / config_.k;
    remaining_nodes_ = config_.n % config_.k;

    for (SInt i = 0; i < num_chunks; i++) {
      GenerateChunks(row);
      row++;
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

  // Globals
  SInt nodes_per_chunk_, remaining_nodes_;

  // Variates
  RNGWrapper<> rng_;

  // I/O
  GeneratorIO<> io_;

  void GenerateChunks(const SInt row) {
    QueryTriangular(config_.m, config_.k, config_.k, row, row, 0, 0, 1);
  }

  void QueryTriangular(const SInt m, const SInt num_rows,
                       const SInt num_columns, const SInt row_id,
                       const SInt column_id, const SInt offset_row,
                       const SInt offset_column, const SInt level) {
    // Stop if there are no edges left
    if (m <= 0) return;

    // Total number of edges;
    SInt n_row = NodesInRows(num_rows, offset_row);
    SInt n_column = NodesInColumns(num_columns, offset_column);
    HPFloat total_edges = NumTriangleEdges(n_row, n_column, config_.self_loops);

    // Base Case if only one chunk is left
    if (num_rows == 1 && num_columns == 1) {
      GenerateTriangularEdges(m, offset_row, offset_column);
      return;
    }
    // Find splitter
    SInt row_splitter = (num_rows + 1) / 2;
    SInt column_splitter = (num_columns + 1) / 2;

    // Compute nodes/edges per quadrant
    SInt ul_nodes_row = NodesInRows(row_splitter, offset_row);
    SInt ll_nodes_row = NodesInRows(num_rows / 2, offset_row + row_splitter);
    SInt ll_nodes_column = NodesInColumns(column_splitter, offset_column);
    HPFloat ul_edges = NumTriangleEdges(ul_nodes_row, ul_nodes_row);
    HPFloat ll_edges = NumRectangleEdges(ll_nodes_row, ll_nodes_column);
    HPFloat lr_edges = NumTriangleEdges(ll_nodes_row, ll_nodes_row);

    // Generate variate for quadrants
    SInt chunk_start = ChunkStart(offset_row, offset_column);
    SInt h = Spooky::Hash(config_.seed + level * config_.n + chunk_start);
    SInt upper_variate =
        (SInt)rng_.GenerateHypergeometric(h, ul_edges, m, total_edges);
    SInt ll_variate = (SInt)rng_.GenerateHypergeometric(
        h, ll_edges, m - upper_variate, ll_edges + lr_edges);

    // Recursive calls for quadrants
    // Row within upper half then column will automatically be in left half
    if (row_id < offset_row + row_splitter) {
      QueryTriangular(upper_variate, row_splitter, column_splitter, row_id,
                      column_id, offset_row, offset_column, level + 1);
      QueryColumnRectangle(ll_variate, num_rows / 2, column_splitter, row_id,
                           offset_row + row_splitter, offset_column, level + 1);
    } else {
      QueryRowRectangle(ll_variate, num_rows / 2, column_splitter, row_id,
                        offset_row + row_splitter, offset_column, level + 1);
      QueryTriangular(m - upper_variate - ll_variate, num_rows / 2,
                      num_columns / 2, row_id, column_id,
                      offset_row + row_splitter,
                      offset_column + column_splitter, level + 1);
    }
  }

  void QueryRowRectangle(const SInt m, const SInt num_rows,
                         const SInt num_columns, const SInt row_id,
                         const SInt offset_row, const SInt offset_column,
                         const SInt level) {
    // Stop if there are no edges left
    if (m <= 0) return;

    // Total number of edges;
    SInt n_row = NodesInRows(num_rows, offset_row);
    SInt n_column = NodesInColumns(num_columns, offset_column);
    HPFloat total_edges = NumRectangleEdges(n_row, n_column);

    // Base Case if only one chunk is left
    if (num_rows == 1 && num_columns == 1) {
      GenerateRectangleEdges(m, row_id, offset_column);
      return;
    }

    // Find splitters
    SInt row_splitter = (num_rows + 1) / 2;
    SInt column_splitter = (num_columns + 1) / 2;

    // Compute nodes/edges per quadrant
    SInt ul_nodes_row = NodesInRows(row_splitter, offset_row);
    SInt ul_nodes_column = NodesInColumns(column_splitter, offset_column);
    SInt ur_nodes_column =
        NodesInColumns(num_columns / 2, offset_column + column_splitter);
    HPFloat ul_edges = NumRectangleEdges(ul_nodes_row, ul_nodes_column);
    HPFloat ur_edges = NumRectangleEdges(ul_nodes_row, ur_nodes_column);

    // Generate variate for upper/lower half
    SInt chunk_start = ChunkStart(offset_row, offset_column);
    SInt h = Spooky::Hash(config_.seed + level * config_.n + chunk_start);
    SInt upper_variate = (SInt)rng_.GenerateHypergeometric(
        h, ul_edges + ur_edges, m, total_edges);

    // Recursive calls for quadrants
    // Row in upper half
    if (row_id < offset_row + row_splitter) {
      // Generate variate for upper left quadrant
      SInt ul_variate = (SInt)rng_.GenerateHypergeometric(
          h, ul_edges, upper_variate, ul_edges + ur_edges);

      QueryRowRectangle(ul_variate, row_splitter, column_splitter, row_id,
                        offset_row, offset_column, level + 1);
      QueryRowRectangle(upper_variate - ul_variate, row_splitter,
                        num_columns / 2, row_id, offset_row,
                        offset_column + column_splitter, level + 1);
    } else {  // lower half
      // Compute nodes/edges per quadrant
      SInt ll_nodes_row = NodesInRows(num_rows / 2, offset_row + row_splitter);
      HPFloat ll_edges = NumRectangleEdges(ll_nodes_row, ul_nodes_column);
      HPFloat lr_edges = NumRectangleEdges(ll_nodes_row, ur_nodes_column);

      // Generate variate for lower left quadrant
      SInt ll_variate = (SInt)rng_.GenerateHypergeometric(
          h, ll_edges, m - upper_variate, ll_edges + lr_edges);

      QueryRowRectangle(ll_variate, num_rows / 2, column_splitter, row_id,
                        offset_row + row_splitter, offset_column, level + 1);
      QueryRowRectangle(m - upper_variate - ll_variate, num_rows / 2,
                        num_columns / 2, row_id, offset_row + row_splitter,
                        offset_column + column_splitter, level + 1);
    }
  }

  void QueryColumnRectangle(const SInt m, const SInt num_rows,
                            const SInt num_columns, const SInt column_id,
                            const SInt offset_row, const SInt offset_column,
                            const SInt level) {
    // Stop if there are no edges left
    if (m <= 0) return;

    // Total number of edges;
    SInt n_row = NodesInRows(num_rows, offset_row);
    SInt n_column = NodesInColumns(num_columns, offset_column);
    HPFloat total_edges = (HPFloat)n_row * n_column;

    // Base Case if only one chunk is left
    if (num_rows == 1 && num_columns == 1) {
      if (offset_row == offset_column) return;
      GenerateRectangleEdges(m, offset_row, column_id);
      return;
    }

    // Find splitters
    SInt row_splitter = (num_rows + 1) / 2;
    SInt column_splitter = (num_columns + 1) / 2;

    // Compute nodes/edges per quadrant
    SInt ul_nodes_row = NodesInRows(row_splitter, offset_row);
    SInt ul_nodes_column = NodesInColumns(column_splitter, offset_column);
    SInt ur_nodes_column =
        NodesInColumns(num_columns / 2, offset_column + column_splitter);
    SInt ll_nodes_row = NodesInRows(num_rows / 2, offset_row + row_splitter);
    HPFloat ul_edges = NumRectangleEdges(ul_nodes_row, ul_nodes_column);
    HPFloat ur_edges = NumRectangleEdges(ul_nodes_row, ur_nodes_column);
    HPFloat ll_edges = NumRectangleEdges(ll_nodes_row, ul_nodes_column);
    HPFloat lr_edges = NumRectangleEdges(ll_nodes_row, ur_nodes_column);

    // Generate variate for upper/lower half
    SInt chunk_start = ChunkStart(offset_row, offset_column);
    SInt h = Spooky::Hash(config_.seed + level * config_.n + chunk_start);
    SInt upper_variate = (SInt)rng_.GenerateHypergeometric(
        h, ul_edges + ur_edges, m, total_edges);
    SInt ul_variate = (SInt)rng_.GenerateHypergeometric(
        h, ul_edges, upper_variate, ul_edges + ur_edges);
    SInt ll_variate = (SInt)rng_.GenerateHypergeometric(
        h, ll_edges, m - upper_variate, ll_edges + lr_edges);

    // Recursive calls for quadrants
    // Column in left half
    if (column_id < offset_column + column_splitter) {
      QueryColumnRectangle(ul_variate, row_splitter, column_splitter, column_id,
                           offset_row, offset_column, level + 1);
      QueryColumnRectangle(ll_variate, num_rows / 2, column_splitter, column_id,
                           offset_row + row_splitter, offset_column, level + 1);
      // Column in right half
    } else {
      QueryColumnRectangle(upper_variate - ul_variate, row_splitter,
                           num_columns / 2, column_id, offset_row,
                           offset_column + column_splitter, level + 1);
      QueryColumnRectangle(m - upper_variate - ll_variate, num_rows / 2,
                           num_columns / 2, column_id,
                           offset_row + row_splitter,
                           offset_column + column_splitter, level + 1);
    }
  }

  void GenerateTriangularEdges(const SInt m, const SInt row_id,
                               const SInt column_id) {
    // No self loops -> skip diagonal entries
    SInt offset_row = OffsetInRow(row_id);
    SInt offset_column = OffsetInColumn(column_id);
    if (!config_.self_loops) offset_row += 1;

    // Number edges
    SInt n_row = NodesInRow(row_id);
    SInt n_column = NodesInColumn(column_id);
    HPFloat total_edges = NumTriangleEdges(n_row, n_column, config_.self_loops);

    // Sample from [1, total_edges]
    SInt h =
        Spooky::Hash(config_.seed + (((row_id + 1) * row_id) / 2) + column_id);
    rng_.GenerateSample(h, total_edges, m, [&](SInt sample) {
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

  void GenerateRectangleEdges(const SInt m, const SInt row_id,
                              const SInt column_id) {
    SInt offset_row = OffsetInRow(row_id);
    SInt offset_column = OffsetInColumn(column_id);

    // Sample from [1, num_edges]
    SInt n_row = NodesInRow(row_id);
    SInt n_column = NodesInColumn(column_id);
    HPFloat total_edges = NumRectangleEdges(n_row, n_column);

    // Sample from [1, total_edges]
    SInt h =
        Spooky::Hash(config_.seed + (((row_id + 1) * row_id) / 2) + column_id);
    rng_.GenerateSample(h, total_edges, m, [&](SInt sample) {
      SInt i = (sample - 1) / n_column;
      SInt j = (sample - 1) % n_column;
#ifdef OUTPUT_EDGES
      io_.PushEdge(i + offset_row, j + offset_column);
#else 
      io_.UpdateDist(i + offset_row);
      io_.UpdateDist(j + offset_column);
#endif
    });
  }

  inline SInt NodesInRows(const SInt rows, const SInt offset) const {
    return nodes_per_chunk_ * rows +
           std::max(remaining_nodes_ - offset, (SInt)0);
  }

  inline SInt NodesInColumns(const SInt columns, const SInt offset) const {
    return nodes_per_chunk_ * columns +
           std::max(remaining_nodes_ - offset, (SInt)0);
  }

  inline SInt NodesInRow(const SInt row) const {
    return nodes_per_chunk_ + (row < config_.k);
  }

  inline SInt NodesInColumn(const SInt column) const {
    return nodes_per_chunk_ + (column < config_.k);
  }

  inline SInt OffsetInRow(const SInt row) const {
    return nodes_per_chunk_ * row + std::min(remaining_nodes_, row);
  }

  inline SInt OffsetInColumn(const SInt column) const {
    return nodes_per_chunk_ * column + std::min(remaining_nodes_, column);
  }

  inline SInt ChunkStart(const SInt row, const SInt column) const {
    return (((row + 1) * row) / 2) + column;
  }

  inline HPFloat NumTriangleEdges(const HPFloat row, const HPFloat column,
                                  bool loops = false) const {
    return (loops && config_.self_loops) ? row * (column + 1) / 2
                                         : row * (column - 1) / 2;
  }

  inline HPFloat NumRectangleEdges(const HPFloat row,
                                   const HPFloat column) const {
    return row * column;
  }
};

#endif
