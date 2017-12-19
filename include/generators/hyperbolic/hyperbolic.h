/******************************************************************************
 * hyperbolic.h
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

#ifndef _HYPERBOLIC_H_
#define _HYPERBOLIC_H_

#include <google/dense_hash_map>
#include <iostream>
#include <limits>
#include <tuple>
#include <vector>

#include "definitions.h"
#include "generator_config.h"
#include "generator_io.h"
#include "geometry.h"
#include "rng_wrapper.h"
#include "sorted_mersenne.h"
#include "hash.hpp"
#include "methodD.hpp"

class Hyperbolic {
 public:
  // n, min_r, max_r, offset
  using Annulus = std::tuple<SInt, LPFloat, LPFloat, SInt>;
  // n, min_phi, max_phi, generated, offset
  using Chunk = std::tuple<SInt, LPFloat, LPFloat, bool, SInt>;
  // n, min_phi, max_phi, generated, generated
  using Cell = std::tuple<SInt, LPFloat, LPFloat, bool, SInt>;
  // phi, r, x, y, gamma, id
  using Vertex = std::tuple<LPFloat, LPFloat, LPFloat, LPFloat, LPFloat, SInt>;

  Hyperbolic(const PGeneratorConfig &config, const PEID rank)
      : config_(config), rank_(rank), rng_(config), io_(config) {
    MPI_Comm_size(MPI_COMM_WORLD, &size_);

    // Globals
    alpha_ = (config_.plexp - 1) / 2;
    // alpha_ = 1;
    target_r_ = PGGeometry::GetTargetRadius(
        config_.n, config_.n * config_.avg_degree / 2, alpha_);
    cosh_target_r_ = cosh(target_r_);
    pdm_target_r_ = (cosh_target_r_ - 1) / 2;
    clique_thres_ = target_r_ / 2.0;

    // PE-specific
    total_annuli_ = floor(alpha_ * target_r_ / log(2));
    SInt chunks_per_pe = config_.k / size_;
    SInt leftover_chunks = config_.k % size_;
    local_chunks_ = chunks_per_pe + ((SInt)rank_ < leftover_chunks);

    // Compute local chunk range
    local_chunk_start_ = local_chunks_ * rank_ +
                         ((SInt)rank_ >= leftover_chunks ? leftover_chunks : 0);
    local_chunk_end_ = local_chunk_start_ + local_chunks_;

    LPFloat phi_per_chunk = 2 * M_PI / config_.k;
    pe_min_phi_ = local_chunk_start_ * phi_per_chunk;
    pe_max_phi_ = local_chunk_end_ * phi_per_chunk;

    // Init data structures
    annuli_.reserve(total_annuli_);
    chunks_.set_empty_key(total_annuli_ * config_.k);

    // Compute number of cells_
    SInt total_cells = 0;
    cells_per_annulus_.resize(total_annuli_, std::numeric_limits<SInt>::max());
    for (SInt i = 0; i < total_annuli_; ++i) {
      global_cell_ids_.emplace_back(total_cells);
      total_cells += GridSizeForAnnulus(i) * size_;
    }
    cells_.set_empty_key(total_cells + 1);
    vertices_.set_empty_key(total_cells + 1);
    boundaries_.reserve(total_annuli_);

    // Epsilon comparison
    chunk_eps_ = phi_per_chunk / 1000;
    cell_eps_ = (2 * M_PI / GridSizeForAnnulus(total_annuli_ - 1)) / 1000;
    point_eps_ = std::numeric_limits<LPFloat>::epsilon();

    // I/O (Debug)
    // edge_file = fopen((config_.debug_output + std::to_string(rank_)).c_str(), "w"); fprintf(edge_file, "target_r_ %f\n", target_r_);
  }

  void Generate() {
    // Annuli and chunks
    ComputeAnnuli(config_.n);

    // Compute local chunks
    for (SInt i = 0; i < total_annuli_; ++i) {
      for (SInt j = local_chunk_start_; j < local_chunk_end_; ++j)
        ComputeChunk(i, j);
    }

    // Local points
    for (SInt i = 0; i < total_annuli_; ++i) {
      for (SInt j = local_chunk_start_; j < local_chunk_end_; ++j) {
        GenerateCells(i, j);
        for (SInt k = 0; k < GridSizeForAnnulus(i); ++k)
          GenerateVertices(i, j, k);
      }
    }

    // Local edges
    for (SInt i = 0; i < total_annuli_; ++i) {
      for (SInt j = local_chunk_start_; j < local_chunk_end_; ++j) {
        GenerateEdges(i, j);
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
  PEID rank_, size_;

  // Variates
  RNGWrapper<VarGen<LPFloat>, VarGen<>> rng_;
  Mersenne mersenne;
  SortedMersenne sorted_mersenne;

  // I/O
  GeneratorIO<> io_;
  // FILE* edge_file;

  // Constants and variables
  LPFloat alpha_, target_r_, cosh_target_r_, pdm_target_r_;
  LPFloat pe_min_phi_, pe_max_phi_;
  LPFloat clique_thres_;
  SInt local_chunks_;
  SInt local_chunk_start_, local_chunk_end_;
  SInt total_annuli_;
  // Eps
  LPFloat chunk_eps_, cell_eps_, point_eps_;
  // State
  SInt current_annulus_, current_chunk_, current_cell_;
  LPFloat current_min_phi_, current_max_phi_;

  // Data structures
  std::vector<Annulus> annuli_;
  google::dense_hash_map<SInt, Chunk> chunks_;
  google::dense_hash_map<SInt, Cell> cells_;
  google::dense_hash_map<SInt, std::vector<Vertex>> vertices_;

  // Avoid costly recomputations
  std::vector<SInt> global_cell_ids_;
  std::vector<SInt> cells_per_annulus_;
  std::vector<std::pair<LPFloat, LPFloat>> boundaries_;

  void ComputeAnnuli(SInt n) {
    LPFloat min_r = 0;
    LPFloat total_area = PGGeometry::RadiusToHyperbolicArea(alpha_ * target_r_);
    SInt level = 1;
    SInt offset = 0;

    for (SInt i = 1; i < total_annuli_ + 1; i++) {
      // Distribute points
      LPFloat max_r = i * target_r_ / total_annuli_;
      LPFloat ring_area = PGGeometry::RadiusToHyperbolicArea(alpha_ * max_r) -
                          PGGeometry::RadiusToHyperbolicArea(alpha_ * min_r);

      // Variate
      SInt h = sampling::Spooky::hash(config_.seed + level * total_annuli_ + i);
      SInt n_annulus =
          (SInt)rng_.GenerateBinomial(h, n, (LPFloat)ring_area / total_area);

      // Push annuli_
      annuli_.emplace_back(n_annulus, min_r, max_r, offset);
      boundaries_.emplace_back(cosh(min_r), sinh(min_r));
      // fprintf(edge_file, "a %llu %f %f\n", n_annulus, min_r, max_r);
      min_r = max_r;
      n -= n_annulus;
      offset += n_annulus;
      total_area -= ring_area;
    }
    if (config_.thres > 0) 
      clique_thres_ = std::get<2>(annuli_[config_.thres]);
  }

  void ComputeChunk(const SInt annulus_id, const SInt chunk_id) {
    auto &annulus = annuli_[annulus_id];
    ComputeChunk(annulus_id, chunk_id, std::get<0>(annulus), config_.k, 0,
                 2 * M_PI, 0, 1, std::get<3>(annulus));
  }

  void ComputeChunk(const SInt annulus_id, const SInt chunk_id, const SInt n,
                    const SInt k, const LPFloat min_phi, const LPFloat max_phi,
                    const SInt chunk_start, const SInt level, const SInt offset) {
    // Base case
    if (k == 1) {
      chunks_[ComputeGlobalChunkId(annulus_id, chunk_start)] =
          std::make_tuple(n, min_phi, max_phi, false, offset);
      // fprintf(edge_file, "c %llu %f %f %f %f\n", n, min_phi, max_phi, std::get<1>(annuli_[annulus_id]), std::get<2>(annuli_[annulus_id]));
      return;
    }

    // Compute point distribution
    SInt midk = (k + 1) / 2;

    // Generate variate
    SInt h = sampling::Spooky::hash(config_.seed + level * config_.k + chunk_start + annulus_id);
    SInt splitter_variate = rng_.GenerateBinomial(h, n, (LPFloat)midk / k);

    // Compute splitter
    LPFloat middlePhi = (max_phi - min_phi) * ((LPFloat)midk / k) + min_phi;

    // Recurse
    if (chunk_id < chunk_start + midk)
      ComputeChunk(annulus_id, chunk_id, splitter_variate, midk, min_phi,
                   middlePhi, chunk_start, level + 1, offset);
    else
      ComputeChunk(annulus_id, chunk_id, n - splitter_variate, k - midk,
                   middlePhi, max_phi, chunk_start + midk, level + 1, offset + splitter_variate);
  }

  void GenerateCells(const SInt annulus_id, SInt chunk_id) {
    bool clique = false;
    auto &annulus = annuli_[annulus_id];
    if (std::get<1>(annulus) < clique_thres_) {
      chunk_id = local_chunk_start_;
      clique = true;
    }

    // Lazily compute chunk
    SInt global_chunk_id = ComputeGlobalChunkId(annulus_id, chunk_id);
    if (chunks_.find(global_chunk_id) == end(chunks_))
      ComputeChunk(annulus_id, chunk_id);
    auto &chunk = chunks_[global_chunk_id];

    // Stop if cell distribution already generated
    if (std::get<3>(chunk)) return;

    SInt n, offset, seed = 0;
    LPFloat min_phi, max_phi;

    // Retrieve parameters
    if (clique) {
      n = std::get<0>(annulus);
      offset = std::get<3>(annulus);
      min_phi = 0.0;
      max_phi = 2 * M_PI;
      seed = config_.seed + annulus_id + config_.n;
    } else {
      n = std::get<0>(chunk);
      offset = std::get<4>(chunk);
      min_phi = std::get<1>(chunk);
      max_phi = std::get<2>(chunk);
    }

    LPFloat total_phi = max_phi - min_phi;
    LPFloat grid_phi = total_phi / GridSizeForAnnulus(annulus_id);
    for (SInt i = 0; i < GridSizeForAnnulus(annulus_id); ++i) {
      // Variate
      if (!clique)
        seed = config_.seed + annulus_id * config_.k + chunk_id + i + config_.n;
      SInt h = sampling::Spooky::hash(seed);
      SInt n_cell = (SInt)rng_.GenerateBinomial(h, n, (LPFloat)grid_phi / total_phi);

      SInt global_cell_id = ComputeGlobalCellId(annulus_id, chunk_id, i);
      cells_[global_cell_id] = std::make_tuple(n_cell, min_phi + (grid_phi * i), min_phi + (grid_phi * (i + 1)), false, offset);
      // fprintf(edge_file, "g %llu %f %f %f %f\n", n_cell, min_phi + (grid_phi * i), min_phi + (grid_phi * (i+1)), std::get<1>(annuli_[annulus_id]), std::get<2>(annuli_[annulus_id]));
      n -= n_cell;
      offset += n_cell;
      total_phi -= grid_phi;
    }
    std::get<3>(chunk) = true;
  }

  void GenerateVertices(const SInt annulus_id, SInt chunk_id,
                        const SInt cell_id) {
    bool clique = false;
    auto &annulus = annuli_[annulus_id];
    if (std::get<1>(annulus) < clique_thres_) {
      chunk_id = local_chunk_start_;
      clique = true;
    }

    // Lazily compute chunk
    SInt global_chunk_id = ComputeGlobalChunkId(annulus_id, chunk_id);
    if (chunks_.find(global_chunk_id) == end(chunks_))
      ComputeChunk(annulus_id, chunk_id);
    auto &chunk = chunks_[global_chunk_id];

    // Lazily compute cells distribution
    if (!std::get<3>(chunk)) GenerateCells(annulus_id, chunk_id);

    // Check if cell was generated
    SInt global_cell_id = ComputeGlobalCellId(annulus_id, chunk_id, cell_id);
    auto &cell = cells_[global_cell_id];
    if (std::get<3>(cell)) return;

    // Compute vertex distribution
    SInt n = std::get<0>(cell);
    SInt offset = std::get<4>(cell);
    LPFloat min_phi = std::get<1>(cell);
    LPFloat max_phi = std::get<2>(cell);
    LPFloat min_r = std::get<1>(annulus);
    LPFloat max_r = std::get<2>(annulus);

    SInt seed = 0;
    if (clique)
      seed = config_.seed +
             annulus_id * config_.k * GridSizeForAnnulus(annulus_id);
    else
      seed = config_.seed +
             annulus_id * config_.k * GridSizeForAnnulus(annulus_id) +
             chunk_id * GridSizeForAnnulus(annulus_id) + cell_id + config_.n;

    SInt h = sampling::Spooky::hash(seed);
    mersenne.RandomInit(h);
    sorted_mersenne.RandomInit(h, n);
    LPFloat mincdf = cosh(alpha_ * min_r);
    LPFloat maxcdf = cosh(alpha_ * max_r);
    std::vector<Vertex> &cell_vertices = vertices_[global_cell_id];
    cell_vertices.reserve(n);
    for (SInt i = 0; i < n; i++) {
      // Compute coordinates
      LPFloat angle = sorted_mersenne.Random() * (max_phi - min_phi) + min_phi;
      LPFloat radius =
          acosh(mersenne.Random() * (maxcdf - mincdf) + mincdf) / alpha_;
      // fprintf(edge_file, "p %f %f %d\n", radius, angle, rank_);

      // Perform pdm transformation
      LPFloat inv_len = (cosh(radius) + 1.0) / 2.0;
      LPFloat pdm_radius = sqrt(1.0 - 1.0 / inv_len);
      LPFloat x = pdm_radius * sin(angle);
      LPFloat y = pdm_radius * cos(angle);
      LPFloat gamma = 1.0 / (1.0 - pdm_radius * pdm_radius);
      cell_vertices.emplace_back(angle, radius, x, y, gamma, offset + i);
    }
    std::get<3>(cell) = true;
  }

  void GenerateEdges(const SInt annulus_id, const SInt chunk_id) {
    current_annulus_ = annulus_id;
    current_chunk_ = chunk_id;
    for (SInt cell_id = 0; cell_id < GridSizeForAnnulus(annulus_id);
         ++cell_id) {
      SInt global_cell_id = ComputeGlobalCellId(annulus_id, chunk_id, cell_id);
      if (std::get<0>(cells_[global_cell_id]) == 0) continue;
      current_cell_ = cell_id;
      for (SInt i = 0; i < vertices_[global_cell_id].size(); ++i) {
        // const Vertex &v = cell_vertices[i];
        // Need copy because of hash movement
        const Vertex v = vertices_[global_cell_id][i];
        if (pe_min_phi_ > std::get<0>(v) || pe_max_phi_ < std::get<0>(v))
          continue;
        // fprintf(edge_file, "qp %f %f %f\n", std::get<1>(v), std::get<0>(v), target_r_);
        QueryBoth(annulus_id, chunk_id, cell_id, v);
      }
    }
  }

  void QueryBoth(const SInt annulus_id, const SInt chunk_id, const SInt cell_id,
                 const Vertex &q) {
    Query(annulus_id, chunk_id, cell_id, q);
    if (config_.query_both && annulus_id > 0) {
      auto &chunk = chunks_[ComputeGlobalChunkId(annulus_id - 1, chunk_id)];
      LPFloat min_chunk_phi = std::get<1>(chunk);
      LPFloat max_chunk_phi = std::get<2>(chunk);
      LPFloat grid_phi =
          (max_chunk_phi - min_chunk_phi) / GridSizeForAnnulus(annulus_id - 1);
      SInt next_cell_id = floor((std::get<0>(q) - min_chunk_phi) / grid_phi);
      Query(annulus_id - 1, chunk_id, next_cell_id, q, false);
    }
  }

  void Query(const SInt annulus_id, const SInt chunk_id, const SInt cell_id,
             const Vertex &q, bool search_down = true) {
    // Boundaries
    auto &annulus = annuli_[annulus_id];
    auto current_bounds = GetBoundaryPhis(std::get<0>(q), std::get<1>(q), annulus_id);
    current_min_phi_ = std::get<0>(current_bounds);
    current_max_phi_ = std::get<1>(current_bounds);
    
    LPFloat min_cell_phi =
        std::get<1>(cells_[ComputeGlobalCellId(annulus_id, chunk_id, cell_id)]);
    LPFloat max_cell_phi =
        std::get<2>(cells_[ComputeGlobalCellId(annulus_id, chunk_id, cell_id)]);

    // Iterate over cell
    GenerateGridEdges(annulus_id, chunk_id, cell_id, q);

    if (std::get<1>(annulus) >= clique_thres_ && std::max(TotalGridSizeForAnnulus(annulus_id), config_.k) > 1) {
      // Continue right
      if (current_min_phi_ < min_cell_phi || OutOfBounds(current_min_phi_)) {
        SInt next_chunk_id = chunk_id;
        if (cell_id == 0)
          next_chunk_id = (chunk_id + config_.k - 1) % config_.k;
        SInt next_cell_id = (cell_id + GridSizeForAnnulus(annulus_id) - 1) %
                            GridSizeForAnnulus(annulus_id);
        GenerateVertices(annulus_id, next_chunk_id, next_cell_id);
        QueryRightNeighbor(annulus_id, next_chunk_id, next_cell_id, q, 
                           std::abs(min_cell_phi - 0.0) < cell_eps_);
      }

      // Continue left
      if (current_max_phi_ > max_cell_phi || OutOfBounds(current_max_phi_)) {
        SInt next_chunk_id = chunk_id;
        if (cell_id == GridSizeForAnnulus(annulus_id) - 1)
          next_chunk_id = (chunk_id + config_.k + 1) % config_.k;
        SInt next_cell_id = (cell_id + GridSizeForAnnulus(annulus_id) + 1) %
                            GridSizeForAnnulus(annulus_id);
        GenerateVertices(annulus_id, next_chunk_id, next_cell_id);
        QueryLeftNeighbor(annulus_id, next_chunk_id, next_cell_id, q, 
                          std::abs(max_cell_phi - 2 * M_PI) < cell_eps_);
      }
    }

    // Continue with next annuli_
    SInt next_annulus;
    if (search_down)
      next_annulus = annulus_id + 1;
    else
      next_annulus = annulus_id - 1;

    if (next_annulus >= total_annuli_ || (LONG)next_annulus < 0) return;

    // Find next cell
    auto &chunk = chunks_[ComputeGlobalChunkId(next_annulus, chunk_id)];
    LPFloat min_chunk_phi = std::get<1>(chunk);
    LPFloat max_chunk_phi = std::get<2>(chunk);
    LPFloat grid_phi =
        (max_chunk_phi - min_chunk_phi) / GridSizeForAnnulus(next_annulus);
    SInt next_cell_id = floor((std::get<0>(q) - min_chunk_phi) / grid_phi);

    Query(next_annulus, chunk_id, next_cell_id, q, search_down);
  }

  void QueryRightNeighbor(const SInt annulus_id, SInt chunk_id, SInt cell_id,
                          const Vertex &q, bool phase) {
    while (true) {
      // Boundaries
      if (phase && current_min_phi_ < 0.0) current_min_phi_ += 2 * M_PI;
      if (phase && (OutOfBounds(current_min_phi_) ||
                    std::get<1>(annuli_[annulus_id]) < clique_thres_))
        return;

      auto &cell = cells_[ComputeGlobalCellId(annulus_id, chunk_id, cell_id)];
      LPFloat min_cell_phi = std::get<1>(cell);

      // Iterate over cell
      GenerateGridEdges(annulus_id, chunk_id, cell_id, q);

      phase = phase || std::abs(min_cell_phi - 0.0) < cell_eps_;
      if (current_min_phi_ < min_cell_phi || OutOfBounds(current_min_phi_)) {
        SInt next_chunk_id = chunk_id;
        if (cell_id == 0)
          next_chunk_id = (chunk_id + config_.k - 1) % config_.k;
        SInt next_cell_id = (cell_id + GridSizeForAnnulus(annulus_id) - 1) %
                            GridSizeForAnnulus(annulus_id);
        GenerateVertices(annulus_id, next_chunk_id, next_cell_id);
        cell_id = next_cell_id;
        chunk_id = next_chunk_id;
        continue;
      }
      return;
    }
  }

  void QueryLeftNeighbor(const SInt annulus_id, SInt chunk_id, SInt cell_id,
                         const Vertex &q, bool phase) {
    while (true) {
      // Boundaries
      if (phase && current_max_phi_ >= 2 * M_PI) current_max_phi_ -= 2 * M_PI;
      if (phase && (OutOfBounds(current_max_phi_) ||
                    std::get<1>(annuli_[annulus_id]) < clique_thres_))
        return;

      auto &cell = cells_[ComputeGlobalCellId(annulus_id, chunk_id, cell_id)];
      LPFloat max_cell_phi = std::get<2>(cell);

      // Iterate over cell
      GenerateGridEdges(annulus_id, chunk_id, cell_id, q);

      phase = phase || std::abs(max_cell_phi - 2 * M_PI) < cell_eps_;
      if (current_max_phi_ > max_cell_phi || OutOfBounds(current_max_phi_)) {
        SInt next_chunk_id = chunk_id;
        if (cell_id == GridSizeForAnnulus(annulus_id) - 1)
          next_chunk_id = (chunk_id + config_.k + 1) % config_.k;
        SInt next_cell_id = (cell_id + GridSizeForAnnulus(annulus_id) + 1) %
                            GridSizeForAnnulus(annulus_id);
        GenerateVertices(annulus_id, next_chunk_id, next_cell_id);
        cell_id = next_cell_id;
        chunk_id = next_chunk_id;
        continue;
      }
      return;
    }
  }

  void GenerateGridEdges(const SInt annulus_id, const SInt chunk_id,
                         const SInt cell_id, const Vertex &q) {
    // Check if vertices not generated
    SInt global_cell_id = ComputeGlobalCellId(annulus_id, chunk_id, cell_id);
    GenerateVertices(annulus_id, chunk_id, cell_id);

    // Gather vertices
    const std::vector<Vertex> &cell_vertices = vertices_[global_cell_id];
    // Same cell
    // if (current_annulus_ == annulus_id && current_chunk_ == chunk_id && current_cell_ == cell_id) {
    //   for (SInt j = 0; j < cell_vertices.size(); ++j) {
    //     const Vertex &v = cell_vertices[j];
    //     // Skip if larger angle or same angle and larger radius
    //     // TODO: Use in sequential? (slower but no duplicates)
    //     if (std::get<0>(v) > std::get<0>(q) || (std::abs(std::get<0>(v) - std::get<0>(q)) < point_eps_ && std::get<1>(v) < std::get<1>(q))) continue; 
    //     if (std::abs(std::get<1>(v) - std::get<1>(q)) < point_eps_ && std::abs(std::get<0>(v) - std::get<0>(q)) < point_eps_) continue;
    //     // Generate edge
    //     if (PGGeometry::HyperbolicDistance(q, v) <= pdm_target_r_) {
    //       // fprintf(edge_file, "e %f %f %f %f %d\n", std::get<1>(q), std::get<0>(q), std::get<1>(v), std::get<0>(v), rank_);
    //       // io_.PushEdge(std::get<0>(q), std::get<1>(q), std::get<0>(v), std::get<1>(v));
    //       io_.UpdateDist(std::get<5>(q));
    //       io_.UpdateDist(std::get<5>(v));
    //     }
    //   }
    // }
    // Different cells
    // else {
    for (SInt j = 0; j < cell_vertices.size(); ++j) {
      const Vertex &v = cell_vertices[j];
      if (PGGeometry::HyperbolicDistance(q, v) <= pdm_target_r_) {
        // fprintf(edge_file, "e %f %f %f %f %d\n", std::get<1>(q), std::get<0>(q), std::get<1>(v), std::get<0>(v), rank_);
#ifdef OUTPUT_EDGES
        io_.PushEdge(std::get<5>(q), std::get<5>(v));
#else
        io_.UpdateDist(std::get<5>(q));
        io_.UpdateDist(std::get<5>(v));
#endif
      }
    }
    // }
  }

  inline std::pair<LPFloat,LPFloat> GetBoundaryPhis(const LPFloat boundary_phi,
                                                    const LPFloat boundary_r,
                                                    const SInt annulus_id) const {
    auto &boundary = boundaries_[annulus_id];
    LPFloat cosh_min_r = std::get<0>(boundary);
    LPFloat sinh_min_r = std::get<1>(boundary);
    LPFloat diff = acos((cosh(boundary_r) * cosh_min_r - cosh_target_r_) / (sinh(boundary_r) * sinh_min_r));
    LPFloat lower_bound = boundary_phi - diff;
    LPFloat upper_bound = boundary_phi + diff;
    LPFloat min_phi = std::min(lower_bound, upper_bound);
    LPFloat max_phi = std::max(lower_bound, upper_bound);
    return std::make_pair(min_phi, max_phi);
  }

  inline bool OutOfBounds(const LPFloat num) const {
    return (num < -2*M_PI || num > 2*M_PI);
  }

  inline SInt ComputeGlobalChunkId(const SInt annulus, const SInt chunk) const {
    return annulus * config_.k + chunk;
  }

  inline SInt ComputeGlobalCellId(const SInt annulus, const SInt chunk,
                                const SInt cell) {
    return global_cell_ids_[annulus] + chunk * GridSizeForAnnulus(annulus) + cell;
  }

  SInt TotalGridSizeForAnnulus(const SInt annulus_id) {
    if (cells_per_annulus_[annulus_id] != std::numeric_limits<SInt>::max())
      return cells_per_annulus_[annulus_id];
    LPFloat min_r = annulus_id * target_r_ / total_annuli_;
    LPFloat max_r = (annulus_id + 1) * target_r_ / total_annuli_;
    LPFloat ring_area = PGGeometry::RadiusToHyperbolicArea(alpha_ * max_r) -
                        PGGeometry::RadiusToHyperbolicArea(alpha_ * min_r);
    LPFloat total_area = PGGeometry::RadiusToHyperbolicArea(alpha_ * target_r_);

    SInt exp_points = config_.n * (LPFloat)ring_area / total_area;
    SInt cells = exp_points / config_.hyp_base;

    SInt result = std::max<SInt>(1, cells);
    cells_per_annulus_[annulus_id] = result;
    return result;
  }

  inline SInt GridSizeForAnnulus(const SInt annulus_id) {
    return std::max<SInt>(1, TotalGridSizeForAnnulus(annulus_id) / size_);
  }
};

#endif
