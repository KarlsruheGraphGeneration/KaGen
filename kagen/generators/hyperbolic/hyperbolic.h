/*******************************************************************************
 * include/generators/hyperbolic/hyperbolic.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include <google/dense_hash_map>
#include <iostream>
#include <limits>
#include <tuple>
#include <vector>

#include "hash.hpp"
#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/generator.h"
#include "kagen/tools/geometry.h"
#include "kagen/tools/mersenne.h"
#include "kagen/tools/rng_wrapper.h"
#include "kagen/tools/sorted_mersenne.h"
#include "methodD.hpp"

namespace kagen {
class Hyperbolic : public Generator {
public:
    // n, min_r, max_r, generated, offset
    using Annulus = std::tuple<SInt, LPFloat, LPFloat, bool, SInt>;
    // n, min_phi, max_phi, offset
    using Chunk = std::tuple<SInt, LPFloat, LPFloat, SInt>;
    // n, min_phi, max_phi, generated, generated
    using Cell = std::tuple<SInt, LPFloat, LPFloat, bool, SInt>;
    // phi, r, x, y, gamma, id
    using Vertex = std::tuple<LPFloat, LPFloat, LPFloat, LPFloat, LPFloat, SInt>;

    Hyperbolic(const PGeneratorConfig& config, PEID rank, PEID size);

    int Requirements() const final;

    bool AlmostUndirected() const final;

    bool InvalidVertexRangeIfEmpty() const final;

protected:
    void GenerateImpl() final;

private:
    // Config
    const PGeneratorConfig& config_;
    PEID                    rank_, size_;

    // Variates
    RNGWrapper     rng_;
    Mersenne       mersenne;
    SortedMersenne sorted_mersenne;

    // Constants and variables
    LPFloat alpha_, target_r_, cosh_target_r_, pdm_target_r_;
    LPFloat pe_min_phi_, pe_max_phi_;
    LPFloat clique_thres_;
    SInt    local_chunks_;
    SInt    local_chunk_start_, local_chunk_end_;
    SInt    total_annuli_;
    SInt    start_node_, num_nodes_;
    // Eps
    LPFloat chunk_eps_, cell_eps_, point_eps_;
    // State
    SInt    current_annulus_, current_chunk_, current_cell_;
    LPFloat current_min_phi_, current_max_phi_;
    SInt    right_processed_chunk_, right_processed_cell_;

    // Data structures
    google::dense_hash_map<SInt, Annulus>             annuli_;
    google::dense_hash_map<SInt, Chunk>               chunks_;
    google::dense_hash_map<SInt, Cell>                cells_;
    google::dense_hash_map<SInt, std::vector<Vertex>> vertices_;

    // Avoid costly recomputations
    std::vector<SInt>                        global_cell_ids_;
    std::vector<SInt>                        cells_per_annulus_;
    std::vector<std::pair<LPFloat, LPFloat>> boundaries_;

    void ComputeAnnuli(SInt chunk_id);

    void ComputeChunk(SInt chunk_id);

    void ComputeChunk(
        SInt chunk_id, SInt n, SInt k, LPFloat min_phi, LPFloat max_phi, SInt chunk_start, SInt level, SInt offset);

    void GenerateCells(SInt annulus_id, SInt chunk_id);

    void GenerateVertices(SInt annulus_id, SInt chunk_id, SInt cell_id);

    void GenerateEdges(SInt annulus_id, SInt chunk_id);

    void QueryBoth(SInt annulus_id, SInt chunk_id, SInt cell_id, const Vertex& q);

    void Query(SInt annulus_id, SInt chunk_id, SInt cell_id, const Vertex& q, bool search_down = true);

    void
    QueryRightNeighbor(SInt annulus_id, SInt chunk_id, SInt cell_id, const Vertex& q, bool phase, bool search_down);

    void QueryLeftNeighbor(SInt annulus_id, SInt chunk_id, SInt cell_id, const Vertex& q, bool phase, bool search_down);

    void GenerateGridEdges(SInt annulus_id, SInt chunk_id, SInt cell_id, const Vertex& q);

    std::pair<LPFloat, LPFloat> GetBoundaryPhis(LPFloat boundary_phi, LPFloat boundary_r, SInt annulus_id) const;

    bool OutOfBounds(LPFloat num) const;

    SInt ComputeGlobalChunkId(SInt annulus, SInt chunk) const;

    SInt ComputeGlobalCellId(SInt annulus, SInt chunk, SInt cell);

    SInt TotalGridSizeForAnnulus(SInt annulus_id);

    SInt GridSizeForAnnulus(SInt annulus_id);

    bool IsLocalChunk(SInt chunk_id) const;
};
} // namespace kagen
