/*******************************************************************************
 * include/generators/hyperbolic/hyperbolic.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/tools/hash_map.h"
#include "kagen/tools/mersenne.h"
#include "kagen/tools/rng_wrapper.h"
#include "kagen/tools/sorted_mersenne.h"

#include <mpi.h>

#include <tuple>
#include <vector>

namespace kagen {
class HyperbolicFactory : public GeneratorFactory {
public:
    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, PEID rank, PEID size, bool output) const override;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const override;
};

template <typename Double>
class Hyperbolic : public virtual Generator, private EdgeListOnlyGenerator {
public:
    // n, min_r, max_r, generated, offset
    using Annulus = std::tuple<SInt, Double, Double, bool, SInt>;
    // n, min_phi, max_phi, offset
    using Chunk = std::tuple<SInt, Double, Double, SInt>;
    // n, min_phi, max_phi, generated, generated
    using Cell = std::tuple<SInt, Double, Double, bool, SInt>;
    // phi, r, x, y, gamma, id
    using Vertex = std::tuple<Double, Double, Double, Double, Double, SInt>;

    Hyperbolic(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void GenerateEdgeList() final;

    void FinalizeEdgeList(MPI_Comm comm) final;

private:
    // Config
    const PGeneratorConfig& config_;
    PEID                    rank_, size_;

    // Variates
    RNGWrapper<>   rng_;
    Mersenne       mersenne;
    SortedMersenne sorted_mersenne;

    // Constants and variables
    Double alpha_, target_r_, cosh_target_r_, pdm_target_r_;
    Double pe_min_phi_, pe_max_phi_;
    Double clique_thres_;
    SInt   local_chunks_;
    SInt   local_chunk_start_, local_chunk_end_;
    SInt   total_annuli_;
    SInt   num_nodes_;
    // Eps
    Double chunk_eps_, cell_eps_, point_eps_;
    // State
    SInt   current_annulus_, current_chunk_, current_cell_;
    Double current_min_phi_, current_max_phi_;
    SInt   right_processed_chunk_, right_processed_cell_;

    // Data structures
    HashMap<SInt, Annulus>             annuli_;
    HashMap<SInt, Chunk>               chunks_;
    HashMap<SInt, Cell>                cells_;
    HashMap<SInt, std::vector<Vertex>> vertices_;

    // Avoid costly recomputations
    std::vector<SInt>                      global_cell_ids_;
    std::vector<SInt>                      cells_per_annulus_;
    std::vector<std::pair<Double, Double>> boundaries_;

    void ComputeAnnuli(SInt chunk_id);

    void ComputeChunk(SInt chunk_id);

    void ComputeChunk(
        SInt chunk_id, SInt n, SInt k, Double min_phi, Double max_phi, SInt chunk_start, SInt level, SInt offset);

    void GenerateCells(SInt annulus_id, SInt chunk_id);

    void GenerateVertices(SInt annulus_id, SInt chunk_id, SInt cell_id);

    void GenerateEdges(SInt annulus_id, SInt chunk_id);

    void QueryBoth(SInt annulus_id, SInt chunk_id, SInt cell_id, const Vertex& q);

    void Query(SInt annulus_id, SInt chunk_id, SInt cell_id, const Vertex& q, bool search_down = true);

    bool QueryRightNeighbor(SInt annulus_id, SInt chunk_id, SInt cell_id, const Vertex& q, bool phase);

    bool QueryLeftNeighbor(SInt annulus_id, SInt chunk_id, SInt cell_id, const Vertex& q, bool phase, bool search_down);

    void GenerateGridEdges(SInt annulus_id, SInt chunk_id, SInt cell_id, const Vertex& q);

    std::pair<Double, Double> GetBoundaryPhis(Double boundary_phi, Double boundary_r, SInt annulus_id) const;

    bool OutOfBounds(Double num) const;

    SInt ComputeGlobalChunkId(SInt annulus, SInt chunk) const;

    SInt ComputeGlobalCellId(SInt annulus, SInt chunk, SInt cell);

    SInt TotalGridSizeForAnnulus(SInt annulus_id);

    SInt GridSizeForAnnulus(SInt annulus_id);

    bool IsLocalChunk(SInt chunk_id) const;
};

using LowPrecisionHyperbolic  = Hyperbolic<LPFloat>;
using HighPrecisionHyperbolic = Hyperbolic<HPFloat>;
} // namespace kagen
