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
#include "kagen/kagen.h"
#include "kagen/tools/hash_map.h"
#include "kagen/tools/mersenne.h"
#include "kagen/tools/rng_wrapper.h"
#include "kagen/tools/sorted_mersenne.h"

#include <mpi.h>

#include <tuple>
#include <vector>

namespace kagen {
class Hyper_HyperbolicFactory : public GeneratorFactory {
public:
    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, PEID rank, PEID size, bool output) const override;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const override;
};

template <typename Double>
struct HyperbolicHyperedgeCenter {
    Double phi;
    Double r;
    SInt   sampled_id;
    SInt   annulus_id;
};

template <typename Double>
class HyperbolicGeometryPolicy;

template <typename Double>
class Hyper_Hyperbolic : public virtual Generator, private CSROnlyGenerator {
    friend class HyperbolicGeometryPolicy<Double>;

public:
    // n, min_r, max_r, generated, offset
    using Annulus = std::tuple<SInt, Double, Double, bool, SInt>;
    // n, min_phi, max_phi, offset
    using Chunk = std::tuple<SInt, Double, Double, SInt>;
    // n, min_phi, max_phi, generated, offset, min_x, max_x, min_y, max_y
    using Cell = std::tuple<SInt, Double, Double, bool, SInt, Double, Double, Double, Double>;
    // phi, r, x, y, gamma, id, cosh_r, sinh_r, cos_phi, sin_phi
    using Vertex = std::tuple<Double, Double, Double, Double, Double, SInt, Double, Double, Double, Double>;
    struct VertexBlock {
        std::vector<Double> r;
        std::vector<SInt>   id;
        std::vector<Double> cosh_r;
        std::vector<Double> sinh_r;
        std::vector<Double> cos_phi;
        std::vector<Double> sin_phi;
        std::vector<Double> phi;

        void reserve(const SInt n) {
            r.reserve(n);
            id.reserve(n);
            cosh_r.reserve(n);
            sinh_r.reserve(n);
            cos_phi.reserve(n);
            sin_phi.reserve(n);
            phi.reserve(n);
        }

        std::size_t size() const {
            return id.size();
        }
    };

    // Used for placement of hyperball centers
    //  m, min_phi, max_phi, offset
    using CenterChunk = std::tuple<SInt, Double, Double, SInt>;

    // m, min_r, max_r, generated, offset
    using CenterAnnulus = std::tuple<SInt, Double, Double, bool, SInt>;

    // m, min_phi, max_phi, generated, offset
    using CenterCell = std::tuple<SInt, Double, Double, bool, SInt>;

    Hyper_Hyperbolic(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void GenerateCSR() final;

    void FinalizeCSR(MPI_Comm comm) final;

private:
    // Config
    PGeneratorConfig config_;
    PEID             rank_, size_;

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
    SInt                  current_annulus_, current_chunk_, current_cell_;
    Double                current_min_phi_, current_max_phi_;
    SInt                  right_processed_chunk_, right_processed_cell_;
    std::vector<SInt>     current_hyperedge_pins_;
    std::vector<PinRange> current_hyperedge_ranges_;
    Double                current_hyperedge_radius_;
    Double                current_hyperedge_pdm_radius_;

    // Data structures
    HashMap<SInt, Annulus>     annuli_;
    HashMap<SInt, Chunk>       chunks_;
    HashMap<SInt, Cell>        cells_;
    HashMap<SInt, VertexBlock> vertices_;
    // used for hyperedge centers
    HashMap<SInt, CenterChunk>   center_chunks_;
    HashMap<SInt, CenterAnnulus> center_annuli_;
    HashMap<SInt, CenterCell>    center_cells_;

    // Avoid costly recomputations
    std::vector<SInt>                      global_cell_ids_;
    std::vector<SInt>                      cells_per_annulus_;
    std::vector<std::pair<Double, Double>> boundaries_;
    std::vector<Double>                    minimum_radii_by_center_annulus_;
    std::vector<Double>                    annulus_min_r_;
    std::vector<Double>                    annulus_max_r_;
    std::vector<Double>                    annulus_min_cosh_;
    std::vector<Double>                    annulus_min_sinh_;
    std::vector<Double>                    annulus_max_cosh_;
    std::vector<Double>                    annulus_max_sinh_;

    struct DebugCenter {
        SInt   sampled_id;
        Double phi;
        Double r;
        Double radius;
    };

    struct CenterSamplingRegion {
        Double min_phi;
        Double max_phi;
        Double min_cdf;
        Double max_cdf;
        SInt   offset;
        SInt   count;
    };

    struct SampledVertex {
        Double r;
        Double phi;
        Double cosh_r;
        Double sinh_r;
        Double cos_phi;
        Double sin_phi;
        Double x;
        Double y;
    };

    std::vector<DebugCenter> debug_centers_;

    void ComputeAnnuli(SInt chunk_id);

    void ComputeChunk(SInt chunk_id);

    void GenerateCells(SInt annulus_id, SInt chunk_id);

    void GenerateVertices(SInt annulus_id, SInt chunk_id, SInt cell_id);

    void ComputeCenterAnnuli(SInt chunk_id);

    void ComputeCenterChunk(SInt chunk_id);

    void GenerateCenterCells(SInt annulus_id, SInt chunk_id);

    void GenerateHyperedges(SInt annulus_id, SInt chunk_id);

    template <typename ChunkMap, typename AnnulusMap>
    void ComputeAnnuliInto(const ChunkMap& chunks, AnnulusMap& annuli, SInt chunk_id, SInt seed_offset);

    template <typename ChunkMap>
    void ComputeChunkInto(
        ChunkMap& chunks, SInt chunk_id, SInt total_objects, SInt num_chunks, Double min_phi, Double max_phi,
        SInt chunk_start, SInt level, SInt offset, SInt seed_offset);

    template <typename ChunkMap, typename AnnulusMap, typename CellMap>
    void GenerateCellsInto(
        SInt annulus_id, SInt chunk_id, ChunkMap& chunks, AnnulusMap& annuli, CellMap& cells, SInt seed_offset);

    bool OutOfBounds(Double num) const;

    SInt ComputeGlobalChunkId(SInt annulus, SInt chunk) const;

    SInt ComputeGlobalCellId(SInt annulus, SInt chunk, SInt cell);

    SInt TotalGridSizeForAnnulus(SInt annulus_id);

    SInt GridSizeForAnnulus(SInt annulus_id);

    bool IsLocalChunk(SInt chunk_id) const;

    void BeginHyperedge(const HyperbolicHyperedgeCenter<Double>& center, Mersenne& mersenne);

    void EndHyperedge();

    void PushHyperedgePin(SInt pin);

    void PushHyperedgeRange(SInt begin, SInt end);

    void PushHyperedgeCompressed(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges);

    Double MinimumRadius(const HyperbolicHyperedgeCenter<Double>& center) const;

    Double MaximumRadius(const HyperbolicHyperedgeCenter<Double>& center) const;

    Double SampleRadius(const HyperbolicHyperedgeCenter<Double>& center);

    void PrecomputeMinimumRadii();

    Double Radius(const HyperbolicHyperedgeCenter<Double>& center, Mersenne& mersenne);

    Double FindRadiusForExpectedPins(const HyperbolicHyperedgeCenter<Double>& center, Double desired_pins);

    HyperbolicHyperedgeCenter<Double>
    SampleCenter(SInt annulus_id, SInt sampled_center_id, const CenterSamplingRegion& region);

    void SeedHyperedgeRNG(SInt sampled_center_id);

    Hyper_Hyperbolic<Double>::CenterSamplingRegion
    BuildCenterSamplingRegion(const CenterAnnulus& center_annulus, const CenterCell& center_cell) const;

    SInt VertexCellSeed(SInt annulus_id, SInt chunk_id, SInt cell_id);

    SampledVertex SampleVertex(Double min_phi, Double max_phi, Double min_cdf, Double max_cdf);

    void AppendVertex(VertexBlock& block, SInt id, const SampledVertex& vertex);
};
template <typename Double>
void Hyper_Hyperbolic<Double>::FinalizeCSR(MPI_Comm /*comm*/) {}

using LowPrecisionHyperHyperbolic  = Hyper_Hyperbolic<LPFloat>;
using HighPrecisionHyperHyperbolic = Hyper_Hyperbolic<HPFloat>;
} // namespace kagen
