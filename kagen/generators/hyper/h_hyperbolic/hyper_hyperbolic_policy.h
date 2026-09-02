#pragma once

#include "kagen/generators/hyper/h_hyperbolic/circular_interval.h"
#include "kagen/generators/hyper/h_hyperbolic/hyper_hyperbolic.h"
#include "kagen/generators/hyper/h_hyperbolic/poincare_geometry.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"
#include "kagen/tools/geometry.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <list>
#include <optional>
#include <string>
#include <vector>

namespace kagen {

template <typename Double>
class HyperbolicGeometryPolicy {
public:
    using GeneratorT = Hyper_Hyperbolic<Double>;
    using Center     = HyperbolicHyperedgeCenter<Double>;
    using Vertex     = typename GeneratorT::Vertex;
    struct Cell {
        SInt annulus_id;
        SInt chunk_id;
        SInt cell_id;
        SInt global_cell_id;

        Double min_r;
        Double max_r;
        Double min_phi;
        Double max_phi;

        Double min_x;
        Double max_x;
        Double min_y;
        Double max_y;
    };
    struct CellRegion {
        SInt first_annulus;
        SInt last_annulus; // inclusive

        Double min_r;
        Double max_r;

        Double min_phi;
        Double max_phi;

        Double min_x;
        Double max_x;
        Double min_y;
        Double max_y;
    };

    explicit HyperbolicGeometryPolicy(GeneratorT& generator) : gen_(generator), rng_(RNGWrapper(generator.config_)) {
        current_annulus_half_angle_.resize(gen_.total_annuli_, Double{-1.0});
        EnsureAnnulusMidpoints();
    }

    void AddCenter(
        const Center& center, std::vector<SInt>& /*pins*/
    ) const;

    bool HierarchicalCandidateCells(
        const Center& center, Double radius, std::vector<Cell>& cells, std::vector<PinRange>& ranges);

    void CandidateCells(const Center& center, Double radius, std::vector<Cell>& cells);

    Double Radius(const Center& /*unused*/) const;

    CellBallRelation ClassifyCell(const Center& /*center*/, Double /*radius*/, const Cell& cell) const;

    CellBallRelation ClassifyRegion(const CellRegion& region) const;

    CellRegion MakeRegion(SInt first_annulus, SInt last_annulus, Double min_phi, Double max_phi) const;

    std::pair<CellRegion, CellRegion> SplitRegionRadially(const CellRegion& region) const;

    std::pair<CellRegion, CellRegion> SplitRegionAngularly(const CellRegion& region) const;

    bool IsLeaf(const CellRegion& region) const;

    Double CellCoverage(const Center& center, Double /*hyperball_radius*/, const Cell& cell) const;

    SInt AddWholeCell(const Cell& cell, std::vector<PinRange>& ranges) const;

    SInt AddPartialCellRange(
        const Center& /*center*/, const Cell& cell, const Double coverage, std::vector<SInt>& /*pins*/,
        std::vector<PinRange>& ranges) const;

    SInt AddPartialCellFloyd(
        const Center& /*center*/, const Cell& cell, const Double coverage, std::vector<SInt>& pins,
        std::vector<PinRange>& /*ranges*/) const;

    SInt
    AddPartialCellExact(const Center& /*center*/, Double /*radius*/, const Cell& cell, std::vector<SInt>& pins) const;

    void EmitHyperedge(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges);

    std::string CenterToString(const Center& center) const;

    bool ShouldApproximatePartialCell(const Cell& cell) const;

private:
    // ===== Query state =====
    void CacheQueryState(const Center& center, Double radius);

    // ===== Candidate search =====

    struct CellAnnulusRegion {
        SInt annulus_id;

        // Half-open global cell interval [first_cell, end_cell)
        SInt first_cell;
        SInt end_cell;

        Double min_r;
        Double max_r;

        Double min_phi;
        Double max_phi;

        Double min_x;
        Double max_x;
        Double min_y;
        Double max_y;
    };
    struct CandidateCollector {
        HyperbolicGeometryPolicy& policy;
        std::unordered_set<SInt>  seen_candidate_cells_;
        explicit CandidateCollector(HyperbolicGeometryPolicy& policy) : policy(policy) {}

        void AddLeafCell(const CellRegion& region, std::vector<Cell>& cells) {
            const SInt annulus_id = region.first_annulus;

            const auto [first_cell, last_cell] =
                gen().GlobalCellRangeForAngularInterval(annulus_id, region.min_phi, region.max_phi);

            for (SInt global_cell = first_cell; global_cell <= last_cell; ++global_cell) {
                PushGlobalCell(annulus_id, global_cell, cells);
            }
        }

        void AddLeafCell(const CellAnnulusRegion& region, std::vector<Cell>& cells) {
            assert(region.end_cell - region.first_cell == 1);

            PushGlobalCell(region.annulus_id, region.first_cell, cells);
        }

        GeneratorT& gen() {
            return policy.gen_;
        }

        void Collect(
            const Center& center, const Double radius, std::vector<Cell>& cells, std::vector<PinRange>& inside_ranges) {
            cells.clear();
            inside_ranges.clear();
            seen_candidate_cells_.clear();

            std::fill(
                policy.current_annulus_half_angle_.begin(), policy.current_annulus_half_angle_.end(), Double{-1.0});

#ifdef KAGEN_ENABLE_HIERARCHICAL_CELLS
            CollectHierarchical(center, radius, cells, inside_ranges);
#else
            CollectFlat(center, radius, cells);
#endif
        }

        void CollectHierarchical(
            const Center& /*center*/, const Double /*radius*/, std::vector<Cell>& cells,
            std::vector<PinRange>& ranges) {
            seen_candidate_cells_.clear();

            std::fill(
                policy.current_annulus_half_angle_.begin(), policy.current_annulus_half_angle_.end(), Double{-1.0});

            for (SInt annulus_id = 0; annulus_id < gen().total_annuli_; ++annulus_id) {
                const SInt total_cells = gen().global_cells_per_annulus_[annulus_id];

                if (total_cells == 0) {
                    continue;
                }

                const CellAnnulusRegion root = policy.MakeCellAnnulusRegion(annulus_id, 0, total_cells);

                policy.TraverseCandidateRegion(root, cells, ranges, *this);
            }
        }

        void CollectFlat(const Center& center, const Double radius, std::vector<Cell>& cells) {
            cells.clear();
            seen_candidate_cells_.clear();

            std::fill(
                policy.current_annulus_half_angle_.begin(), policy.current_annulus_half_angle_.end(), Double{-1.0});

            if (gen().total_annuli_ <= 0) {
                return;
            }

            //
            // Start in the annulus containing the hyperball center.
            //
            const Double annulus_width = gen().target_r_ / static_cast<Double>(gen().total_annuli_);

            SInt center_annulus = static_cast<SInt>(std::floor(center.r / annulus_width));

            center_annulus = std::clamp<SInt>(center_annulus, 0, gen().total_annuli_ - 1);

            //
            // Center annulus.
            //
            AddCandidateCellsInAnnulus(center, center_annulus, cells);

            //
            // Walk outward.
            //
            for (SInt annulus_id = center_annulus + 1; annulus_id < gen().total_annuli_; ++annulus_id) {
                const Double min_r = gen().annulus_min_r_[annulus_id];

                //
                // Closest possible point in this annulus is already
                // farther away radially than the hyperball radius.
                //
                if (min_r - center.r > radius) {
                    break;
                }

                AddCandidateCellsInAnnulus(center, annulus_id, cells);
            }

            //
            // Walk inward.
            //
            for (SInt annulus_id = center_annulus; annulus_id > 0;) {
                --annulus_id;

                const Double max_r = gen().annulus_max_r_[annulus_id];

                if (center.r - max_r > radius) {
                    break;
                }

                AddCandidateCellsInAnnulus(center, annulus_id, cells);
            }
        }

        void AddOccupiedCellsInAngularInterval(
            const SInt annulus_id, const Double min_phi, const Double max_phi, std::vector<Cell>& cells) {
            const auto parts = circular_interval::Split(min_phi, max_phi, gen().cell_eps_);

            for (int p = 0; p < parts.count; ++p) {
                const Double q_begin = parts.parts[p].first;

                const Double q_end = parts.parts[p].second;

                if (!(q_begin < q_end)) {
                    continue;
                }

                const auto [first, last] = gen().GlobalCellRangeForAngularInterval(annulus_id, q_begin, q_end);

                for (SInt global_cell = first; global_cell <= last; ++global_cell) {
                    PushGlobalCell(annulus_id, global_cell, cells);
                }
            }
        }

        void AddCandidateCellsInAnnulus(const Center& center, const SInt annulus_id, std::vector<Cell>& cells) {
            const Double half_angle = policy.AllowedHalfAngleForAnnulus(center.r, annulus_id);

            policy.current_annulus_half_angle_[annulus_id] = half_angle;

            if (!(half_angle > Double{0.0})) {
                return;
            }

            //
            // Entire angular extent of this annulus.
            //
            if (half_angle >= Double{M_PI}) {
                const SInt total_cells = gen().global_cells_per_annulus_[annulus_id];

                for (SInt global_cell = 0; global_cell < total_cells; ++global_cell) {
                    PushGlobalCell(annulus_id, global_cell, cells);
                }

                return;
            }

            const Double center_phi = circular_interval::NormalizePhi(center.phi);

            AddOccupiedCellsInAngularInterval(annulus_id, center_phi - half_angle, center_phi + half_angle, cells);
        }
        void PushGlobalCell(const SInt annulus_id, const SInt global_cell, std::vector<Cell>& cells) {
            const auto [chunk_id, local_cell_id] = gen().GlobalCellToChunkCell(annulus_id, global_cell);

            const SInt global_cell_id = gen().ComputeGlobalCellId(annulus_id, chunk_id, local_cell_id);

            const CellBlock& block = policy.GetCellBlock(annulus_id, chunk_id);

            if (local_cell_id < 0 || static_cast<std::size_t>(local_cell_id) >= block.cells.size()) {
                return;
            }

            const auto& stored_cell = block.cells[static_cast<std::size_t>(local_cell_id)];

            if (std::get<0>(stored_cell) <= 0) {
                return;
            }

            if (!seen_candidate_cells_.insert(global_cell_id).second) {
                return;
            }

            PushCandidateCell(annulus_id, chunk_id, local_cell_id, cells);
        }

        void AddCellsInAngularInterval(
            const SInt annulus_id, const Double min_phi, const Double max_phi, std::vector<Cell>& cells) {
            const auto query_parts = circular_interval::Split(min_phi, max_phi, gen().cell_eps_);

            for (int part = 0; part < query_parts.count; ++part) {
                const Double q_begin = query_parts.parts[part].first;

                const Double q_end = query_parts.parts[part].second;

                if (!(q_begin < q_end)) {
                    continue;
                }

                const auto [first_global_cell, last_global_cell] =
                    gen().GlobalCellRangeForAngularInterval(annulus_id, q_begin, q_end);

                for (SInt global_cell = first_global_cell; global_cell <= last_global_cell; ++global_cell) {
                    PushGlobalCell(annulus_id, global_cell, cells);
                }
            }
        }

        void
        PushCandidateCell(const SInt annulus_id, const SInt chunk_id, const SInt cell_id, std::vector<Cell>& cells) {
            const SInt global_cell_id = gen().ComputeGlobalCellId(annulus_id, chunk_id, cell_id);

            const CellBlock& block   = policy.GetCellBlock(annulus_id, chunk_id);
            const auto&      annulus = block.annulus;

            assert(cell_id >= 0);
            assert(static_cast<std::size_t>(cell_id) < block.cells.size());

            const auto& stored_cell = block.cells[static_cast<std::size_t>(cell_id)];

            assert(std::get<0>(stored_cell) > 0);
            const Double min_r   = std::get<1>(annulus);
            const Double max_r   = std::get<2>(annulus);
            const Double min_phi = std::get<1>(stored_cell);
            const Double max_phi = std::get<2>(stored_cell);

            cells.push_back(
                Cell{
                    .annulus_id     = annulus_id,
                    .chunk_id       = chunk_id,
                    .cell_id        = cell_id,
                    .global_cell_id = global_cell_id,
                    .min_r          = min_r,
                    .max_r          = max_r,
                    .min_phi        = min_phi,
                    .max_phi        = max_phi,
                    .min_x          = std::get<5>(stored_cell),
                    .max_x          = std::get<6>(stored_cell),
                    .min_y          = std::get<7>(stored_cell),
                    .max_y          = std::get<8>(stored_cell),
                });
        }
    };

    CellAnnulusRegion MakeCellAnnulusRegion(SInt annulus_id, SInt first_cell, SInt end_cell) const;

    std::pair<CellAnnulusRegion, CellAnnulusRegion> SplitCellAnnulusRegion(const CellAnnulusRegion& region) const;

    bool IsLeaf(const CellAnnulusRegion& region) const;

    CellBallRelation ClassifyRegion(const CellAnnulusRegion& region) const;

    void EmitInsideRegion(const CellAnnulusRegion& region, std::vector<PinRange>& inside_ranges) const;

    void TraverseCandidateRegion(
        const CellAnnulusRegion& region, std::vector<Cell>& cells, std::vector<PinRange>& inside_ranges,
        CandidateCollector& collector);

    void TraverseCandidateRegion(
        const CellRegion& region, std::vector<Cell>& cells, std::vector<PinRange>& inside_ranges,
        CandidateCollector& collector);

    Double AllowedHalfAngleForAnnulus(Double /*center_r*/, SInt annulus_id) const;

    Double AllowedHalfAngleForCachedRadius(Double query_cosh, Double query_sinh) const;

    // ===== Cell classification =====
    bool CellAABBOutsideBall(const Cell& cell) const;

    bool HyperbolicCellInside(const Cell& cell) const;

    bool ShouldTryInside(const Cell& cell) const;

    Double MaxAngularDistanceToInterval(Double phi, Double min_phi, Double max_phi) const;

    // ===== Coverage estimation =====
    Double AllowedHalfAngleAtRadius(Double query_r, Double query_cosh, Double query_sinh) const;

    // ===== Angular / distance helpers =====
    Double MinAngularDistanceToInterval(Double phi, Double min_phi, Double max_phi) const;

    Double CoshDistanceWithDelta(Double query_r, Double delta_phi) const;

    // ===== Approx Partial-Cell Sampling =====
    struct PartialCellSampleInfo {
        SInt global_cell_id;
        SInt size;
        SInt offset;
        SInt k;
    };
    std::optional<PartialCellSampleInfo> GetPartialCellSampleInfo(const Cell& cell, Double coverage) const;

    // ===== Exact vertex checks =====
    struct CachedExactCell {
        typename GeneratorT::VertexBlock vertices;
    };

    const Hyper_Hyperbolic<Double>::VertexBlock& ExactVertices(const Cell& cell) const;

    const CachedExactCell& ExactRemoteCell(const Cell& cell) const;

    bool IsLocalCell(const Cell& cell) const;
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    void RecordRemoteAccess(const SInt global_cell_id) const;

public:
    void PrintExactCacheStats() const;

private:
#endif

    std::size_t ExactCellBytes(const GeneratorT::VertexBlock& vertices) const;

    bool IsInsideHyperbolicBallFast(const Vertex& vertex) const;

    SInt AddExactVerticesFast(const typename GeneratorT::VertexBlock& vertices, std::vector<SInt>& pins) const;

    SInt AddExactVerticesChecked(const typename GeneratorT::VertexBlock& vertices, std::vector<SInt>& pins) const;

    struct CellBlock {
        SInt                                                 annulus_id;
        SInt                                                 chunk_id;
        std::vector<typename Hyper_Hyperbolic<Double>::Cell> cells;
        typename GeneratorT::Annulus                         annulus;
    };
    struct CellBlockKey {
        SInt annulus_id;
        SInt chunk_id;

        bool operator==(const CellBlockKey& other) const {
            return annulus_id == other.annulus_id && chunk_id == other.chunk_id;
        }
    };

    struct CellBlockKeyHash {
        std::size_t operator()(const CellBlockKey& key) const noexcept {
            std::size_t hash = std::hash<SInt>{}(key.annulus_id);
            hash ^= std::hash<SInt>{}(key.chunk_id) + 0x9e3779b9U + (hash << 6U) + (hash >> 2U);
            return hash;
        }
    };

    using CellBlockLRU = std::list<CellBlockKey>;

    struct CachedCellBlock {
        CellBlock              block;
        CellBlockLRU::iterator lru_position;
        std::size_t            bytes;
    };

    const CellBlock& GetCellBlock(SInt annulus_id, SInt chunk_id) const;

    CellBlock BuildCellBlock(SInt annulus_id, SInt chunk_id) const;

    std::size_t CellBlockBytes(const CellBlock& block) const;

    // ===== Initialization =====
    void EnsureAnnulusMidpoints();

    // ===== Data members =====
    GeneratorT&          gen_;
    mutable RNGWrapper<> rng_;

    // Caches
    Double                          center_cosh_r_{};
    Double                          center_sinh_r_{};
    Double                          radius_cosh_{};
    poincare_geometry::Ball<Double> ball_{};
    Double                          center_cos_phi_{};
    Double                          center_sin_phi_{};
    Double                          center_phi_{};
    Double                          center_r_{};
    Double                          center_x_{};
    Double                          center_y_{};
    Vertex                          center_vertex_{};
    Double                          center_gamma_{};

    mutable std::unordered_set<SInt>                                            floyd_scratch_;
    mutable std::vector<Double>                                                 current_annulus_half_angle_;
    mutable std::unordered_map<CellBlockKey, CachedCellBlock, CellBlockKeyHash> cell_blocks_;
    mutable CellBlockLRU                                                        cell_block_lru_;
    mutable std::size_t                                                         cell_block_cached_bytes_ = 0;

    static constexpr std::size_t cell_block_cache_budget_ = 64ULL * 1024ULL * 1024ULL;

    struct AnnulusMidpoint {
        Double r;
        Double cosh_r;
        Double sinh_r;
    };

    std::vector<AnnulusMidpoint>                                annulus_mid_;
    mutable std::unordered_map<SInt, CachedExactCell>           exact_vertices_;
    mutable std::list<SInt>                                     exact_lru_;
    mutable std::unordered_map<SInt, std::list<SInt>::iterator> exact_lru_pos_;
    mutable std::size_t                                         exact_remote_cached_bytes_ = 0;
    static constexpr std::size_t                                exact_remote_cache_budget_ = 256ULL * 1024ULL * 1024ULL;

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    mutable SInt exact_remote_cache_hits_      = 0;
    mutable SInt exact_remote_cache_misses_    = 0;
    mutable SInt exact_remote_cached_vertices_ = 0;

    mutable SInt                           exact_remote_access_counter_ = 0;
    mutable std::unordered_map<SInt, SInt> exact_remote_last_access_;

    mutable SInt exact_remote_reuse_count_        = 0;
    mutable SInt exact_remote_reuse_distance_sum_ = 0;
    mutable SInt exact_remote_reuse_distance_max_ = 0;

    mutable SInt exact_remote_reuse_distance_le_1_  = 0;
    mutable SInt exact_remote_reuse_distance_le_4_  = 0;
    mutable SInt exact_remote_reuse_distance_le_16_ = 0;
    mutable SInt exact_remote_reuse_distance_gt_16_ = 0;
#endif
};

} // namespace kagen