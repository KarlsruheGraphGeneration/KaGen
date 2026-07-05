#pragma once

#include "kagen/generators/hyper/h_hyperbolic/circular_interval.h"
#include "kagen/generators/hyper/h_hyperbolic/hyper_hyperbolic.h"
#include "kagen/generators/hyper/h_hyperbolic/poincare_geometry.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
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

    explicit HyperbolicGeometryPolicy(GeneratorT& generator, const SInt annulus_id, const SInt chunk_id)
        : gen_(generator),
          rng_(RNGWrapper(generator.config_)) {
        current_annulus_half_angle_.resize(gen_.total_annuli_, Double{-1.0});
        EnsureAnnulusMidpoints();
    }

    void AddCenter(const Center& /*unused*/, std::vector<SInt>& /*unused*/) const {}
    void CandidateCells(const Center& center, const Double radius, std::vector<Cell>& cells) {
        CacheQueryState(center, radius);

        CandidateCollector collector{*this};
        collector.Collect(center, radius, cells);
    }
    Double Radius(const Center& /*unused*/) const {
        return gen_.current_hyperedge_radius_;
    }

    CellBallRelation ClassifyCell(const Center& /*center*/, const Double /*radius*/, const Cell& cell) const {
        if (CellAABBOutsideBall(cell)) {
            return CellBallRelation::OUTSIDE;
        }

        if (ShouldTryInside(cell) && HyperbolicCellInside(cell)) {
            return CellBallRelation::INSIDE;
        }

        return CellBallRelation::PARTIAL;
    }
    Double CellCoverage(const Center& center, const Double hyperball_radius, const Cell& cell) const {
        const CellBallRelation relation = ClassifyCell(center, hyperball_radius, cell);

        if (relation == CellBallRelation::OUTSIDE) {
            return Double{0.0};
        }

        if (relation == CellBallRelation::INSIDE) {
            return Double{1.0};
        }

        const Double min_phi = cell.min_phi;
        const Double max_phi = cell.max_phi;

        const Double cell_phi_width = max_phi - min_phi;

        if (cell_phi_width <= Double{0.0}) {
            return Double{0.0};
        }

        Double half_angle = current_annulus_half_angle_[cell.annulus_id];

        if (half_angle < Double{0.0}) {
            const auto& mid = annulus_mid_[cell.annulus_id];
            half_angle      = AllowedHalfAngleAtRadius(mid.r, mid.cosh_r, mid.sinh_r);
        }

        Double overlap_phi = Double{0.0};

        if (half_angle >= Double{M_PI}) {
            overlap_phi = cell_phi_width;
        } else if (half_angle > Double{0.0}) {
            overlap_phi = circular_interval::CircularOverlap(
                min_phi, max_phi, center.phi - half_angle, center.phi + half_angle, gen_.cell_eps_);
        }

        return std::clamp(overlap_phi / cell_phi_width, Double{0.0}, Double{1.0});
    }

    SInt AddWholeCell(const Cell& cell, std::vector<PinRange>& ranges) const {
        const auto it = gen_.cells_.find(cell.global_cell_id);
        if (it == gen_.cells_.end()) {
            return 0;
        }

        const auto& stored_cell = it->second;
        const SInt  size        = std::get<0>(stored_cell);
        const SInt  offset      = std::get<4>(stored_cell);

        if (size > 0) {
            ranges.push_back({.begin = offset, .end = offset + size});
        }

        return size;
    }
    SInt AddPartialCellRange(
        const Center& center, const Cell& cell, const Double coverage, std::vector<SInt>& pins,
        std::vector<PinRange>& ranges) const {
        const auto info = GetPartialCellSampleInfo(center, cell, coverage);
        if (!info) {
            return 0;
        }

        auto sampled = getRandomPinRange(info->size, info->k, info->offset, info->seed, gen_.mersenne);
        ranges.insert(ranges.end(), sampled);

        return info->k;
    }
    SInt AddPartialCellFloyd(
        const Center& center, const Cell& cell, const Double coverage, std::vector<SInt>& pins,
        std::vector<PinRange>& ranges) const {
        const auto info = GetPartialCellSampleInfo(center, cell, coverage);
        if (!info) {
            return 0;
        }

        SInt seed = info->seed;

        FloydSampleGeometricAppend(info->offset, info->size, info->k, rng_, seed, pins, floyd_scratch_);
        return info->k;
    }
    SInt AddPartialCellExact(const Center& center, const Double /*radius*/, const Cell& cell, std::vector<SInt>& pins) {
        const auto* vertices = EnsureVerticesForCell(cell);
        if (vertices == nullptr) {
            return 0;
        }

        return gen_.config_.debug ? AddExactVerticesChecked(*vertices, pins) : AddExactVerticesFast(*vertices, pins);
    }

    void EmitHyperedge(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges) {
        gen_.PushHyperedgeCompressed(pins, ranges);
    }

    std::string CenterToString(const Center& center) const {
        std::ostringstream out;

        out << "phi=" << center.phi << ";r=" << center.r;

        return out.str();
    }

private:
    // ===== Query state =====
    void CacheQueryState(const Center& center, const Double radius) {
        center_cosh_r_  = std::cosh(center.r);
        center_sinh_r_  = std::sinh(center.r);
        radius_cosh_    = std::cosh(radius);
        ball_           = poincare_geometry::MakeBall(center.r, center.phi, radius);
        center_cos_phi_ = std::cos(center.phi);
        center_sin_phi_ = std::sin(center.phi);
        center_phi_     = center.phi;
        center_r_       = center.r;
    }

    // ===== Candidate search =====
    struct CandidateCollector {
        HyperbolicGeometryPolicy& policy;
        GeneratorT&               gen() {
            return policy.gen_;
        }

        void Collect(const Center& center, const Double radius, std::vector<Cell>& cells) {
            cells.clear();

            const Double r_min = std::max(Double{0.0}, center.r - radius);
            const Double r_max = std::min(gen().target_r_, center.r + radius);

            if (r_max < r_min) {
                return;
            }

            const Double annulus_width = gen().target_r_ / static_cast<Double>(gen().total_annuli_);

            SInt first_annulus = static_cast<SInt>(std::floor(r_min / annulus_width));
            SInt last_annulus  = static_cast<SInt>(std::floor(r_max / annulus_width));

            first_annulus = std::clamp<SInt>(first_annulus, 0, gen().total_annuli_ - 1);
            last_annulus  = std::clamp<SInt>(last_annulus, 0, gen().total_annuli_ - 1);

            std::fill(
                policy.current_annulus_half_angle_.begin(), policy.current_annulus_half_angle_.end(), Double{-1.0});

            for (SInt annulus_id = first_annulus; annulus_id <= last_annulus; ++annulus_id) {
                AddCandidateCellsInAnnulus(center, annulus_id, cells);
            }
        }

        void AddCandidateCellsInAnnulus(const Center& center, const SInt annulus_id, std::vector<Cell>& cells) {
            const Double half_angle                        = policy.AllowedHalfAngleForAnnulus(center.r, annulus_id);
            policy.current_annulus_half_angle_[annulus_id] = half_angle;

            if (half_angle <= 0.0) {
                return;
            }

            AddCellsInAngularInterval(annulus_id, center.phi - half_angle, center.phi + half_angle, cells);
        }

        void AddCellsInAngularInterval(
            const SInt annulus_id, const Double min_phi, const Double max_phi, std::vector<Cell>& cells) {
            const SInt grid_size   = gen().GridSizeForAnnulus(annulus_id);
            const auto query_parts = circular_interval::Split(min_phi, max_phi, gen().cell_eps_);

            for (SInt chunk_id = 0; chunk_id < gen().config_.k; ++chunk_id) {
                const auto& chunk = gen().chunks_[chunk_id];

                const Double chunk_min_phi = std::get<1>(chunk);
                const Double chunk_max_phi = std::get<2>(chunk);

                bool chunk_overlaps = false;

                for (int p = 0; p < query_parts.count; ++p) {
                    const Double q_begin = query_parts.parts[p].first;
                    const Double q_end   = query_parts.parts[p].second;

                    if (std::max(q_begin, chunk_min_phi) < std::min(q_end, chunk_max_phi)) {
                        chunk_overlaps = true;
                        break;
                    }
                }

                if (!chunk_overlaps) {
                    continue;
                }

                gen().GenerateCells(annulus_id, chunk_id);

                const Double cell_width = (chunk_max_phi - chunk_min_phi) / grid_size;

                for (int p = 0; p < query_parts.count; ++p) {
                    const Double q_begin = query_parts.parts[p].first;
                    const Double q_end   = query_parts.parts[p].second;

                    const Double overlap_begin = std::max(q_begin, chunk_min_phi);
                    const Double overlap_end   = std::min(q_end, chunk_max_phi);

                    if (overlap_begin >= overlap_end) {
                        continue;
                    }

                    SInt first_cell = static_cast<SInt>(std::floor((overlap_begin - chunk_min_phi) / cell_width));

                    SInt last_cell = static_cast<SInt>(std::ceil((overlap_end - chunk_min_phi) / cell_width)) - 1;

                    first_cell = std::clamp<SInt>(first_cell, 0, grid_size - 1);
                    last_cell  = std::clamp<SInt>(last_cell, 0, grid_size - 1);

                    for (SInt cell_id = first_cell; cell_id <= last_cell; ++cell_id) {
                        PushCandidateCell(annulus_id, chunk_id, cell_id, cells);
                    }
                }
            }
        }

        void
        PushCandidateCell(const SInt annulus_id, const SInt chunk_id, const SInt cell_id, std::vector<Cell>& cells) {
            const SInt global_chunk_id = gen().ComputeGlobalChunkId(annulus_id, chunk_id);
            const SInt global_cell_id  = gen().ComputeGlobalCellId(annulus_id, chunk_id, cell_id);

            const auto& annulus     = gen().annuli_.find(global_chunk_id)->second;
            const auto& stored_cell = gen().cells_.find(global_cell_id)->second;

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
    Double AllowedHalfAngleForAnnulus(const Double center_r, const SInt annulus_id) const {
        Double reach = std::max(
            AllowedHalfAngleForCachedRadius(gen_.annulus_min_cosh_[annulus_id], gen_.annulus_min_sinh_[annulus_id]),
            AllowedHalfAngleForCachedRadius(gen_.annulus_max_cosh_[annulus_id], gen_.annulus_max_sinh_[annulus_id]));

        const Double min_r = gen_.annulus_min_r_[annulus_id];
        const Double max_r = gen_.annulus_max_r_[annulus_id];

        if (center_r >= min_r && center_r <= max_r) {
            reach = std::max(reach, AllowedHalfAngleForCachedRadius(center_cosh_r_, center_sinh_r_));
        }

        return reach;
    }
    Double AllowedHalfAngleForCachedRadius(const Double query_cosh, const Double query_sinh) const {
        const Double denom = center_sinh_r_ * query_sinh;

        if (denom <= std::numeric_limits<Double>::epsilon()) {
            return Double{M_PI};
        }

        const Double arg = ((center_cosh_r_ * query_cosh) - radius_cosh_) / denom;

        if (arg <= Double{-1.0}) {
            return Double{M_PI};
        }

        if (arg >= Double{1.0}) {
            return Double{0.0};
        }

        return std::acos(arg);
    }

    // ===== Cell classification =====
    bool CellAABBOutsideBall(const Cell& cell) const {
        const poincare_geometry::AABB<Double> box{
            .min_x = cell.min_x,
            .max_x = cell.max_x,
            .min_y = cell.min_y,
            .max_y = cell.max_y,
        };

        return poincare_geometry::AABBOutsideBall(box, ball_);
    }
    bool HyperbolicCellInside(const Cell& cell) const {
        const Double max_delta = MaxAngularDistanceToInterval(center_phi_, cell.min_phi, cell.max_phi);

        const Double cos_delta = std::cos(max_delta);

        const Double cosh_d_min_r =
            (center_cosh_r_ * std::cosh(cell.min_r)) - (center_sinh_r_ * std::sinh(cell.min_r) * cos_delta);

        if (cosh_d_min_r > radius_cosh_) {
            return false;
        }

        const Double cosh_d_max_r =
            (center_cosh_r_ * std::cosh(cell.max_r)) - (center_sinh_r_ * std::sinh(cell.max_r) * cos_delta);

        return cosh_d_max_r <= radius_cosh_;
    }
    bool ShouldTryInside(const Cell& cell) const {
        const auto it = gen_.cells_.find(cell.global_cell_id);
        if (it == gen_.cells_.end()) {
            return false;
        }

        const SInt n = std::get<0>(it->second);

        if (n < 256) {
            return false;
        }

        // Optional: only if the ball is large enough to plausibly cover whole cells.
        if (gen_.current_hyperedge_radius_ < (cell.max_r - cell.min_r)) {
            return false;
        }

        return true;
    }
    Double MaxAngularDistanceToInterval(const Double phi, const Double min_phi, const Double max_phi) const {
        const Double antipode = circular_interval::NormalizePhi(phi + Double{M_PI});

        if (circular_interval::AngleInInterval(antipode, min_phi, max_phi)) {
            return Double{M_PI};
        }

        return std::max(
            circular_interval::AngularDistance(phi, min_phi), circular_interval::AngularDistance(phi, max_phi));
    }

    // ===== Coverage estimation =====
    Double AllowedHalfAngleAtRadius(const Double query_r, const Double query_cosh, const Double query_sinh) const {
        const Double denominator = center_sinh_r_ * query_sinh;

        if (denominator <= std::numeric_limits<Double>::epsilon()) {
            const Double radial_distance = std::abs(center_r_ - query_r);
            return radial_distance <= gen_.current_hyperedge_radius_ ? M_PI : 0.0;
        }

        Double arg = ((center_cosh_r_ * query_cosh) - radius_cosh_) / denominator;

        if (arg <= -1.0) {
            return M_PI;
        }
        if (arg >= 1.0) {
            return 0.0;
        }

        return std::acos(arg);
    }

    // ===== Angular / distance helpers =====
    Double MinAngularDistanceToInterval(Double phi, Double min_phi, Double max_phi) const {
        if (circular_interval::AngleInInterval(phi, min_phi, max_phi)) {
            return 0.0;
        }

        return std::min(
            circular_interval::AngularDistance(phi, min_phi), circular_interval::AngularDistance(phi, max_phi));
    }
    Double CoshDistanceWithDelta(Double query_r, Double delta_phi) const {
        return (center_cosh_r_ * std::cosh(query_r)) - (center_sinh_r_ * std::sinh(query_r) * std::cos(delta_phi));
    }

    // ===== Approx Partial-Cell Sampling =====
    struct PartialCellSampleInfo {
        SInt global_cell_id;
        SInt size;
        SInt offset;
        SInt k;
        SInt seed;
    };
    std::optional<PartialCellSampleInfo>
    GetPartialCellSampleInfo(const Center& center, const Cell& cell, const Double coverage) const {
        const SInt global_cell_id = cell.global_cell_id;

        const auto it = gen_.cells_.find(global_cell_id);
        if (it == gen_.cells_.end()) {
            return std::nullopt;
        }

        const auto& stored_cell = it->second;

        const SInt size   = std::get<0>(stored_cell);
        const SInt offset = std::get<4>(stored_cell);
        const SInt k      = static_cast<SInt>(std::floor(static_cast<Double>(size) * coverage));

        if (k <= 0) {
            return std::nullopt;
        }

        const SInt seed =
            sampling::Spooky::hash(gen_.config_.seed + (131 * center.sampled_id) + (9973 * global_cell_id));

        return PartialCellSampleInfo{
            .global_cell_id = global_cell_id,
            .size           = size,
            .offset         = offset,
            .k              = k,
            .seed           = seed,
        };
    }

    // ===== Exact vertex checks =====
    const typename GeneratorT::VertexBlock* EnsureVerticesForCell(const Cell& cell) {
        auto it = gen_.vertices_.find(cell.global_cell_id);

        if (it == gen_.vertices_.end()) {
            gen_.GenerateVertices(cell.annulus_id, cell.chunk_id, cell.cell_id);
            it = gen_.vertices_.find(cell.global_cell_id);
        }

        if (it == gen_.vertices_.end()) {
            return nullptr;
        }

        return &it->second;
    }

    bool IsInsideHyperbolicBallFast(
        const Double v_cosh_r, const Double v_sinh_r, const Double v_cos_phi, const Double v_sin_phi) const {
        const Double term_a       = center_cosh_r_ * v_cosh_r;
        const Double term_b_scale = center_sinh_r_ * v_sinh_r;
        const Double cos_delta    = (center_cos_phi_ * v_cos_phi) + (center_sin_phi_ * v_sin_phi);

        const Double pdm_fast = (term_a - (term_b_scale * cos_delta) - Double{1.0}) / Double{2.0};

        return pdm_fast <= gen_.current_hyperedge_pdm_radius_;
    }
    SInt AddExactVerticesFast(const typename GeneratorT::VertexBlock& vertices, std::vector<SInt>& pins) const {
        SInt accepted = 0;

        for (std::size_t i = 0; i < vertices.id.size(); ++i) {
            if (IsInsideHyperbolicBallFast(
                    vertices.cosh_r[i], vertices.sinh_r[i], vertices.cos_phi[i], vertices.sin_phi[i])) {
                pins.push_back(vertices.id[i]);
                ++accepted;
            }
        }

        return accepted;
    }
    bool IsInsideHyperbolicBallChecked(
        const Double phi, const Double v_cosh_r, const Double v_sinh_r, const Double v_cos_phi,
        const Double v_sin_phi) const {
        const Double term_a       = center_cosh_r_ * v_cosh_r;
        const Double term_b_scale = center_sinh_r_ * v_sinh_r;
        const Double cos_delta    = (center_cos_phi_ * v_cos_phi) + (center_sin_phi_ * v_sin_phi);

        const Double pdm_fast = (term_a - (term_b_scale * cos_delta) - Double{1.0}) / Double{2.0};

        constexpr Double rel_eps = Double{4096.0} * std::numeric_limits<Double>::epsilon();

        const Double expr_scale = std::max<Double>(
            Double{1.0}, std::abs(term_a) + std::abs(term_b_scale) + std::abs(gen_.current_hyperedge_pdm_radius_));

        const Double tol = rel_eps * expr_scale;

        if (std::abs(pdm_fast - gen_.current_hyperedge_pdm_radius_) <= tol) {
            return HyperbolicDistanceReference(phi, v_cosh_r, v_sinh_r) <= gen_.current_hyperedge_pdm_radius_;
        }

        return pdm_fast <= gen_.current_hyperedge_pdm_radius_;
    }
    SInt AddExactVerticesChecked(const typename GeneratorT::VertexBlock& vertices, std::vector<SInt>& pins) const {
        SInt accepted = 0;

        for (std::size_t i = 0; i < vertices.id.size(); ++i) {
            if (IsInsideHyperbolicBallChecked(
                    vertices.phi[i], vertices.cosh_r[i], vertices.sinh_r[i], vertices.cos_phi[i],
                    vertices.sin_phi[i])) {
                pins.push_back(vertices.id[i]);
                ++accepted;
            }
        }

        return accepted;
    }
    Double HyperbolicDistanceReference(const Double phi, const Double cosh_r, const Double sinh_r) const {
        Double delta_phi = std::abs(center_phi_ - phi);
        delta_phi        = std::min(delta_phi, Double{2.0 * M_PI} - delta_phi);

        const Double cosh_dist = (center_cosh_r_ * cosh_r) - (center_sinh_r_ * sinh_r * std::cos(delta_phi));

        return (cosh_dist - Double{1.0}) / Double{2.0};
    }

    // ===== Initialization =====
    void EnsureAnnulusMidpoints() {
        if (!annulus_mid_.empty()) {
            return;
        }

        annulus_mid_.resize(gen_.total_annuli_);

        for (SInt i = 0; i < gen_.total_annuli_; ++i) {
            const Double r = (gen_.annulus_min_r_[i] + gen_.annulus_max_r_[i]) * 0.5;

            annulus_mid_[i] = {
                .r      = r,
                .cosh_r = std::cosh(r),
                .sinh_r = std::sinh(r),
            };
        }
    }

    // ===== Data members =====
    GeneratorT&          gen_;
    mutable RNGWrapper<> rng_;

    // Caches
    Double                           center_cosh_r_{};
    Double                           center_sinh_r_{};
    Double                           radius_cosh_{};
    poincare_geometry::Ball<Double>  ball_{};
    Double                           center_cos_phi_{};
    Double                           center_sin_phi_{};
    Double                           center_phi_{};
    Double                           center_r_{};
    mutable std::unordered_set<SInt> floyd_scratch_;
    mutable std::vector<Double>      current_annulus_half_angle_;

    struct AnnulusMidpoint {
        Double r;
        Double cosh_r;
        Double sinh_r;
    };

    std::vector<AnnulusMidpoint> annulus_mid_;
};

} // namespace kagen
