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

    CellBallRelation ClassifyCell(const Center& center, const Double radius, const Cell& cell) const {
        if (CellAABBOutsideBall(cell)) {
            return CellBallRelation::OUTSIDE;
        }

        if (!IsLocalCell(cell)) {
            return CellBallRelation::PARTIAL;
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
    SInt AddPartialCellExact(
        const Center& center, const Double /*radius*/, const Cell& cell, std::vector<SInt>& pins) const {
        const auto& vertices = ExactVertices(cell);

        return gen_.config_.debug ? AddExactVerticesChecked(vertices, pins) : AddExactVerticesFast(vertices, pins);
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
            const SInt grid_size  = gen().GridSizeForAnnulus(annulus_id);
            const SInt num_chunks = gen().config_.k;

            if (num_chunks <= 0 || grid_size <= 0) {
                return;
            }

            const Double two_pi = Double{2.0 * M_PI};

            const Double chunk_width = two_pi / static_cast<Double>(num_chunks);

            const auto query_parts = circular_interval::Split(min_phi, max_phi, gen().cell_eps_);

            for (int part = 0; part < query_parts.count; ++part) {
                const Double q_begin = query_parts.parts[part].first;

                const Double q_end = query_parts.parts[part].second;

                if (!(q_begin < q_end)) {
                    continue;
                }

                // q_end is treated as exclusive.
                const Double q_end_inside = std::nextafter(q_end, q_begin);

                SInt first_chunk = static_cast<SInt>(std::floor(q_begin / chunk_width));

                SInt last_chunk = static_cast<SInt>(std::floor(q_end_inside / chunk_width));

                first_chunk = std::clamp<SInt>(first_chunk, 0, num_chunks - 1);

                last_chunk = std::clamp<SInt>(last_chunk, 0, num_chunks - 1);

                for (SInt chunk_id = first_chunk; chunk_id <= last_chunk; ++chunk_id) {
                    if (gen().chunks_.find(chunk_id) == std::end(gen().chunks_)) {
                        gen().ComputeChunk(chunk_id);
                        gen().ComputeAnnuli(chunk_id);
                    }

                    const auto& chunk = gen().chunks_[chunk_id];

                    const Double chunk_min_phi = std::get<1>(chunk);

                    const Double chunk_max_phi = std::get<2>(chunk);

                    const Double overlap_begin = std::max(q_begin, chunk_min_phi);

                    const Double overlap_end = std::min(q_end, chunk_max_phi);

                    if (!(overlap_begin < overlap_end)) {
                        continue;
                    }

                    gen().GenerateCells(annulus_id, chunk_id);

                    const Double cell_width = (chunk_max_phi - chunk_min_phi) / static_cast<Double>(grid_size);

                    if (!(cell_width > Double{0.0})) {
                        continue;
                    }

                    SInt first_cell = static_cast<SInt>(std::floor((overlap_begin - chunk_min_phi) / cell_width));

                    const Double overlap_end_inside = std::nextafter(overlap_end, overlap_begin);

                    SInt last_cell = static_cast<SInt>(std::floor((overlap_end_inside - chunk_min_phi) / cell_width));

                    first_cell = std::clamp<SInt>(first_cell, 0, grid_size - 1);

                    last_cell = std::clamp<SInt>(last_cell, 0, grid_size - 1);

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
    struct CachedExactCell {
        typename GeneratorT::VertexBlock vertices;
    };

    const typename GeneratorT::VertexBlock& ExactVertices(const Cell& cell) const {
        if (IsLocalCell(cell)) {
            gen_.GenerateVertices(cell.annulus_id, cell.chunk_id, cell.cell_id);

            static const typename GeneratorT::VertexBlock empty;

            const auto it = gen_.vertices_.find(cell.global_cell_id);
            if (it == gen_.vertices_.end()) {
                return empty;
            }

            return it->second;
        }

        return ExactRemoteCell(cell).vertices;
    }

    const CachedExactCell& ExactRemoteCell(const Cell& cell) const {
        RecordRemoteAccess(cell.global_cell_id);

        auto it = exact_vertices_.find(cell.global_cell_id);
        if (it != exact_vertices_.end()) {
            ++exact_remote_cache_hits_;
            exact_lru_.splice(exact_lru_.begin(), exact_lru_, exact_lru_pos_[cell.global_cell_id]);
            return it->second;
        }

        ++exact_remote_cache_misses_;

        auto [inserted_it, inserted] = exact_vertices_.emplace(cell.global_cell_id, CachedExactCell{});

        exact_lru_.push_front(cell.global_cell_id);
        exact_lru_pos_[cell.global_cell_id] = exact_lru_.begin();

        auto& cached = inserted_it->second;

        gen_.GenerateVertices(cell.annulus_id, cell.chunk_id, cell.cell_id, cached.vertices);

        exact_remote_cached_vertices_ += static_cast<SInt>(cached.vertices.size());

        while (exact_vertices_.size() > exact_remote_cache_limit_) {
            const SInt victim = exact_lru_.back();

            const auto victim_it = exact_vertices_.find(victim);
            if (victim_it != exact_vertices_.end()) {
                exact_remote_cached_vertices_ -= static_cast<SInt>(victim_it->second.vertices.size());
                exact_vertices_.erase(victim_it);
            }

            exact_lru_.pop_back();
            exact_lru_pos_.erase(victim);
        }

        return cached;
    }

    bool IsLocalCell(const Cell& cell) const {
        return gen_.IsLocalChunk(cell.chunk_id);
    }

    void RecordRemoteAccess(const SInt global_cell_id) const {
        const SInt t = ++exact_remote_access_counter_;

        auto it = exact_remote_last_access_.find(global_cell_id);
        if (it == exact_remote_last_access_.end()) {
            exact_remote_last_access_[global_cell_id] = t;
            return;
        }

        const SInt distance = t - it->second;
        it->second          = t;

        ++exact_remote_reuse_count_;
        exact_remote_reuse_distance_sum_ += distance;
        exact_remote_reuse_distance_max_ = std::max(exact_remote_reuse_distance_max_, distance);

        if (distance <= 1) {
            ++exact_remote_reuse_distance_le_1_;
        } else if (distance <= 4) {
            ++exact_remote_reuse_distance_le_4_;
        } else if (distance <= 16) {
            ++exact_remote_reuse_distance_le_16_;
        } else {
            ++exact_remote_reuse_distance_gt_16_;
        }
    }

public:
    void PrintExactCacheStats() const {
        std::cerr << " exact_remote_hits=" << exact_remote_cache_hits_
                  << " exact_remote_misses=" << exact_remote_cache_misses_
                  << " exact_remote_cached_vertices=" << exact_remote_cached_vertices_ << '\n';

        const double avg_reuse =
            (exact_remote_reuse_count_ != 0u)
                ? static_cast<double>(exact_remote_reuse_distance_sum_) / static_cast<double>(exact_remote_reuse_count_)
                : 0.0;

        std::cerr << " remote_reuse_count=" << exact_remote_reuse_count_ << " remote_reuse_avg=" << avg_reuse
                  << " remote_reuse_max=" << exact_remote_reuse_distance_max_
                  << " reuse<=1=" << exact_remote_reuse_distance_le_1_
                  << " reuse<=4=" << exact_remote_reuse_distance_le_4_
                  << " reuse<=16=" << exact_remote_reuse_distance_le_16_
                  << " reuse>16=" << exact_remote_reuse_distance_gt_16_ << '\n';
    }

private:
    bool IsInsideHyperbolicBallFast(
        const Double v_r, const Double v_sinh_r, const Double v_cos_phi, const Double v_sin_phi) const {
        // Radial contribution:
        //
        // (cosh(center_r - v_r) - 1) / 2
        //     = sinh^2((center_r - v_r) / 2)
        const Double half_radial_difference = (center_r_ - v_r) / Double{2.0};

        const Double sinh_half_radial_difference = std::sinh(half_radial_difference);

        const Double radial_term = sinh_half_radial_difference * sinh_half_radial_difference;

        // Angular contribution:
        //
        // sin^2(delta_phi / 2)
        //     = ((cos(center_phi) - cos(vertex_phi))^2
        //        + (sin(center_phi) - sin(vertex_phi))^2) / 4
        //
        // This is more accurate than computing 1 - cos(delta_phi) when
        // delta_phi is very small.
        const Double delta_cos = center_cos_phi_ - v_cos_phi;
        const Double delta_sin = center_sin_phi_ - v_sin_phi;

        const Double sin_half_delta_sq = ((delta_cos * delta_cos) + (delta_sin * delta_sin)) / Double{4.0};

        const Double angular_term = center_sinh_r_ * v_sinh_r * sin_half_delta_sq;

        const Double pdm_distance = radial_term + angular_term;

        return pdm_distance <= gen_.current_hyperedge_pdm_radius_;
    }
    SInt AddExactVerticesFast(const typename GeneratorT::VertexBlock& vertices, std::vector<SInt>& pins) const {
        SInt accepted = 0;

        for (std::size_t i = 0; i < vertices.id.size(); ++i) {
            if (IsInsideHyperbolicBallFast(
                    vertices.r[i], vertices.sinh_r[i], vertices.cos_phi[i], vertices.sin_phi[i])) {
                pins.push_back(vertices.id[i]);
                ++accepted;
            }
        }

        return accepted;
    }
    bool IsInsideHyperbolicBallChecked(
        const Double phi, const Double v_r, const Double v_cosh_r, const Double v_sinh_r, const Double v_cos_phi,
        const Double v_sin_phi) const {
        const bool stable_inside = IsInsideHyperbolicBallFast(v_r, v_sinh_r, v_cos_phi, v_sin_phi);

#ifndef NDEBUG
        const bool reference_inside =
            HyperbolicDistanceReference(phi, v_cosh_r, v_sinh_r) <= gen_.current_hyperedge_pdm_radius_;

        if (stable_inside != reference_inside) {
            std::cerr << std::setprecision(std::numeric_limits<Double>::max_digits10)
                      << "STABLE/REFERENCE DISAGREEMENT\n"
                      << "  center_r=" << center_r_ << '\n'
                      << "  center_phi=" << center_phi_ << '\n'
                      << "  vertex_r=" << v_r << '\n'
                      << "  vertex_phi=" << phi << '\n'
                      << "  radius=" << gen_.current_hyperedge_radius_ << '\n'
                      << "  stable_inside=" << stable_inside << '\n'
                      << "  reference_inside=" << reference_inside << '\n';
        }
#endif

        return stable_inside;
    }
    SInt AddExactVerticesChecked(const typename GeneratorT::VertexBlock& vertices, std::vector<SInt>& pins) const {
        SInt accepted = 0;

        for (std::size_t i = 0; i < vertices.id.size(); ++i) {
            if (IsInsideHyperbolicBallChecked(
                    vertices.phi[i], vertices.r[i], vertices.cosh_r[i], vertices.sinh_r[i], vertices.cos_phi[i],
                    vertices.sin_phi[i])) {
                pins.push_back(vertices.id[i]);
                ++accepted;
            }

            if (gen_.config_.debug) {
                const bool fast_inside = IsInsideHyperbolicBallFast(
                    vertices.cosh_r[i], vertices.sinh_r[i], vertices.cos_phi[i], vertices.sin_phi[i]);

                const bool checked_inside = IsInsideHyperbolicBallChecked(
                    vertices.phi[i], vertices.r[i], vertices.cosh_r[i], vertices.sinh_r[i], vertices.cos_phi[i],
                    vertices.sin_phi[i]);

                if (fast_inside != checked_inside) {
                    std::cerr << std::setprecision(std::numeric_limits<Double>::max_digits10)
                              << "FAST/CHECKED DISAGREEMENT\n"
                              << "  vertex_id=" << vertices.id[i] << '\n'
                              << "  vertex_r=" << vertices.r[i] << '\n'
                              << "  vertex_phi=" << vertices.phi[i] << '\n'
                              << "  center_r=" << center_r_ << '\n'
                              << "  center_phi=" << center_phi_ << '\n'
                              << "  hyperball_radius=" << gen_.current_hyperedge_radius_ << '\n'
                              << "  fast_inside=" << fast_inside << '\n'
                              << "  checked_inside=" << checked_inside << '\n';
                }
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

    std::vector<AnnulusMidpoint>                                annulus_mid_;
    mutable std::unordered_map<SInt, CachedExactCell>           exact_vertices_;
    mutable std::list<SInt>                                     exact_lru_;
    mutable std::unordered_map<SInt, std::list<SInt>::iterator> exact_lru_pos_;
    std::size_t                                                 exact_remote_cache_limit_ = 1024;

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
};

} // namespace kagen
