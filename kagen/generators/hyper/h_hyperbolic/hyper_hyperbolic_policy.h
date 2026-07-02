#pragma once

#include "kagen/generators/hyper/h_hyperbolic/hyper_hyperbolic.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
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

    struct PoincareBall {
        Double cx;
        Double cy;
        Double radius;
        Double radius_sq;
    };

    struct AABB {
        Double min_x;
        Double max_x;
        Double min_y;
        Double max_y;
    };
    Double NormalizePhi(Double phi) const {
        const Double two_pi = 2.0 * M_PI;

        while (phi < 0.0) {
            phi += two_pi;
        }

        while (phi >= two_pi) {
            phi -= two_pi;
        }

        return phi;
    }

    bool AngleInInterval(Double phi, Double min_phi, Double max_phi) const {
        phi     = NormalizePhi(phi);
        min_phi = NormalizePhi(min_phi);
        max_phi = NormalizePhi(max_phi);

        if (min_phi <= max_phi) {
            return min_phi <= phi && phi <= max_phi;
        }

        // interval wraps around 2π, e.g. [5.9, 0.2]
        return phi >= min_phi || phi <= max_phi;
    }

    AABB
    ComputePoincareCellAABB(const Double min_r, const Double max_r, const Double min_phi, const Double max_phi) const {
        const bool   full_circle = std::abs((max_phi - min_phi) - Double{2.0 * M_PI}) <= gen_.cell_eps_;
        const Double rho_min     = std::tanh(min_r / Double{2.0});
        const Double rho_max     = std::tanh(max_r / Double{2.0});

        std::vector<Double> phis = {min_phi, max_phi};

        const Double critical_angles[] = {
            Double{0.0},
            Double{M_PI / 2.0},
            Double{M_PI},
            Double{3.0 * M_PI / 2.0},
        };

        for (const Double angle: critical_angles) {
            if (full_circle || AngleInInterval(angle, min_phi, max_phi)) {
                phis.push_back(angle);
            }
        }

        Double min_x = std::numeric_limits<Double>::infinity();
        Double max_x = -std::numeric_limits<Double>::infinity();
        Double min_y = std::numeric_limits<Double>::infinity();
        Double max_y = -std::numeric_limits<Double>::infinity();

        for (const Double phi: phis) {
            for (const Double rho: {rho_min, rho_max}) {
                const Double x = rho * std::sin(phi);
                const Double y = rho * std::cos(phi);

                min_x = std::min(min_x, x);
                max_x = std::max(max_x, x);
                min_y = std::min(min_y, y);
                max_y = std::max(max_y, y);
            }
        }

        return {.min_x = min_x, .max_x = max_x, .min_y = min_y, .max_y = max_y};
    }

    explicit HyperbolicGeometryPolicy(GeneratorT& generator, const SInt annulus_id, const SInt chunk_id)
        : gen_(generator),
          start_annulus_id_(annulus_id),
          start_chunk_id_(chunk_id),
          rng_(RNGWrapper(generator.config_)) {
        current_annulus_half_angle_.resize(gen_.total_annuli_, Double{-1.0});
        EnsureAnnulusMidpoints();
    }

    void SetStartCell(const SInt annulus_id, const SInt chunk_id, const SInt cell_id) {
        start_annulus_id_ = annulus_id;
        start_chunk_id_   = chunk_id;
        start_cell_id_    = cell_id;
    }

    void AddCenter(const Center& /*unused*/, std::vector<SInt>& /*unused*/) const {}

    Double Radius(const Center& /*unused*/) const {
        return gen_.current_hyperedge_radius_;
    }

    void AddCellsInAngularInterval(
        const SInt annulus_id, const Double min_phi, const Double max_phi, std::vector<Cell>& cells) {
        const SInt grid_size   = gen_.GridSizeForAnnulus(annulus_id);
        const auto query_parts = SplitCircularInterval(min_phi, max_phi);

        for (SInt chunk_id = 0; chunk_id < gen_.config_.k; ++chunk_id) {
            const auto& chunk = gen_.chunks_[chunk_id];

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

            gen_.GenerateCells(annulus_id, chunk_id);

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
    void AddCandidateCellsInAnnulus(
        const Center& center, const Double /*radius*/, const SInt annulus_id, std::vector<Cell>& cells) {
        const Double half_angle                 = AllowedHalfAngleForAnnulus(center.r, annulus_id);
        current_annulus_half_angle_[annulus_id] = half_angle;

        if (half_angle <= 0.0) {
            return;
        }

        AddCellsInAngularInterval(annulus_id, center.phi - half_angle, center.phi + half_angle, cells);
    }

    void CandidateCells(const Center& center, const Double radius, std::vector<Cell>& cells) {
        center_cosh_r_  = std::cosh(center.r);
        center_sinh_r_  = std::sinh(center.r);
        radius_cosh_    = std::cosh(radius);
        ball_           = MakePoincareBall(center, radius);
        center_cos_phi_ = std::cos(center.phi);
        center_sin_phi_ = std::sin(center.phi);
        center_phi_     = center.phi;
        center_r_       = center.r;

        cells.clear();

        const Double r_min = std::max(Double{0.0}, center.r - radius);
        const Double r_max = std::min(gen_.target_r_, center.r + radius);

        if (r_max < r_min) {
            return;
        }

        const Double annulus_width = gen_.target_r_ / static_cast<Double>(gen_.total_annuli_);

        SInt first_annulus = static_cast<SInt>(std::floor(r_min / annulus_width));
        SInt last_annulus  = static_cast<SInt>(std::floor(r_max / annulus_width));

        first_annulus = std::clamp<SInt>(first_annulus, 0, gen_.total_annuli_ - 1);
        last_annulus  = std::clamp<SInt>(last_annulus, 0, gen_.total_annuli_ - 1);

        std::fill(current_annulus_half_angle_.begin(), current_annulus_half_angle_.end(), Double{-1.0});

        for (SInt annulus_id = first_annulus; annulus_id <= last_annulus; ++annulus_id) {
            AddCandidateCellsInAnnulus(center, radius, annulus_id, cells);
        }
    }

    Double AngularDistance(Double a, Double b) const {
        a = NormalizePhi(a);
        b = NormalizePhi(b);

        Double d = std::abs(a - b);
        return std::min(d, Double{2.0 * M_PI} - d);
    }

    Double CoshDistanceWithDelta(Double query_r, Double delta_phi) const {
        return (center_cosh_r_ * std::cosh(query_r)) - (center_sinh_r_ * std::sinh(query_r) * std::cos(delta_phi));
    }

    Double MinAngularDistanceToInterval(Double phi, Double min_phi, Double max_phi) const {
        if (AngleInInterval(phi, min_phi, max_phi)) {
            return 0.0;
        }

        return std::min(AngularDistance(phi, min_phi), AngularDistance(phi, max_phi));
    }

    bool CellsOutside(
        const Center& center, const Double /*hyperball_radius*/, const Double min_r, const Double max_r,
        const Double min_phi, const Double max_phi) const {
        const Double min_delta = MinAngularDistanceToInterval(center.phi, min_phi, max_phi);

        Double       best_r    = center.r;
        const Double candidate = center_sinh_r_ * std::cos(min_delta) / center_cosh_r_;

        if (candidate > 0.0 && candidate < 1.0) {
            best_r = std::atanh(candidate);
        }

        best_r = std::clamp(best_r, min_r, max_r);

        const Double min_cosh_dist = CoshDistanceWithDelta(best_r, min_delta);
        return min_cosh_dist > radius_cosh_;
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

    CellBallRelation ClassifyCell(const Center& /*center*/, const Double /*radius*/, const Cell& cell) const {
        if (PoincareCellOutside(cell)) {
            return CellBallRelation::OUTSIDE;
        }

        if (ShouldTryInside(cell) && HyperbolicCellInside(cell)) {
            return CellBallRelation::INSIDE;
        }

        return CellBallRelation::PARTIAL;
    }

    Double IntervalOverlap(Double a_begin, Double a_end, Double b_begin, Double b_end) const {
        const Double begin = std::max(a_begin, b_begin);
        const Double end   = std::min(a_end, b_end);

        return std::max<Double>(0.0, end - begin);
    }

    struct SplitInterval {
        std::pair<Double, Double> parts[2];
        int                       count;
    };

    SplitInterval SplitCircularInterval(Double begin, Double end) const {
        const Double two_pi = Double{2.0 * M_PI};
        const Double width  = end - begin;

        SplitInterval result{};
        result.parts[0] = {Double{0.0}, Double{0.0}};
        result.parts[1] = {Double{0.0}, Double{0.0}};
        result.count    = 0;

        if (std::abs(std::abs(width) - two_pi) <= gen_.cell_eps_ || std::abs(width) > two_pi) {
            result.parts[0] = {Double{0.0}, two_pi};
            result.count    = 1;
            return result;
        }

        begin = NormalizePhi(begin);
        end   = NormalizePhi(end);

        if (begin <= end) {
            result.parts[0] = {begin, end};
            result.count    = 1;
            return result;
        }

        result.parts[0] = {begin, two_pi};
        result.parts[1] = {Double{0.0}, end};
        result.count    = 2;
        return result;
    }

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

    Double CircularIntervalOverlap(Double a_begin, Double a_end, Double b_begin, Double b_end) const {
        Double overlap = Double{0.0};

        const auto a_parts = SplitCircularInterval(a_begin, a_end);
        const auto b_parts = SplitCircularInterval(b_begin, b_end);

        for (int i = 0; i < a_parts.count; ++i) {
            for (int j = 0; j < b_parts.count; ++j) {
                overlap += IntervalOverlap(
                    a_parts.parts[i].first, a_parts.parts[i].second, b_parts.parts[j].first, b_parts.parts[j].second);
            }
        }

        return overlap;
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
            overlap_phi = CircularIntervalOverlap(min_phi, max_phi, center.phi - half_angle, center.phi + half_angle);
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
        const SInt global_cell_id = cell.global_cell_id;

        const auto it = gen_.cells_.find(global_cell_id);
        if (it == gen_.cells_.end()) {
            return 0;
        }

        const auto& stored_cell = it->second;

        const SInt size   = std::get<0>(stored_cell);
        const SInt offset = std::get<4>(stored_cell);
        const SInt k      = static_cast<SInt>(std::floor(static_cast<Double>(size) * coverage));

        if (k <= 0) {
            return 0;
        }

        SInt seed = sampling::Spooky::hash(gen_.config_.seed + (131 * center.sampled_id) + (9973 * global_cell_id));

        auto sampled = getRandomPinRange(size, k, offset, seed, gen_.mersenne);

        ranges.insert(ranges.end(), sampled);

        return k;
    }

    SInt AddPartialCellFloyd(
        const Center& center, const Cell& cell, const Double coverage, std::vector<SInt>& pins,
        std::vector<PinRange>& ranges) const {
        const SInt global_cell_id = cell.global_cell_id;

        const auto it = gen_.cells_.find(global_cell_id);
        if (it == gen_.cells_.end()) {
            return 0;
        }

        const auto& stored_cell = it->second;

        const SInt size   = std::get<0>(stored_cell);
        const SInt offset = std::get<4>(stored_cell);
        const SInt k      = static_cast<SInt>(std::floor(static_cast<Double>(size) * coverage));

        if (k <= 0) {
            return 0;
        }

        SInt seed = sampling::Spooky::hash(gen_.config_.seed + (131 * center.sampled_id) + (9973 * global_cell_id));

        FloydSampleGeometricAppend(offset, size, k, rng_, seed, pins, floyd_scratch_);

        return k;
    }

    SInt AddPartialCellExact(const Center& center, const Double /*radius*/, const Cell& cell, std::vector<SInt>& pins) {
        auto it = gen_.vertices_.find(cell.global_cell_id);

        if (it == gen_.vertices_.end()) {
            gen_.GenerateVertices(cell.annulus_id, cell.chunk_id, cell.cell_id);
            it = gen_.vertices_.find(cell.global_cell_id);
        }

        if (it == gen_.vertices_.end()) {
            return 0;
        }

        const auto&       vertices       = it->second;
        SInt              vertex_counter = 0;
        const std::size_t n              = vertices.id.size();

        const Double* __restrict cosh_r  = vertices.cosh_r.data();
        const Double* __restrict sinh_r  = vertices.sinh_r.data();
        const Double* __restrict cos_phi = vertices.cos_phi.data();
        const Double* __restrict sin_phi = vertices.sin_phi.data();
        const SInt* __restrict ids       = vertices.id.data();
        const Double* __restrict phi     = vertices.phi.data();

        if (gen_.config_.debug) {
            for (std::size_t i = 0; i < n; ++i) {
                if (IsInsideHyperbolicBallChecked(phi[i], cosh_r[i], sinh_r[i], cos_phi[i], sin_phi[i])) {
                    pins.push_back(ids[i]);
                    ++vertex_counter;
                }
            }
        } else {
            for (std::size_t i = 0; i < n; ++i) {
                if (IsInsideHyperbolicBallFast(cosh_r[i], sinh_r[i], cos_phi[i], sin_phi[i])) {
                    pins.push_back(ids[i]);
                    ++vertex_counter;
                }
            }
        }
        return vertex_counter;
    }

    void EmitHyperedge(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges) {
        gen_.PushHyperedgeCompressed(pins, ranges);
    }

    std::string CenterToString(const Center& center) const {
        std::ostringstream out;

        out << "phi=" << center.phi << ";r=" << center.r;

        return out.str();
    }

    PoincareBall MakePoincareBall(const Center& center, const Double rho) const {
        const Double s  = std::tanh(center.r / Double{2.0});
        const Double x0 = s * std::sin(center.phi);
        const Double y0 = s * std::cos(center.phi);

        const Double t     = std::tanh(rho / Double{2.0});
        const Double denom = Double{1.0} - (s * s * t * t);

        const Double scale            = (Double{1.0} - (t * t)) / denom;
        const Double euclidean_radius = ((Double{1.0} - (s * s)) * t) / denom;

        const Double cx = scale * x0;
        const Double cy = scale * y0;

        return {
            .cx        = cx,
            .cy        = cy,
            .radius    = euclidean_radius,
            .radius_sq = euclidean_radius * euclidean_radius,
        };
    }

private:
    GeneratorT&          gen_;
    SInt                 start_annulus_id_;
    SInt                 start_chunk_id_;
    SInt                 start_cell_id_{};
    mutable RNGWrapper<> rng_;

    // Caches
    Double                           center_cosh_r_{};
    Double                           center_sinh_r_{};
    Double                           radius_cosh_{};
    PoincareBall                     ball_{};
    Double                           center_cos_phi_{};
    Double                           center_sin_phi_{};
    Double                           center_phi_{};
    Double                           center_r_{};
    mutable std::unordered_set<SInt> floyd_scratch_;
    mutable std::vector<Double>      current_annulus_half_angle_;
    struct CachedExactCell {
        std::vector<Vertex> vertices_by_x;
    };

    mutable std::unordered_map<SInt, CachedExactCell> exact_vertices_;

    struct CoverageSample {
        Double r;
        Double cosh_r;
        Double sinh_r;
        Double weight;
    };
    static constexpr int                                              kCoverageSamples = 4;
    mutable std::vector<std::array<CoverageSample, kCoverageSamples>> coverage_samples_;

    struct AnnulusMidpoint {
        Double r;
        Double cosh_r;
        Double sinh_r;
    };

    std::vector<AnnulusMidpoint> annulus_mid_;

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

    void EnsureCoverageSamples() const {
        if (!coverage_samples_.empty()) {
            return;
        }

        coverage_samples_.resize(gen_.total_annuli_);

        const Double alpha = (gen_.config_.plexp - 1.0) / 2.0;

        for (SInt annulus_id = 0; annulus_id < gen_.total_annuli_; ++annulus_id) {
            const Double min_r = gen_.annulus_min_r_[annulus_id];
            const Double max_r = gen_.annulus_max_r_[annulus_id];
            const Double dr    = (max_r - min_r) / kCoverageSamples;

            for (int i = 0; i < kCoverageSamples; ++i) {
                const Double r = min_r + ((static_cast<Double>(i) + 0.5) * dr);

                coverage_samples_[annulus_id][i] = {
                    .r      = r,
                    .cosh_r = std::cosh(r),
                    .sinh_r = std::sinh(r),
                    .weight = std::sinh(alpha * r),
                };
            }
        }
    }

    void
    PushCandidateCell(const SInt annulus_id, const SInt chunk_id, const SInt cell_id, std::vector<Cell>& cells) const {
        const SInt global_chunk_id = gen_.ComputeGlobalChunkId(annulus_id, chunk_id);
        const SInt global_cell_id  = gen_.ComputeGlobalCellId(annulus_id, chunk_id, cell_id);

        const auto& annulus     = gen_.annuli_.find(global_chunk_id)->second;
        const auto& stored_cell = gen_.cells_.find(global_cell_id)->second;

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

    void QueryCollect(
        const SInt annulus_id, const SInt chunk_id, const SInt cell_id, const Center& center, const bool search_down,
        std::vector<Cell>& cells) {
        auto& annulus = gen_.annuli_[gen_.ComputeGlobalChunkId(annulus_id, chunk_id)];
        auto  bounds  = GetBoundaryPhis(center.phi, center.r, annulus_id);

        gen_.current_min_phi_ = bounds.first;
        gen_.current_max_phi_ = bounds.second;

        const auto& stored_cell = gen_.cells_[gen_.ComputeGlobalCellId(annulus_id, chunk_id, cell_id)];

        const Double min_cell_phi = std::get<1>(stored_cell);
        const Double max_cell_phi = std::get<2>(stored_cell);

        gen_.right_processed_chunk_ = chunk_id;
        gen_.right_processed_cell_  = cell_id;

        if (search_down || !gen_.IsLocalChunk(chunk_id)) {
            PushCandidateCell(annulus_id, chunk_id, cell_id, cells);
        }

        if (std::get<1>(annulus) >= gen_.clique_thres_
            && std::max(gen_.TotalGridSizeForAnnulus(annulus_id), gen_.config_.k) > 1) {
            if (gen_.current_min_phi_ < min_cell_phi
                || (gen_.OutOfBounds(gen_.current_min_phi_) && !(std::abs(min_cell_phi - 0.0) < gen_.cell_eps_))) {
                SInt next_chunk_id = chunk_id;

                if (cell_id == 0) {
                    next_chunk_id = (chunk_id + gen_.config_.k - 1) % gen_.config_.k;
                }

                const SInt next_cell_id =
                    (cell_id + gen_.GridSizeForAnnulus(annulus_id) - 1) % gen_.GridSizeForAnnulus(annulus_id);

                gen_.GenerateCells(annulus_id, next_chunk_id);

                QueryRightNeighborCollect(
                    annulus_id, next_chunk_id, next_cell_id, center, std::abs(min_cell_phi - 0.0) < gen_.cell_eps_,
                    search_down, cells);
            }

            if (gen_.current_max_phi_ > max_cell_phi
                || (gen_.OutOfBounds(gen_.current_max_phi_)
                    && !(std::abs(max_cell_phi - (2 * M_PI)) < gen_.cell_eps_))) {
                SInt next_chunk_id = chunk_id;

                if (cell_id == gen_.GridSizeForAnnulus(annulus_id) - 1) {
                    next_chunk_id = (chunk_id + gen_.config_.k + 1) % gen_.config_.k;
                }

                const SInt next_cell_id =
                    (cell_id + gen_.GridSizeForAnnulus(annulus_id) + 1) % gen_.GridSizeForAnnulus(annulus_id);

                gen_.GenerateCells(annulus_id, next_chunk_id);

                QueryLeftNeighborCollect(
                    annulus_id, next_chunk_id, next_cell_id, center,
                    std::abs(max_cell_phi - (2 * M_PI)) < gen_.cell_eps_, search_down, cells);
            }
        }

        const SSInt next_annulus =
            search_down ? static_cast<SSInt>(annulus_id + 1) : static_cast<SSInt>(annulus_id - 1);

        if (next_annulus < 0 || next_annulus >= static_cast<SSInt>(gen_.total_annuli_)) {
            return;
        }

        auto& chunk = gen_.chunks_[chunk_id];

        const Double min_chunk_phi = std::get<1>(chunk);
        const Double max_chunk_phi = std::get<2>(chunk);
        const Double grid_phi      = (max_chunk_phi - min_chunk_phi) / gen_.GridSizeForAnnulus(next_annulus);

        const SInt next_cell_id = static_cast<SInt>(std::floor((center.phi - min_chunk_phi) / grid_phi));

        QueryCollect(next_annulus, chunk_id, next_cell_id, center, search_down, cells);
    }

    void QueryRightNeighborCollect(
        const SInt annulus_id, SInt chunk_id, SInt cell_id, const Center& /*center*/, bool phase,
        const bool search_down, std::vector<Cell>& cells) {
        while (true) {
            if (phase && gen_.current_min_phi_ < 0.0) {
                gen_.current_min_phi_ += 2 * M_PI;
            }

            if (phase && gen_.OutOfBounds(gen_.current_min_phi_)) {
                return;
            }

            auto& stored_cell = gen_.cells_[gen_.ComputeGlobalCellId(annulus_id, chunk_id, cell_id)];

            const Double min_cell_phi = std::get<1>(stored_cell);

            gen_.right_processed_chunk_ = chunk_id;
            gen_.right_processed_cell_  = cell_id;

            if (search_down || !gen_.IsLocalChunk(chunk_id)) {
                PushCandidateCell(annulus_id, chunk_id, cell_id, cells);
            }

            phase = phase || std::abs(min_cell_phi - 0.0) < gen_.cell_eps_;

            if (gen_.current_min_phi_ < min_cell_phi || gen_.OutOfBounds(gen_.current_min_phi_)) {
                SInt next_chunk_id = chunk_id;

                if (cell_id == 0) {
                    next_chunk_id = (chunk_id + gen_.config_.k - 1) % gen_.config_.k;
                }

                const SInt next_cell_id =
                    (cell_id + gen_.GridSizeForAnnulus(annulus_id) - 1) % gen_.GridSizeForAnnulus(annulus_id);

                gen_.GenerateCells(annulus_id, next_chunk_id);

                cell_id  = next_cell_id;
                chunk_id = next_chunk_id;
                continue;
            }

            return;
        }
    }

    void QueryLeftNeighborCollect(
        const SInt annulus_id, SInt chunk_id, SInt cell_id, const Center& /*center*/, bool phase,
        const bool search_down, std::vector<Cell>& cells) {
        while (true) {
            if (phase && gen_.current_max_phi_ >= 2 * M_PI) {
                gen_.current_max_phi_ -= 2 * M_PI;
            }

            if (phase && gen_.OutOfBounds(gen_.current_max_phi_)) {
                return;
            }

            if (chunk_id == gen_.right_processed_chunk_ && cell_id == gen_.right_processed_cell_) {
                return;
            }

            auto& stored_cell = gen_.cells_[gen_.ComputeGlobalCellId(annulus_id, chunk_id, cell_id)];

            const Double max_cell_phi = std::get<2>(stored_cell);

            if (search_down || !gen_.IsLocalChunk(chunk_id)) {
                PushCandidateCell(annulus_id, chunk_id, cell_id, cells);
            }

            phase = phase || std::abs(max_cell_phi - (2 * M_PI)) < gen_.cell_eps_;

            if (gen_.current_max_phi_ > max_cell_phi || gen_.OutOfBounds(gen_.current_max_phi_)) {
                SInt next_chunk_id = chunk_id;

                if (cell_id == gen_.GridSizeForAnnulus(annulus_id) - 1) {
                    next_chunk_id = (chunk_id + gen_.config_.k + 1) % gen_.config_.k;
                }

                const SInt next_cell_id =
                    (cell_id + gen_.GridSizeForAnnulus(annulus_id) + 1) % gen_.GridSizeForAnnulus(annulus_id);

                gen_.GenerateCells(annulus_id, next_chunk_id);

                cell_id  = next_cell_id;
                chunk_id = next_chunk_id;
                continue;
            }

            return;
        }
    }

    Double MaxAngularDistanceToInterval(const Double phi, const Double min_phi, const Double max_phi) const {
        const Double antipode = NormalizePhi(phi + Double{M_PI});

        if (AngleInInterval(antipode, min_phi, max_phi)) {
            return Double{M_PI};
        }

        return std::max(AngularDistance(phi, min_phi), AngularDistance(phi, max_phi));
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

    std::pair<Double, Double>
    GetBoundaryPhis(const Double boundary_phi, const Double /*boundary_r*/, const SInt annulus_id) const {
        const auto& boundary = gen_.boundaries_[annulus_id];

        const Double cosh_min_r      = std::get<0>(boundary);
        const Double sinh_min_r      = std::get<1>(boundary);
        const Double cosh_boundary_r = center_cosh_r_;
        const Double sinh_boundary_r = center_sinh_r_;

        if (std::abs(sinh_boundary_r * sinh_min_r) <= std::numeric_limits<Double>::epsilon()) {
            return {boundary_phi - M_PI, boundary_phi + M_PI};
        }

        Double arg = ((cosh_boundary_r * cosh_min_r) - radius_cosh_) / (sinh_boundary_r * sinh_min_r);
        arg        = std::clamp(arg, Double{-1.0}, Double{1.0});

        const Double diff = std::acos(arg);

        return {boundary_phi - diff, boundary_phi + diff};
    }

    bool PoincareCellOutside(const Cell& cell) const {
        const Double closest_x = std::clamp(ball_.cx, cell.min_x, cell.max_x);
        const Double closest_y = std::clamp(ball_.cy, cell.min_y, cell.max_y);

        const Double dx = ball_.cx - closest_x;
        const Double dy = ball_.cy - closest_y;

        return (dx * dx) + (dy * dy) > ball_.radius_sq;
    }

    Double HyperbolicDistanceReference(const Double phi, const Double cosh_r, const Double sinh_r) const {
        Double delta_phi = std::abs(center_phi_ - phi);
        delta_phi        = std::min(delta_phi, Double{2.0 * M_PI} - delta_phi);

        const Double cosh_dist = (center_cosh_r_ * cosh_r) - (center_sinh_r_ * sinh_r * std::cos(delta_phi));

        return (cosh_dist - Double{1.0}) / Double{2.0};
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

    bool IsInsideHyperbolicBallFast(
        const Double v_cosh_r, const Double v_sinh_r, const Double v_cos_phi, const Double v_sin_phi) const {
        const Double term_a       = center_cosh_r_ * v_cosh_r;
        const Double term_b_scale = center_sinh_r_ * v_sinh_r;
        const Double cos_delta    = (center_cos_phi_ * v_cos_phi) + (center_sin_phi_ * v_sin_phi);

        const Double pdm_fast = (term_a - (term_b_scale * cos_delta) - Double{1.0}) / Double{2.0};

        return pdm_fast <= gen_.current_hyperedge_pdm_radius_;
    }
};

} // namespace kagen
