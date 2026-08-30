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

    explicit HyperbolicGeometryPolicy(GeneratorT& generator, const SInt annulus_id, const SInt chunk_id)
        : gen_(generator),
          rng_(RNGWrapper(generator.config_)) {
        current_annulus_half_angle_.resize(gen_.total_annuli_, Double{-1.0});
        EnsureAnnulusMidpoints();
    }

    void AddCenter(
        const Center& center, std::vector<SInt>& /*pins*/
    ) const {
        const SInt sampling_seed = sampling::Spooky::hash(gen_.config_.seed + (131 * center.sampled_id));

        switch (gen_.config_.partial_cell_mode) {
            case PartialCellMode::EstimateByCoverageRange:
                gen_.mersenne.RandomInit(sampling_seed);
                break;

            case PartialCellMode::GenerateAndCheck:
                break;

            case PartialCellMode::EstimateByCoverageFloyd:
                rng_.SeedUniformStream(sampling_seed);
                break;
        }
    }

    void CandidateCells(const Center& center, const Double radius, std::vector<Cell>& cells) {
        CacheQueryState(center, radius);

        CandidateCollector collector{*this};
        collector.Collect(center, radius, cells);
        if (collector.duplicate_candidates_ != 0) {
            std::cout << "duplicate_candidates_" << collector.duplicate_candidates_
                      << " at hyperedge center: " << center.phi << "," << center.r << "," << center.sampled_id << ","
                      << center.annulus_id << "\n";
        }
    }
    Double Radius(const Center& /*unused*/) const {
        return gen_.current_hyperedge_radius_;
    }

    CellBallRelation ClassifyCell(const Center& /*center*/, const Double /*radius*/, const Cell& cell) const {
        static std::uint64_t calls   = 0;
        static std::uint64_t outside = 0;
        static std::uint64_t inside  = 0;
        static std::uint64_t partial = 0;

        CellBallRelation result;

        if (CellAABBOutsideBall(cell)) {
            result = CellBallRelation::OUTSIDE;
            ++outside;
        } else if (ShouldTryInside(cell) && HyperbolicCellInside(cell)) {
            result = CellBallRelation::INSIDE;
            ++inside;
        } else {
            result = CellBallRelation::PARTIAL;
            ++partial;
        }

        ++calls;

        if (calls % 100000 == 0) {
            std::cerr << "[cell classification]"
                      << " calls=" << calls << " inside=" << static_cast<double>(inside) / static_cast<double>(calls)
                      << " partial=" << static_cast<double>(partial) / static_cast<double>(calls)
                      << " outside=" << static_cast<double>(outside) / static_cast<double>(calls) << '\n';
        }

        return result;
    }

    CellBallRelation ClassifyRegion(const CellRegion& region) const {
        const poincare_geometry::AABB<Double> box{
            .min_x = region.min_x,
            .max_x = region.max_x,
            .min_y = region.min_y,
            .max_y = region.max_y,
        };

        if (poincare_geometry::AABBOutsideBall(box, ball_)) {
            return CellBallRelation::OUTSIDE;
        }

        const Double max_delta = MaxAngularDistanceToInterval(center_phi_, region.min_phi, region.max_phi);

        const Double cos_delta = std::cos(max_delta);

        const Double min_cosh = std::cosh(region.min_r);

        const Double min_sinh = std::sinh(region.min_r);

        const Double max_cosh = std::cosh(region.max_r);

        const Double max_sinh = std::sinh(region.max_r);

        const Double cosh_d_min_r = (center_cosh_r_ * min_cosh) - (center_sinh_r_ * min_sinh * cos_delta);

        if (cosh_d_min_r > radius_cosh_) {
            return CellBallRelation::PARTIAL;
        }

        const Double cosh_d_max_r = (center_cosh_r_ * max_cosh) - (center_sinh_r_ * max_sinh * cos_delta);

        if (cosh_d_max_r <= radius_cosh_) {
            return CellBallRelation::INSIDE;
        }

        return CellBallRelation::PARTIAL;
    }

    CellRegion
    MakeRegion(const SInt first_annulus, const SInt last_annulus, const Double min_phi, const Double max_phi) const {
        if (first_annulus < 0 || last_annulus >= gen_.total_annuli_ || first_annulus > last_annulus) {
            throw std::out_of_range("MakeRegion: invalid annulus range");
        }

        if (!(min_phi < max_phi)) {
            throw std::logic_error("MakeRegion: invalid angular interval");
        }

        const Double min_r = gen_.annulus_min_r_[first_annulus];

        const Double max_r = gen_.annulus_max_r_[last_annulus];

        const auto box = poincare_geometry::ComputeCellAABB(min_r, max_r, min_phi, max_phi, gen_.cell_eps_);

        return CellRegion{
            .first_annulus = first_annulus,
            .last_annulus  = last_annulus,

            .min_r = min_r,
            .max_r = max_r,

            .min_phi = min_phi,
            .max_phi = max_phi,

            .min_x = box.min_x,
            .max_x = box.max_x,
            .min_y = box.min_y,
            .max_y = box.max_y,
        };
    }

    std::pair<CellRegion, CellRegion> SplitRegionRadially(const CellRegion& region) const {
        if (region.first_annulus >= region.last_annulus) {
            throw std::logic_error("SplitRegionRadially: cannot split single-annulus region");
        }

        const SInt mid_annulus = region.first_annulus + (region.last_annulus - region.first_annulus) / 2;

        CellRegion inner = MakeRegion(region.first_annulus, mid_annulus, region.min_phi, region.max_phi);

        CellRegion outer = MakeRegion(mid_annulus + 1, region.last_annulus, region.min_phi, region.max_phi);

        return {std::move(inner), std::move(outer)};
    }

    std::pair<CellRegion, CellRegion> SplitRegionAngularly(const CellRegion& region) const {
        const Double mid_phi = (region.min_phi + region.max_phi) / Double{2.0};

        if (!(region.min_phi < mid_phi && mid_phi < region.max_phi)) {
            throw std::logic_error("SplitRegionAngularly: cannot split angular interval");
        }

        CellRegion left = MakeRegion(region.first_annulus, region.last_annulus, region.min_phi, mid_phi);

        CellRegion right = MakeRegion(region.first_annulus, region.last_annulus, mid_phi, region.max_phi);

        return {std::move(left), std::move(right)};
    }

    void TraverseRegion(
        const CellRegion& region, std::vector<CellRegion>& leaf_regions, std::vector<CellRegion>& inside_regions,
        std::uint64_t& regions_visited, std::uint64_t& leaves_visited, std::uint64_t& outside_regions,
        std::uint64_t& inside_regions_count, std::uint64_t& partial_regions) const {
        const CellBallRelation relation = ClassifyRegion(region);
        switch (relation) {
            case CellBallRelation::OUTSIDE:
                ++outside_regions;
                break;

            case CellBallRelation::INSIDE:
                ++inside_regions_count;
                break;

            case CellBallRelation::PARTIAL:
                ++partial_regions;
                break;
        }

        ++regions_visited;
        if (relation == CellBallRelation::OUTSIDE) {
            return;
        }

        if (relation == CellBallRelation::INSIDE) {
            inside_regions.push_back(region);
            return;
        }

        const bool single_annulus = region.first_annulus == region.last_annulus;

        const Double angular_width = region.max_phi - region.min_phi;

        if (!single_annulus) {
            const auto [inner, outer] = SplitRegionRadially(region);

            TraverseRegion(
                inner, leaf_regions, inside_regions, regions_visited, leaves_visited, outside_regions,
                inside_regions_count, partial_regions);

            TraverseRegion(
                outer, leaf_regions, inside_regions, regions_visited, leaves_visited, outside_regions,
                inside_regions_count, partial_regions);

            return;
        }

        // From here on, the region is exactly one annulus.
        const Double leaf_width = gen_.CellWidthForChunkAnnulus(region.first_annulus, 0);

        if (angular_width <= leaf_width) {
            leaf_regions.push_back(region);
            ++leaves_visited;
            return;
        }

        const auto [left, right] = SplitRegionAngularly(region);

        TraverseRegion(
            left, leaf_regions, inside_regions, regions_visited, leaves_visited, outside_regions, inside_regions_count,
            partial_regions);

        TraverseRegion(
            right, leaf_regions, inside_regions, regions_visited, leaves_visited, outside_regions, inside_regions_count,
            partial_regions);
    }

    Double CellCoverage(const Center& center, const Double /*hyperball_radius*/, const Cell& cell) const {
        const Double cell_phi_width = cell.max_phi - cell.min_phi;

        if (cell_phi_width <= Double{0.0}) {
            return Double{0.0};
        }

        // 8-point Gauss-Legendre rule on [-1, 1].
        static constexpr double nodes[] = {
            -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
            0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363,
        };

        static constexpr double weights[] = {
            0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
            0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763,
        };

        const Double u_min = std::cosh(gen_.alpha_ * cell.min_r);
        const Double u_max = std::cosh(gen_.alpha_ * cell.max_r);

        const Double u_mid  = (u_min + u_max) / Double{2.0};
        const Double u_half = (u_max - u_min) / Double{2.0};

        Double weighted_overlap = Double{0.0};

        for (std::size_t i = 0; i < 8; ++i) {
            const Double u = u_mid + (u_half * static_cast<Double>(nodes[i]));

            const Double r = std::acosh(u) / gen_.alpha_;

            const Double cosh_r = std::cosh(r);
            const Double sinh_r = std::sinh(r);

            const Double half_angle = AllowedHalfAngleAtRadius(r, cosh_r, sinh_r);

            Double overlap = Double{0.0};

            if (half_angle >= Double{M_PI}) {
                overlap = cell_phi_width;
            } else if (half_angle > Double{0.0}) {
                overlap = circular_interval::CircularOverlap(
                    cell.min_phi, cell.max_phi, center.phi - half_angle, center.phi + half_angle, gen_.cell_eps_);
            }

            weighted_overlap += static_cast<Double>(weights[i]) * overlap;
        }

        // Division by 2 is the normalization of Gauss-Legendre weights
        // when computing the average over [u_min, u_max].
        const Double mean_overlap = weighted_overlap / Double{2.0};

        return std::clamp(mean_overlap / cell_phi_width, Double{0.0}, Double{1.0});
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
        const Center& /*center*/, const Cell& cell, const Double coverage, std::vector<SInt>& /*pins*/,
        std::vector<PinRange>& ranges) const {
        const auto info = GetPartialCellSampleInfo(cell, coverage);
        if (!info) {
            return 0;
        }

        auto sampled = getRandomPinRange(info->size, info->k, info->offset, gen_.mersenne);
        ranges.push_back(sampled);

        return info->k;
    }
    SInt AddPartialCellFloyd(
        const Center& /*center*/, const Cell& cell, const Double coverage, std::vector<SInt>& pins,
        std::vector<PinRange>& /*ranges*/) const {
        const auto info = GetPartialCellSampleInfo(cell, coverage);
        if (!info) {
            return 0;
        }

        FloydSampleGeometricAppend(info->offset, info->size, info->k, rng_, pins, floyd_scratch_);
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

    bool ShouldApproximatePartialCell(const Cell& cell) const {
        return !IsLocalCell(cell) && !exact_vertices_.contains(cell.global_cell_id);
    }

private:
    // ===== Query state =====
    void CacheQueryState(const Center& center, const Double radius) {
        center_cosh_r_ = std::cosh(center.r);
        center_sinh_r_ = std::sinh(center.r);
        radius_cosh_   = std::cosh(radius);

        center_cos_phi_ = std::cos(center.phi);
        center_sin_phi_ = std::sin(center.phi);

        center_phi_ = center.phi;
        center_r_   = center.r;

        const Double center_inv_len = (center_cosh_r_ + Double{1.0}) / Double{2.0};

        const Double center_pdm_radius = std::sqrt(Double{1.0} - (Double{1.0} / center_inv_len));

        center_x_ = center_pdm_radius * center_sin_phi_;
        center_y_ = center_pdm_radius * center_cos_phi_;

        center_vertex_ = Vertex{
            center.phi,
            center.r,
            center_pdm_radius * center_sin_phi_,
            center_pdm_radius * center_cos_phi_,
            center_inv_len,
            SInt{0},
            center_cosh_r_,
            center_sinh_r_,
            center_cos_phi_,
            center_sin_phi_};

        center_gamma_ = center_inv_len;

        ball_ = poincare_geometry::MakeBall(center.r, center.phi, radius);
    }

    // ===== Candidate search =====
    struct CandidateCollector {
        HyperbolicGeometryPolicy& policy;
        std::unordered_set<SInt>  seen_candidate_cells_;
        GeneratorT&               gen() {
            return policy.gen_;
        }

        void Collect(const Center& center, const Double radius, std::vector<Cell>& cells) {
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
            const auto& occupied = gen().nonempty_cells_per_annulus_[annulus_id];

            if (occupied.empty()) {
                return;
            }

            const auto parts = circular_interval::Split(min_phi, max_phi, gen().cell_eps_);

            for (int p = 0; p < parts.count; ++p) {
                const Double q_begin = parts.parts[p].first;

                const Double q_end = parts.parts[p].second;

                if (!(q_begin < q_end)) {
                    continue;
                }

                const auto [first, last] = gen().GlobalCellRangeForAngularInterval(annulus_id, q_begin, q_end);

                auto begin = std::lower_bound(occupied.begin(), occupied.end(), first);

                auto end = std::upper_bound(begin, occupied.end(), last);

                for (auto it = begin; it != end; ++it) {
                    PushGlobalCell(annulus_id, *it, cells);
                }
            }
        }

        void AddCandidateCellsInAnnulus(const Center& center, const SInt annulus_id, std::vector<Cell>& cells) {
            const Double half_angle = policy.AllowedHalfAngleForAnnulus(center.r, annulus_id);

            policy.current_annulus_half_angle_[annulus_id] = half_angle;

            if (!(half_angle > Double{0.0})) {
                return;
            }

            const auto& occupied = gen().nonempty_cells_per_annulus_[annulus_id];

            if (occupied.empty()) {
                return;
            }

            //
            // Entire angular extent of this annulus.
            //
            if (half_angle >= Double{M_PI}) {
                for (const SInt global_cell: occupied) {
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

            const auto it = gen().cells_.find(global_cell_id);

            if (it == gen().cells_.end()) {
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
            const SInt global_chunk_id = gen().ComputeGlobalChunkId(annulus_id, chunk_id);
            const SInt global_cell_id  = gen().ComputeGlobalCellId(annulus_id, chunk_id, cell_id);

            const auto& annulus     = gen().annuli_.find(global_chunk_id)->second;
            const auto& stored_cell = gen().cells_.find(global_cell_id)->second;

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
    Double AllowedHalfAngleForAnnulus(const Double /*center_r*/, const SInt annulus_id) const {
        const Double min_r = gen_.annulus_min_r_[annulus_id];
        const Double max_r = gen_.annulus_max_r_[annulus_id];

        Double reach = std::max(
            AllowedHalfAngleForCachedRadius(gen_.annulus_min_cosh_[annulus_id], gen_.annulus_min_sinh_[annulus_id]),
            AllowedHalfAngleForCachedRadius(gen_.annulus_max_cosh_[annulus_id], gen_.annulus_max_sinh_[annulus_id]));

        // The angular reach is not generally maximized at center_r.  Its only
        // interior stationary point satisfies cosh(query_r) = cosh(center_r) / cosh(radius).
        const Double critical_cosh = center_cosh_r_ / radius_cosh_;
        if (critical_cosh >= Double{1.0}) {
            const Double critical_r = std::acosh(critical_cosh);
            if (critical_r >= min_r && critical_r <= max_r) {
                reach = std::max(reach, AllowedHalfAngleForCachedRadius(critical_cosh, std::sinh(critical_r)));
            }
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

        const SInt a = cell.annulus_id;

        const Double cosh_d_min_r =
            (center_cosh_r_ * gen_.annulus_min_cosh_[a]) - (center_sinh_r_ * gen_.annulus_min_sinh_[a] * cos_delta);

        if (cosh_d_min_r > radius_cosh_) {
            return false;
        }

        const Double cosh_d_max_r =
            (center_cosh_r_ * gen_.annulus_max_cosh_[a]) - (center_sinh_r_ * gen_.annulus_max_sinh_[a] * cos_delta);

        return cosh_d_max_r <= radius_cosh_;
    }
    bool ShouldTryInside(const Cell& cell) const {
        return gen_.current_hyperedge_radius_ >= (cell.max_r - cell.min_r);
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
    };
    std::optional<PartialCellSampleInfo> GetPartialCellSampleInfo(const Cell& cell, const Double coverage) const {
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

        return PartialCellSampleInfo{
            .global_cell_id = global_cell_id,
            .size           = size,
            .offset         = offset,
            .k              = k,
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
    bool HierarchicalCandidateCells(
        const Center& center, LPFloat radius, std::vector<Cell>& partial_cells, std::vector<PinRange>& ranges) const {
        return false;
    }

private:
    bool IsInsideHyperbolicBallFast(const Vertex& vertex) const {
        return PGGeometry<Double>::HyperbolicDistance(center_vertex_, vertex) <= gen_.current_hyperedge_pdm_radius_;
    }

    SInt AddExactVerticesFast(const typename GeneratorT::VertexBlock& vertices, std::vector<SInt>& pins) const {
        SInt accepted = 0;

        for (std::size_t i = 0; i < vertices.id.size(); ++i) {
            const Vertex vertex{vertices.phi[i],     vertices.r[i],      vertices.x[i],      vertices.y[i],
                                vertices.gamma[i],   vertices.id[i],     vertices.cosh_r[i], vertices.sinh_r[i],
                                vertices.cos_phi[i], vertices.sin_phi[i]};

            if (IsInsideHyperbolicBallFast(vertex)) {
                pins.push_back(vertices.id[i]);
                ++accepted;
            }
        }

        return accepted;
    }

    SInt AddExactVerticesChecked(const typename GeneratorT::VertexBlock& vertices, std::vector<SInt>& pins) const {
        return AddExactVerticesFast(vertices, pins);
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