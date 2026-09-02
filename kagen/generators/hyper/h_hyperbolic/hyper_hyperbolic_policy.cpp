
#include "kagen/generators/hyper/h_hyperbolic/hyper_hyperbolic_policy.h"

#include "kagen/generators/hyper/h_hyperbolic/hyper_hyperbolic.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"

#include <unistd.h>

namespace kagen {

template <typename Double>
void HyperbolicGeometryPolicy<Double>::AddCenter(
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

template <typename Double>
bool HyperbolicGeometryPolicy<Double>::HierarchicalCandidateCells(
    const Center& center, const Double radius, std::vector<Cell>& cells, std::vector<PinRange>& ranges) {
    CacheQueryState(center, radius);

    cells.clear();

    const std::size_t ranges_before = ranges.size();

    CandidateCollector collector{*this};
    collector.CollectHierarchical(center, radius, cells, ranges);

    return ranges.size() > ranges_before;
}

template <typename Double>
void HyperbolicGeometryPolicy<Double>::CandidateCells(
    const Center& center, const Double radius, std::vector<Cell>& cells) {
    CacheQueryState(center, radius);

    CandidateCollector collector{*this};
    collector.CollectFlat(center, radius, cells);
}

template <typename Double>
Double HyperbolicGeometryPolicy<Double>::Radius(const Center& /*unused*/) const {
    return gen_.current_hyperedge_radius_;
}

template <typename Double>
CellBallRelation HyperbolicGeometryPolicy<Double>::ClassifyCell(
    const Center& /*center*/, const Double /*radius*/, const Cell& cell) const {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    static std::uint64_t calls   = 0;
    static std::uint64_t outside = 0;
    static std::uint64_t inside  = 0;
    static std::uint64_t partial = 0;
#endif
    CellBallRelation result;

    if (CellAABBOutsideBall(cell)) {
        result = CellBallRelation::OUTSIDE;
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        ++outside;
#endif
    } else if (ShouldTryInside(cell) && HyperbolicCellInside(cell)) {
        result = CellBallRelation::INSIDE;
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        ++inside;
#endif
    } else {
        result = CellBallRelation::PARTIAL;
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        ++partial;
#endif
    }
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    ++calls;

    if (calls % 100000 == 0) {
        std::cerr << "[cell classification]"
                  << " calls=" << calls << " inside=" << static_cast<double>(inside) / static_cast<double>(calls)
                  << " partial=" << static_cast<double>(partial) / static_cast<double>(calls)
                  << " outside=" << static_cast<double>(outside) / static_cast<double>(calls) << '\n';
    }
#endif
    return result;
}

template <typename Double>
CellBallRelation HyperbolicGeometryPolicy<Double>::ClassifyRegion(const CellRegion& region) const {
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

template <typename Double>
HyperbolicGeometryPolicy<Double>::CellRegion HyperbolicGeometryPolicy<Double>::MakeRegion(
    const SInt first_annulus, const SInt last_annulus, const Double min_phi, const Double max_phi) const {
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
template <typename Double>
std::pair<typename HyperbolicGeometryPolicy<Double>::CellRegion, typename HyperbolicGeometryPolicy<Double>::CellRegion>
HyperbolicGeometryPolicy<Double>::SplitRegionRadially(const CellRegion& region) const {
    if (region.first_annulus >= region.last_annulus) {
        throw std::logic_error("SplitRegionRadially: cannot split single-annulus region");
    }

    const SInt mid_annulus = region.first_annulus + (region.last_annulus - region.first_annulus) / 2;

    CellRegion inner = MakeRegion(region.first_annulus, mid_annulus, region.min_phi, region.max_phi);

    CellRegion outer = MakeRegion(mid_annulus + 1, region.last_annulus, region.min_phi, region.max_phi);

    return {std::move(inner), std::move(outer)};
}

template <typename Double>
std::pair<typename HyperbolicGeometryPolicy<Double>::CellRegion, typename HyperbolicGeometryPolicy<Double>::CellRegion>
HyperbolicGeometryPolicy<Double>::SplitRegionAngularly(const CellRegion& region) const {
    const Double mid_phi = (region.min_phi + region.max_phi) / Double{2.0};

    if (!(region.min_phi < mid_phi && mid_phi < region.max_phi)) {
        throw std::logic_error("SplitRegionAngularly: cannot split angular interval");
    }

    CellRegion left = MakeRegion(region.first_annulus, region.last_annulus, region.min_phi, mid_phi);

    CellRegion right = MakeRegion(region.first_annulus, region.last_annulus, mid_phi, region.max_phi);

    return {std::move(left), std::move(right)};
}

template <typename Double>
bool HyperbolicGeometryPolicy<Double>::IsLeaf(const CellRegion& region) const {
    if (region.first_annulus != region.last_annulus) {
        return false;
    }

    const auto [first_cell, last_cell] =
        gen_.GlobalCellRangeForAngularInterval(region.first_annulus, region.min_phi, region.max_phi);

    return first_cell == last_cell;
}

template <typename Double>
Double HyperbolicGeometryPolicy<Double>::CellCoverage(
    const Center& center, const Double /*hyperball_radius*/, const Cell& cell) const {
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

template <typename Double>
SInt HyperbolicGeometryPolicy<Double>::AddWholeCell(const Cell& cell, std::vector<PinRange>& ranges) const {
    const CellBlock& block = GetCellBlock(cell.annulus_id, cell.chunk_id);

    if (cell.cell_id < 0 || static_cast<std::size_t>(cell.cell_id) >= block.cells.size()) {
        return 0;
    }

    const auto& stored_cell = block.cells[static_cast<std::size_t>(cell.cell_id)];
    const SInt  size        = std::get<0>(stored_cell);
    const SInt  offset      = std::get<4>(stored_cell);

    if (size > 0) {
        ranges.push_back({.begin = offset, .end = offset + size});
    }

    return size;
}

template <typename Double>
SInt HyperbolicGeometryPolicy<Double>::AddPartialCellRange(
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

template <typename Double>
SInt HyperbolicGeometryPolicy<Double>::AddPartialCellFloyd(
    const Center& /*center*/, const Cell& cell, const Double coverage, std::vector<SInt>& pins,
    std::vector<PinRange>& /*ranges*/) const {
    const auto info = GetPartialCellSampleInfo(cell, coverage);
    if (!info) {
        return 0;
    }

    FloydSampleGeometricAppend(info->offset, info->size, info->k, rng_, pins, floyd_scratch_);
    return info->k;
}

template <typename Double>
SInt HyperbolicGeometryPolicy<Double>::AddPartialCellExact(
    const Center& /*center*/, const Double /*radius*/, const Cell& cell, std::vector<SInt>& pins) const {
    const auto& vertices = ExactVertices(cell);

    return gen_.config_.debug ? AddExactVerticesChecked(vertices, pins) : AddExactVerticesFast(vertices, pins);
}

template <typename Double>
void HyperbolicGeometryPolicy<Double>::EmitHyperedge(
    const std::vector<SInt>& pins, const std::vector<PinRange>& ranges) {
    gen_.PushHyperedgeCompressed(pins, ranges);
}

template <typename Double>
std::string HyperbolicGeometryPolicy<Double>::CenterToString(const Center& center) const {
    std::ostringstream out;

    out << "phi=" << center.phi << ";r=" << center.r;

    return out.str();
}

template <typename Double>
bool HyperbolicGeometryPolicy<Double>::ShouldApproximatePartialCell(const Cell& cell) const {
    return !IsLocalCell(cell) && !exact_vertices_.contains(cell.global_cell_id);
}

// ===== Query state =====
template <typename Double>
void HyperbolicGeometryPolicy<Double>::CacheQueryState(const Center& center, const Double radius) {
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

template <typename Double>
HyperbolicGeometryPolicy<Double>::CellAnnulusRegion HyperbolicGeometryPolicy<Double>::MakeCellAnnulusRegion(
    const SInt annulus_id, const SInt first_cell, const SInt end_cell) const {
    if (annulus_id >= gen_.total_annuli_) {
        throw std::out_of_range("MakeCellAnnulusRegion: invalid annulus");
    }

    const SInt total_cells = gen_.global_cells_per_annulus_[annulus_id];

    if (first_cell >= end_cell || end_cell > total_cells) {
        throw std::logic_error("MakeCellAnnulusRegion: invalid cell interval");
    }

    const Double cell_width = Double{2.0 * M_PI} / static_cast<Double>(total_cells);

    const Double min_phi = static_cast<Double>(first_cell) * cell_width;

    const Double max_phi = static_cast<Double>(end_cell) * cell_width;

    const Double min_r = gen_.annulus_min_r_[annulus_id];

    const Double max_r = gen_.annulus_max_r_[annulus_id];

    const auto box = poincare_geometry::ComputeCellAABB(min_r, max_r, min_phi, max_phi, gen_.cell_eps_);

    return CellAnnulusRegion{
        .annulus_id = annulus_id,
        .first_cell = first_cell,
        .end_cell   = end_cell,
        .min_r      = min_r,
        .max_r      = max_r,
        .min_phi    = min_phi,
        .max_phi    = max_phi,
        .min_x      = box.min_x,
        .max_x      = box.max_x,
        .min_y      = box.min_y,
        .max_y      = box.max_y,
    };
}

template <typename Double>
std::pair<
    typename HyperbolicGeometryPolicy<Double>::CellAnnulusRegion,
    typename HyperbolicGeometryPolicy<Double>::CellAnnulusRegion>
HyperbolicGeometryPolicy<Double>::SplitCellAnnulusRegion(const CellAnnulusRegion& region) const {
    if (region.end_cell - region.first_cell <= 1) {
        throw std::logic_error("SplitCellAnnulusRegion: cannot split leaf region");
    }

    const SInt mid = region.first_cell + ((region.end_cell - region.first_cell) / 2);

    CellAnnulusRegion left = MakeCellAnnulusRegion(region.annulus_id, region.first_cell, mid);

    CellAnnulusRegion right = MakeCellAnnulusRegion(region.annulus_id, mid, region.end_cell);

    return {std::move(left), std::move(right)};
}

template <typename Double>
bool HyperbolicGeometryPolicy<Double>::IsLeaf(const CellAnnulusRegion& region) const {
    return region.end_cell - region.first_cell == 1;
}

template <typename Double>
CellBallRelation HyperbolicGeometryPolicy<Double>::ClassifyRegion(const CellAnnulusRegion& region) const {
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

template <typename Double>
void HyperbolicGeometryPolicy<Double>::EmitInsideRegion(
    const CellAnnulusRegion& region, std::vector<PinRange>& inside_ranges) const {
    for (SInt global_cell = region.first_cell; global_cell < region.end_cell; ++global_cell) {
        const auto [chunk_id, local_cell_id] = gen_.GlobalCellToChunkCell(region.annulus_id, global_cell);

        const CellBlock& block = GetCellBlock(region.annulus_id, chunk_id);

        if (local_cell_id < 0 || static_cast<std::size_t>(local_cell_id) >= block.cells.size()) {
            continue;
        }

        const auto& stored_cell = block.cells[static_cast<std::size_t>(local_cell_id)];

        const SInt size   = std::get<0>(stored_cell);
        const SInt offset = std::get<4>(stored_cell);

        if (size > 0) {
            inside_ranges.push_back({
                .begin = offset,
                .end   = offset + size,
            });
        }
    }
}

template <typename Double>
void HyperbolicGeometryPolicy<Double>::TraverseCandidateRegion(
    const CellAnnulusRegion& region, std::vector<Cell>& cells, std::vector<PinRange>& inside_ranges,
    CandidateCollector& collector) {
    const auto relation = ClassifyRegion(region);

    if (relation == CellBallRelation::OUTSIDE) {
        return;
    }

    if (relation == CellBallRelation::INSIDE) {
        EmitInsideRegion(region, inside_ranges);
        return;
    }

    if (IsLeaf(region)) {
        collector.AddLeafCell(region, cells);
        return;
    }

    auto [left, right] = SplitCellAnnulusRegion(region);

    TraverseCandidateRegion(left, cells, inside_ranges, collector);
    TraverseCandidateRegion(right, cells, inside_ranges, collector);
}

template <typename Double>
void HyperbolicGeometryPolicy<Double>::TraverseCandidateRegion(
    const CellRegion& region, std::vector<Cell>& cells, std::vector<PinRange>& inside_ranges,
    CandidateCollector& collector) {
    const auto relation = ClassifyRegion(region);

    if (relation == CellBallRelation::OUTSIDE) {
        return;
    }

    if (relation == CellBallRelation::INSIDE && region.first_annulus != region.last_annulus) {
        auto [inner, outer] = SplitRegionRadially(region);

        TraverseCandidateRegion(inner, cells, inside_ranges, collector);
        TraverseCandidateRegion(outer, cells, inside_ranges, collector);
        return;
    }

    if (IsLeaf(region)) {
        collector.AddLeafCell(region, cells);
        return;
    }

    if (region.first_annulus != region.last_annulus) {
        auto [inner, outer] = SplitRegionRadially(region);
        TraverseCandidateRegion(inner, cells, inside_ranges, collector);
        TraverseCandidateRegion(outer, cells, inside_ranges, collector);
    } else {
        auto [left, right] = SplitRegionAngularly(region);
        TraverseCandidateRegion(left, cells, inside_ranges, collector);
        TraverseCandidateRegion(right, cells, inside_ranges, collector);
    }
}

template <typename Double>
Double
HyperbolicGeometryPolicy<Double>::AllowedHalfAngleForAnnulus(const Double /*center_r*/, const SInt annulus_id) const {
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

template <typename Double>
Double HyperbolicGeometryPolicy<Double>::AllowedHalfAngleForCachedRadius(
    const Double query_cosh, const Double query_sinh) const {
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
template <typename Double>
bool HyperbolicGeometryPolicy<Double>::CellAABBOutsideBall(const Cell& cell) const {
    const poincare_geometry::AABB<Double> box{
        .min_x = cell.min_x,
        .max_x = cell.max_x,
        .min_y = cell.min_y,
        .max_y = cell.max_y,
    };

    return poincare_geometry::AABBOutsideBall(box, ball_);
}

template <typename Double>
bool HyperbolicGeometryPolicy<Double>::HyperbolicCellInside(const Cell& cell) const {
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

template <typename Double>
bool HyperbolicGeometryPolicy<Double>::ShouldTryInside(const Cell& cell) const {
    return gen_.current_hyperedge_radius_ >= (cell.max_r - cell.min_r);
}

template <typename Double>
Double HyperbolicGeometryPolicy<Double>::MaxAngularDistanceToInterval(
    const Double phi, const Double min_phi, const Double max_phi) const {
    if (max_phi - min_phi >= Double{2.0 * M_PI} - gen_.cell_eps_) {
        return Double{M_PI};
    }

    const Double antipode = circular_interval::NormalizePhi(phi + Double{M_PI});

    if (circular_interval::AngleInInterval(antipode, min_phi, max_phi)) {
        return Double{M_PI};
    }

    return std::max(circular_interval::AngularDistance(phi, min_phi), circular_interval::AngularDistance(phi, max_phi));
}

// ===== Coverage estimation =====

template <typename Double>
Double HyperbolicGeometryPolicy<Double>::AllowedHalfAngleAtRadius(
    const Double query_r, const Double query_cosh, const Double query_sinh) const {
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

template <typename Double>
Double
HyperbolicGeometryPolicy<Double>::MinAngularDistanceToInterval(Double phi, Double min_phi, Double max_phi) const {
    if (circular_interval::AngleInInterval(phi, min_phi, max_phi)) {
        return 0.0;
    }

    return std::min(circular_interval::AngularDistance(phi, min_phi), circular_interval::AngularDistance(phi, max_phi));
}

template <typename Double>
Double HyperbolicGeometryPolicy<Double>::CoshDistanceWithDelta(Double query_r, Double delta_phi) const {
    return (center_cosh_r_ * std::cosh(query_r)) - (center_sinh_r_ * std::sinh(query_r) * std::cos(delta_phi));
}

// ===== Approx Partial-Cell Sampling =====
template <typename Double>
std::optional<typename HyperbolicGeometryPolicy<Double>::PartialCellSampleInfo>
HyperbolicGeometryPolicy<Double>::GetPartialCellSampleInfo(const Cell& cell, const Double coverage) const {
    const SInt global_cell_id = cell.global_cell_id;

    const CellBlock& block = GetCellBlock(cell.annulus_id, cell.chunk_id);

    if (cell.cell_id < 0 || static_cast<std::size_t>(cell.cell_id) >= block.cells.size()) {
        return std::nullopt;
    }

    const auto& stored_cell = block.cells[static_cast<std::size_t>(cell.cell_id)];

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

template <typename Double>
const Hyper_Hyperbolic<Double>::VertexBlock& HyperbolicGeometryPolicy<Double>::ExactVertices(const Cell& cell) const {
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

template <typename Double>
const HyperbolicGeometryPolicy<Double>::CachedExactCell&
HyperbolicGeometryPolicy<Double>::ExactRemoteCell(const Cell& cell) const {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    RecordRemoteAccess(cell.global_cell_id);
#endif
    auto it = exact_vertices_.find(cell.global_cell_id);
    if (it != exact_vertices_.end()) {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        ++exact_remote_cache_hits_;
#endif
        exact_lru_.splice(exact_lru_.begin(), exact_lru_, exact_lru_pos_[cell.global_cell_id]);
        return it->second;
    }
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    ++exact_remote_cache_misses_;
#endif
    auto [inserted_it, inserted] = exact_vertices_.emplace(cell.global_cell_id, CachedExactCell{});

    exact_lru_.push_front(cell.global_cell_id);
    exact_lru_pos_[cell.global_cell_id] = exact_lru_.begin();

    auto& cached = inserted_it->second;

    const CellBlock& block = GetCellBlock(cell.annulus_id, cell.chunk_id);

    assert(cell.cell_id >= 0);
    assert(static_cast<std::size_t>(cell.cell_id) < block.cells.size());

    const auto& stored_cell = block.cells[static_cast<std::size_t>(cell.cell_id)];

    gen_.GenerateVertices(cell.annulus_id, cell.chunk_id, cell.cell_id, block.annulus, stored_cell, cached.vertices);

    exact_remote_cached_bytes_ += ExactCellBytes(cached.vertices);
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    exact_remote_cached_vertices_ += static_cast<SInt>(cached.vertices.size());
#endif
    while (exact_remote_cached_bytes_ > exact_remote_cache_budget_ && exact_vertices_.size() > 1) {
        const SInt victim = exact_lru_.back();

        const auto victim_it = exact_vertices_.find(victim);
        if (victim_it != exact_vertices_.end()) {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
            exact_remote_cached_vertices_ -= static_cast<SInt>(victim_it->second.vertices.size());
#endif
            exact_remote_cached_bytes_ -= ExactCellBytes(victim_it->second.vertices);
            exact_vertices_.erase(victim_it);
        }

        exact_lru_.pop_back();
        exact_lru_pos_.erase(victim);
    }

    return cached;
}

template <typename Double>
bool HyperbolicGeometryPolicy<Double>::IsLocalCell(const Cell& cell) const {
    return gen_.IsLocalChunk(cell.chunk_id);
}

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
template <typename Double>
void HyperbolicGeometryPolicy<Double>::RecordRemoteAccess(const SInt global_cell_id) const {
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

template <typename Double>
void HyperbolicGeometryPolicy<Double>::PrintExactCacheStats() const {
    std::cerr << " exact_remote_hits=" << exact_remote_cache_hits_
              << " exact_remote_misses=" << exact_remote_cache_misses_
              << " exact_remote_cached_vertices=" << exact_remote_cached_vertices_
              << " exact_remote_cached_bytes=" << exact_remote_cached_bytes_
              << " exact_remote_cache_budget=" << exact_remote_cache_budget_ << '\n';

    const double avg_reuse = (exact_remote_reuse_count_ != 0u) ? static_cast<double>(exact_remote_reuse_distance_sum_)
                                                                     / static_cast<double>(exact_remote_reuse_count_)
                                                               : 0.0;

    std::cerr << " remote_reuse_count=" << exact_remote_reuse_count_ << " remote_reuse_avg=" << avg_reuse
              << " remote_reuse_max=" << exact_remote_reuse_distance_max_
              << " reuse<=1=" << exact_remote_reuse_distance_le_1_ << " reuse<=4=" << exact_remote_reuse_distance_le_4_
              << " reuse<=16=" << exact_remote_reuse_distance_le_16_
              << " reuse>16=" << exact_remote_reuse_distance_gt_16_ << '\n';

    std::cerr << " cell_block_count=" << cell_blocks_.size() << " cell_block_cached_bytes=" << cell_block_cached_bytes_
              << " cell_block_cache_budget=" << cell_block_cache_budget_ << '\n';

    std::cerr << "[HyperRHG capacity]"
              << " edge_offsets=" << gen_.graph_.hyperedge_offsets.capacity()
              << " pins=" << gen_.graph_.hyperedge_pins.capacity()
              << " pin_ranges=" << gen_.graph_.hyperedge_ranges.capacity()
              << " range_offsets=" << gen_.graph_.hyperedge_range_offsets.capacity() << '\n';

    std::cerr << "[HyperRHG maps]"
              << " cells_size=" << gen_.cells_.size() << " cells_buckets=" << gen_.cells_.bucket_count()
              << " center_cells_size=" << gen_.center_cells_.size()
              << " center_cells_buckets=" << gen_.center_cells_.bucket_count()
              << " vertices_size=" << gen_.vertices_.size() << " vertices_buckets=" << gen_.vertices_.bucket_count()
              << " annuli_size=" << gen_.annuli_.size() << " center_annuli_size=" << gen_.center_annuli_.size() << '\n';
}
#endif

template <typename Double>
bool HyperbolicGeometryPolicy<Double>::IsInsideHyperbolicBallFast(const Vertex& vertex) const {
    return PGGeometry<Double>::HyperbolicDistance(center_vertex_, vertex) <= gen_.current_hyperedge_pdm_radius_;
}

template <typename Double>
SInt HyperbolicGeometryPolicy<Double>::AddExactVerticesFast(
    const typename GeneratorT::VertexBlock& vertices, std::vector<SInt>& pins) const {
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

template <typename Double>
SInt HyperbolicGeometryPolicy<Double>::AddExactVerticesChecked(
    const typename GeneratorT::VertexBlock& vertices, std::vector<SInt>& pins) const {
    return AddExactVerticesFast(vertices, pins);
}

template <typename Double>
void HyperbolicGeometryPolicy<Double>::EnsureAnnulusMidpoints() {
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

template <typename Double>
std::size_t HyperbolicGeometryPolicy<Double>::ExactCellBytes(const GeneratorT::VertexBlock& vertices) const {
    return (sizeof(Double)
            * (vertices.phi.capacity() + vertices.r.capacity() + vertices.x.capacity() + vertices.y.capacity()
               + vertices.gamma.capacity() + vertices.cosh_r.capacity() + vertices.sinh_r.capacity()
               + vertices.cos_phi.capacity() + vertices.sin_phi.capacity()))
           + (sizeof(SInt) * vertices.id.capacity());
}

template <typename Double>
typename HyperbolicGeometryPolicy<Double>::CellBlock
HyperbolicGeometryPolicy<Double>::BuildCellBlock(const SInt annulus_id, const SInt chunk_id) const {
    typename GeneratorT::Chunk   chunk;
    typename GeneratorT::Annulus annulus;

    if (gen_.IsLocalChunk(chunk_id)) {
        const SInt global_chunk_id = gen_.ComputeGlobalChunkId(annulus_id, chunk_id);

        const auto chunk_it = gen_.chunks_.find(chunk_id);

        const auto annulus_it = gen_.annuli_.find(global_chunk_id);

        assert(chunk_it != gen_.chunks_.end());
        assert(annulus_it != gen_.annuli_.end());

        chunk   = chunk_it->second;
        annulus = annulus_it->second;
    } else {
        const auto metadata = gen_.ReconstructChunkAnnulus(annulus_id, chunk_id);

        chunk   = metadata.chunk;
        annulus = metadata.annulus;
    }

    SInt size   = std::get<0>(annulus);
    SInt offset = std::get<4>(annulus);

    const Double min_phi = std::get<1>(chunk);
    const Double max_phi = std::get<2>(chunk);

    Double remaining_phi = max_phi - min_phi;

    const SInt num_cells = gen_.CellsPerChunkForAnnulus(annulus_id, chunk_id);

    CellBlock block{
        .annulus_id = annulus_id,
        .chunk_id   = chunk_id,
        .cells      = {},
        .annulus    = annulus,
    };

    block.cells.reserve(num_cells);

    if (num_cells == 0) {
        return block;
    }

    const Double cell_width = (max_phi - min_phi) / static_cast<Double>(num_cells);

    for (SInt cell_id = 0; cell_id < num_cells; ++cell_id) {
        const SInt seed = gen_.config_.seed + (annulus_id * gen_.config_.k) + chunk_id + cell_id + size;

        const SInt hash_value = sampling::Spooky::hash(seed);

        const SInt cell_size = gen_.rng_.GenerateBinomial(
            hash_value, size, std::clamp(cell_width / remaining_phi, Double{0.0}, Double{1.0}));

        const SInt global_cell = gen_.ChunkCellToGlobalCell(annulus_id, chunk_id, cell_id);

        const Double global_cell_width =
            Double{2.0 * M_PI} / static_cast<Double>(gen_.global_cells_per_annulus_[annulus_id]);

        const Double cell_min_phi = static_cast<Double>(global_cell) * global_cell_width;

        const Double cell_max_phi = static_cast<Double>(global_cell + 1) * global_cell_width;

        const Double annulus_min_r = std::get<1>(annulus);
        const Double annulus_max_r = std::get<2>(annulus);

        const auto box = poincare_geometry::ComputeCellAABB(
            annulus_min_r, annulus_max_r, cell_min_phi, cell_max_phi, gen_.cell_eps_);

        block.cells.emplace_back(
            cell_size, cell_min_phi, cell_max_phi, false, offset, box.min_x, box.max_x, box.min_y, box.max_y);

        size -= cell_size;
        offset += cell_size;
        remaining_phi -= cell_width;
    }
    return block;
}

template <typename Double>
const typename HyperbolicGeometryPolicy<Double>::CellBlock&
HyperbolicGeometryPolicy<Double>::GetCellBlock(const SInt annulus_id, const SInt chunk_id) const {
    const CellBlockKey key{annulus_id, chunk_id};

    auto it = cell_blocks_.find(key);
    if (it != cell_blocks_.end()) {
        cell_block_lru_.splice(cell_block_lru_.begin(), cell_block_lru_, it->second.lru_position);

        it->second.lru_position = cell_block_lru_.begin();
        return it->second.block;
    }

    CellBlock         block = BuildCellBlock(annulus_id, chunk_id);
    const std::size_t bytes = CellBlockBytes(block);

    cell_block_lru_.push_front(key);

    auto [inserted_it, inserted] = cell_blocks_.emplace(
        key, CachedCellBlock{
                 .block        = std::move(block),
                 .lru_position = cell_block_lru_.begin(),
                 .bytes        = bytes,
             });

    assert(inserted);
    cell_block_cached_bytes_ += bytes;

    while (cell_block_cached_bytes_ > cell_block_cache_budget_ && cell_blocks_.size() > 1) {
        const CellBlockKey victim    = cell_block_lru_.back();
        const auto         victim_it = cell_blocks_.find(victim);

        assert(victim_it != cell_blocks_.end());

        cell_block_cached_bytes_ -= victim_it->second.bytes;
        cell_blocks_.erase(victim_it);
        cell_block_lru_.pop_back();
    }

    return inserted_it->second.block;
}

template <typename Double>
std::size_t HyperbolicGeometryPolicy<Double>::CellBlockBytes(const CellBlock& block) const {
    return sizeof(CellBlock) + (block.cells.capacity() * sizeof(typename GeneratorT::Cell));
}

template class HyperbolicGeometryPolicy<LPFloat>;
template class HyperbolicGeometryPolicy<HPFloat>;

} // namespace kagen