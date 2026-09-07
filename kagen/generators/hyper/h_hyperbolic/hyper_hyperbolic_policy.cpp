
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
    collector.CollectRadialHierarchy(center, radius, cells, ranges);

    return ranges.size() > ranges_before;
}

template <typename Double>
void HyperbolicGeometryPolicy<Double>::CandidateCollector::TraverseSingleAnnulus(
    const SInt annulus_id, std::vector<Cell>& cells, std::vector<PinRange>& ranges) {
    const SInt total_cells = gen().global_cells_per_annulus_[annulus_id];

    if (total_cells <= 0) {
        return;
    }

    const Double half_angle = policy.AllowedHalfAngleForAnnulus(policy.center_r_, annulus_id);

    if (!(half_angle > Double{0.0})) {
        return;
    }

    if (half_angle >= Double{M_PI}) {
        const CellAnnulusRegion root = policy.MakeCellAnnulusRegion(annulus_id, 0, total_cells);

        policy.TraverseCandidateRegion(root, cells, ranges, *this);

        return;
    }

    const auto parts =
        circular_interval::Split(policy.center_phi_ - half_angle, policy.center_phi_ + half_angle, gen().cell_eps_);

    for (int part = 0; part < parts.count; ++part) {
        const Double q_begin = parts.parts[part].first;
        const Double q_end   = parts.parts[part].second;

        if (!(q_begin < q_end)) {
            continue;
        }

        const auto [first_cell, last_cell] = gen().GlobalCellRangeForAngularInterval(annulus_id, q_begin, q_end);

        const SInt first = std::clamp<SInt>(first_cell, 0, total_cells - 1);

        const SInt end = std::clamp<SInt>(last_cell + 1, 1, total_cells);

        if (first >= end) {
            continue;
        }

        const CellAnnulusRegion root = policy.MakeCellAnnulusRegion(annulus_id, first, end);

        policy.TraverseCandidateRegion(root, cells, ranges, *this);
    }
}

template <typename Double>
void HyperbolicGeometryPolicy<Double>::EmitInsideRegion(
    const CellRegion& region, std::vector<PinRange>& inside_ranges) const {
    for (SInt annulus_id = region.first_annulus; annulus_id <= region.last_annulus; ++annulus_id) {
        const auto [first_cell, last_cell] =
            gen_.GlobalCellRangeForAngularInterval(annulus_id, region.min_phi, region.max_phi);

        const SInt end_cell = last_cell + 1;

        const SInt cells_per_chunk = gen_.global_cells_per_annulus_[annulus_id] / gen_.config_.k;

        const SInt first_chunk = first_cell / cells_per_chunk;

        const SInt last_chunk = (end_cell - 1) / cells_per_chunk;

        for (SInt chunk_id = first_chunk; chunk_id <= last_chunk; ++chunk_id) {
            const SInt chunk_begin = chunk_id * cells_per_chunk;

            const SInt chunk_end = chunk_begin + cells_per_chunk;

            const SInt intersection_begin = std::max(first_cell, chunk_begin);

            const SInt intersection_end = std::min(end_cell, chunk_end);

            if (intersection_begin < intersection_end) {
                EmitInsideChunkIntersection(annulus_id, chunk_id, intersection_begin, intersection_end, inside_ranges);
            }
        }
    }
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
std::pair<SInt, SInt>
HyperbolicGeometryPolicy<Double>::ReachableAnnulusRange(const Center& center, const Double radius) const {
    if (gen_.total_annuli_ <= 0) {
        return {1, 0}; // empty interval
    }

    const Double width = gen_.target_r_ / static_cast<Double>(gen_.total_annuli_);

    const Double min_r = std::max<Double>(Double{0.0}, center.r - radius);

    const Double max_r = std::min<Double>(gen_.target_r_, center.r + radius);

    if (min_r > max_r) {
        return {1, 0};
    }

    SInt first = static_cast<SInt>(std::floor(min_r / width));

    SInt last = static_cast<SInt>(std::floor(max_r / width));

    first = std::clamp<SInt>(first, SInt{0}, gen_.total_annuli_ - SInt{1});

    last = std::clamp<SInt>(last, SInt{0}, gen_.total_annuli_ - SInt{1});

    // The replicated inner annuli are handled separately.
    first = std::max<SInt>(first, gen_.replicated_inner_last_annulus_ + SInt{1});

    return {first, last};
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

        std::cerr << "[HRHG hierarchy] cell_region_visits=" << cell_region_visits_
                  << " cell_annulus_region_visits=" << cell_annulus_region_visits_ << '\n';
    }
#endif
    return result;
}

template <typename Double>
CellBallRelation HyperbolicGeometryPolicy<Double>::ClassifyRegion(const CellRegion& region) const {
    if (region.min_r > center_r_ + gen_.current_hyperedge_radius_
        || region.max_r < center_r_ - gen_.current_hyperedge_radius_) {
        return CellBallRelation::OUTSIDE;
    }
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
    if (last_annulus >= gen_.total_annuli_ || first_annulus > last_annulus) {
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
SInt HyperbolicGeometryPolicy<Double>::AddReplicatedInnerVertices(
    const Center& center, const Double radius, std::vector<SInt>& pins) {
    if (gen_.replicated_inner_last_annulus_ < 0 || gen_.replicated_inner_vertices_.id.empty()) {
        return 0;
    }

    //
    // Make sure all cached query geometry corresponds
    // exactly to this center/radius.
    //
    CacheQueryState(center, radius);

    //
    // Entire origin-centered replicated disk lies inside
    // the query hyperball.
    //
    if (center_r_ + gen_.replicated_inner_radius_ <= radius) {
        pins.insert(pins.end(), gen_.replicated_inner_vertices_.id.begin(), gen_.replicated_inner_vertices_.id.end());

        return static_cast<SInt>(gen_.replicated_inner_vertices_.id.size());
    }

    //
    // Quick radial nonintersection.
    //
    if (center_r_ - radius > gen_.replicated_inner_radius_) {
        return 0;
    }

    //
    // Partial overlap: exact-test the bounded replica.
    //
    return gen_.config_.debug ? AddExactVerticesChecked(gen_.replicated_inner_vertices_, pins)
                              : AddExactVerticesFast(gen_.replicated_inner_vertices_, pins);
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
    if (region.min_r > center_r_ + gen_.current_hyperedge_radius_
        || region.max_r < center_r_ - gen_.current_hyperedge_radius_) {
        return CellBallRelation::OUTSIDE;
    }
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
    const SInt cells_per_chunk = gen_.global_cells_per_annulus_[region.annulus_id] / gen_.config_.k;

    const SInt first_chunk = region.first_cell / cells_per_chunk;

    const SInt last_chunk = (region.end_cell - 1) / cells_per_chunk;

    for (SInt chunk_id = first_chunk; chunk_id <= last_chunk; ++chunk_id) {
        const SInt chunk_first = chunk_id * cells_per_chunk;

        const SInt chunk_end = chunk_first + cells_per_chunk;

        const SInt first = std::max(region.first_cell, chunk_first);

        const SInt end = std::min(region.end_cell, chunk_end);

        EmitInsideChunkIntersection(region.annulus_id, chunk_id, first, end, inside_ranges);
    }
}

template <typename Double>
void HyperbolicGeometryPolicy<Double>::EmitInsideChunkIntersection(
    const SInt annulus_id, const SInt chunk_id, const SInt first_global_cell, const SInt end_global_cell,
    std::vector<PinRange>& inside_ranges) const {
    const SInt cells_per_chunk = gen_.global_cells_per_annulus_[annulus_id] / gen_.config_.k;

    const SInt chunk_first = chunk_id * cells_per_chunk;

    const SInt chunk_end = chunk_first + cells_per_chunk;

    //
    // Entire chunk-annulus accepted.
    //
    if (first_global_cell == chunk_first && end_global_cell == chunk_end) {
        typename GeneratorT::Annulus annulus;

        if (gen_.IsLocalChunk(chunk_id)) {
            const auto id = gen_.ComputeGlobalChunkId(annulus_id, chunk_id);

            annulus = gen_.annuli_[id];
        } else {
            annulus = gen_.ReconstructChunkAnnulus(annulus_id, chunk_id).annulus;
        }

        const SInt size   = std::get<0>(annulus);
        const SInt offset = std::get<4>(annulus);

        if (size > 0) {
            inside_ranges.push_back({
                .begin = offset,
                .end   = offset + size,
            });
        }

        return;
    }

    //
    // Only a boundary part of this chunk is accepted.
    // Here we still need the cells.
    //
    const CellBlock& block = GetCellBlock(annulus_id, chunk_id);

    for (SInt global_cell = first_global_cell; global_cell < end_global_cell; ++global_cell) {
        const SInt local_cell = global_cell - chunk_first;

        const auto& stored_cell = block.cells[static_cast<std::size_t>(local_cell)];

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
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    ++cell_annulus_region_visits_;
#endif

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
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    ++cell_region_visits_;
#endif

    const auto relation = ClassifyRegion(region);

    if (relation == CellBallRelation::OUTSIDE) {
        return;
    }

    if (relation == CellBallRelation::INSIDE) {
        EmitInsideRegion(region, inside_ranges);
        return;
    }

    if (region.first_annulus != region.last_annulus) {
        auto [inner, outer] = SplitRegionRadially(region);

        TraverseCandidateRegion(inner, cells, inside_ranges, collector);
        TraverseCandidateRegion(outer, cells, inside_ranges, collector);
        return;
    }

    collector.TraverseSingleAnnulus(region.first_annulus, cells, inside_ranges);
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
void HyperbolicGeometryPolicy<Double>::CandidateCollector::CollectRadialHierarchy(
    const Center& center, const Double radius, std::vector<Cell>& cells, std::vector<PinRange>& ranges) {
    cells.clear();
    seen_candidate_cells_.clear();

    const auto [first_annulus, last_annulus] = policy.ReachableAnnulusRange(center, radius);

    if (first_annulus > last_annulus) {
        return;
    }

    const CellRegion root = policy.MakeRegion(first_annulus, last_annulus, Double{0.0}, Double{2.0 * M_PI});

    policy.TraverseCandidateRegion(root, cells, ranges, *this);
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
    //
    // Hyperbolic distance is at least the absolute difference
    // of the radial coordinates. Hence an annulus that is
    // radially disjoint from the query ball cannot intersect it,
    // irrespective of its angular extent.
    //
    if (cell.min_r > center_r_ + gen_.current_hyperedge_radius_
        || cell.max_r < center_r_ - gen_.current_hyperedge_radius_) {
        return true;
    }

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