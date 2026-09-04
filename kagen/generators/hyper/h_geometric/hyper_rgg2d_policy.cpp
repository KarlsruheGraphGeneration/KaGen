#include "kagen/generators/hyper/h_geometric/hyper_rgg2d_policy.h"

#include "kagen/generators/hyper/h_geometric/hyper_rgg2d.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <optional>
#include <vector>

namespace kagen {

HyperRGG2DPolicy::HyperRGG2DPolicy(HyperRGG2D& generator) : gen_(&generator), rng_(RNGWrapper(generator.config_)) {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    local_exact_last_access_.set_empty_key(std::numeric_limits<SInt>::max());
#endif
}

void HyperRGG2DPolicy::AddCenter(
    const Center& center, std::vector<SInt>& /*pins*/
) const {
    const SInt sampling_seed = sampling::Spooky::hash(gen_->config_.seed + (131 * center.sampled_id));

    switch (gen_->config_.partial_cell_mode) {
        case PartialCellMode::EstimateByCoverageRange:
            gen_->mersenne.RandomInit(sampling_seed);
            break;

        case PartialCellMode::GenerateAndCheck:
            // Exact mode performs no approximate random sampling.
            break;

        case PartialCellMode::EstimateByCoverageFloyd:
            rng_.SeedUniformStream(sampling_seed);
            break;
    }
}

SInt HyperRGG2DPolicy::GetNumVerticeOfCellCoord(const SSInt global_cell_x, const SSInt global_cell_y) {
    const auto cell = TryMakeCell(global_cell_x, global_cell_y);
    if (!cell) {
        return SInt{0};
    }

    const auto cell_it = gen_->cells_.find(cell->global_cell_id);
    if (cell_it == gen_->cells_.end()) {
        return SInt{0};
    }

    return std::get<0>(cell_it->second);
}

double HyperRGG2DPolicy::MinimumRadius(
    const Center& /*unused*/
) {
    return EuclideanRadiusForExpectedHyperedgeSize(2, gen_->config_.n);
}

LPFloat HyperRGG2DPolicy::Radius(const Center& center) const {
    return center.radius;
}

std::optional<HyperRGG2DPolicy::Cell>
HyperRGG2DPolicy::TryMakeCell(const SSInt global_cell_x, const SSInt global_cell_y) const {
    const auto total = static_cast<SSInt>(gen_->SafeTotalCellsPerDim());

    if (global_cell_x < 0 || global_cell_y < 0 || global_cell_x >= total || global_cell_y >= total) {
        return std::nullopt;
    }

    const SInt chunk_x = static_cast<SInt>(global_cell_x) / gen_->cells_per_dim_;
    const SInt chunk_y = static_cast<SInt>(global_cell_y) / gen_->cells_per_dim_;
    const SInt cell_x  = static_cast<SInt>(global_cell_x) % gen_->cells_per_dim_;
    const SInt cell_y  = static_cast<SInt>(global_cell_y) % gen_->cells_per_dim_;

    const SInt chunk_id = gen_->Encode(chunk_x, chunk_y);
    const SInt cell_id  = gen_->EncodeCell(cell_x, cell_y);

    return Cell{
        .chunk_id       = chunk_id,
        .cell_id        = cell_id,
        .global_cell_x  = static_cast<SInt>(global_cell_x),
        .global_cell_y  = static_cast<SInt>(global_cell_y),
        .global_cell_id = gen_->ComputeGlobalCellId(chunk_id, cell_id),
    };
}

void HyperRGG2DPolicy::CandidateCells(const Center& center, const LPFloat radius, std::vector<Cell>& cells) const {
    const SInt total_cells_per_dim = gen_->SafeTotalCellsPerDim();

    const LPFloat cell_size = 1.0 / static_cast<LPFloat>(total_cells_per_dim);

    // Replace these with the actual Center coordinate accessors.
    const LPFloat center_x = center.x;
    const LPFloat center_y = center.y;

    const auto min_cell_x = static_cast<SSInt>(std::floor((center_x - radius) / cell_size));

    const auto max_cell_x = static_cast<SSInt>(std::floor((center_x + radius) / cell_size));

    const auto min_cell_y = static_cast<SSInt>(std::floor((center_y - radius) / cell_size));

    const auto max_cell_y = static_cast<SSInt>(std::floor((center_y + radius) / cell_size));

    for (SSInt global_cell_x = min_cell_x; global_cell_x <= max_cell_x; ++global_cell_x) {
        for (SSInt global_cell_y = min_cell_y; global_cell_y <= max_cell_y; ++global_cell_y) {
            if (auto cell = TryMakeCell(global_cell_x, global_cell_y)) {
                cells.push_back(*cell);
            }
        }
    }
}

std::optional<HyperRGG2DPolicy::StoredCell> HyperRGG2DPolicy::TryGetStoredCell(const Cell& cell) const {
    if (gen_->IsLocalChunk(cell.chunk_id)) {
        const auto it = gen_->cells_.find(cell.global_cell_id);

        if (it == gen_->cells_.end()) {
            return StoredCell{
                .size   = 0,
                .offset = 0,
            };
        }

        return StoredCell{
            .size   = std::get<0>(it->second),
            .offset = std::get<4>(it->second),
        };
    }

    const auto metadata = gen_->ReconstructCellMetadata(cell.chunk_id, cell.cell_id);

    return StoredCell{
        .size   = metadata.size,
        .offset = metadata.offset,
    };
}

CellBallRelation HyperRGG2DPolicy::ClassifyCell(const Center& center, const LPFloat radius, const Cell& cell) const {
    const LPFloat radius_sq = radius * radius;
    CellBounds    bounds    = GetCellBounds(cell);

    const LPFloat closest_x = std::clamp(center.x, bounds.min_x, bounds.max_x);

    const LPFloat closest_y = std::clamp(center.y, bounds.min_y, bounds.max_y);

    const LPFloat dx_min = center.x - closest_x;
    const LPFloat dy_min = center.y - closest_y;

    if ((dx_min * dx_min) + (dy_min * dy_min) > radius_sq) {
        return CellBallRelation::OUTSIDE;
    }

    const LPFloat dx1 = center.x - bounds.min_x;
    const LPFloat dx2 = center.x - bounds.max_x;
    const LPFloat dy1 = center.y - bounds.min_y;
    const LPFloat dy2 = center.y - bounds.max_y;

    const LPFloat max_dist_sq = std::max(dx1 * dx1, dx2 * dx2) + std::max(dy1 * dy1, dy2 * dy2);

    if (max_dist_sq <= radius_sq) {
        return CellBallRelation::INSIDE;
    }

    return CellBallRelation::PARTIAL;
}

double HyperRGG2DPolicy::EstimatedCircleRectCoverage(
    const double center_x, const double center_y, const double min_x, const double max_x, const double min_y,
    const double max_y, const double radius) const {
    const double cell_area = (max_x - min_x) * (max_y - min_y);

    if (cell_area <= 0.0 || radius <= 0.0) {
        return 0.0;
    }

    const double covered_area =
        EstimateCoverageRecursive(center_x, center_y, radius * radius, min_x, max_x, min_y, max_y, coverage_max_depth_);

    return std::clamp(covered_area / cell_area, 0.0, 1.0);
}

double HyperRGG2DPolicy::EstimateCoverageRecursive(
    const double center_x, const double center_y, const double radius_sq, const double min_x, const double max_x,
    const double min_y, const double max_y, const int depth) const {
    const double width  = max_x - min_x;
    const double height = max_y - min_y;
    const double area   = width * height;

    // Closest point in the rectangle to the center.
    const double closest_x = std::clamp(center_x, min_x, max_x);
    const double closest_y = std::clamp(center_y, min_y, max_y);

    const double closest_dx = center_x - closest_x;
    const double closest_dy = center_y - closest_y;

    const double min_distance_sq = (closest_dx * closest_dx) + (closest_dy * closest_dy);

    if (min_distance_sq >= radius_sq) {
        return 0.0;
    }

    // Farthest rectangle corner from the center.
    const double dx_min = center_x - min_x;
    const double dx_max = center_x - max_x;
    const double dy_min = center_y - min_y;
    const double dy_max = center_y - max_y;

    const double max_distance_sq =
        std::max(dx_min * dx_min, dx_max * dx_max) + std::max(dy_min * dy_min, dy_max * dy_max);

    if (max_distance_sq <= radius_sq) {
        return area;
    }

    // Boundary rectangle: approximate at the leaf.
    if (depth == 0) {
        const double mid_x = 0.5 * (min_x + max_x);
        const double mid_y = 0.5 * (min_y + max_y);

        const double dx = mid_x - center_x;
        const double dy = mid_y - center_y;

        return (dx * dx) + (dy * dy) <= radius_sq ? area : 0.0;
    }

    const double mid_x = 0.5 * (min_x + max_x);
    const double mid_y = 0.5 * (min_y + max_y);

    return EstimateCoverageRecursive(center_x, center_y, radius_sq, min_x, mid_x, min_y, mid_y, depth - 1)
           + EstimateCoverageRecursive(center_x, center_y, radius_sq, mid_x, max_x, min_y, mid_y, depth - 1)
           + EstimateCoverageRecursive(center_x, center_y, radius_sq, min_x, mid_x, mid_y, max_y, depth - 1)
           + EstimateCoverageRecursive(center_x, center_y, radius_sq, mid_x, max_x, mid_y, max_y, depth - 1);
}

double HyperRGG2DPolicy::CellCoverage(const Center& center, const LPFloat radius, const Cell& cell) const {
    const CellBounds bounds = GetCellBounds(cell);

    const double coverage =
        EstimatedCircleRectCoverage(center.x, center.y, bounds.min_x, bounds.max_x, bounds.min_y, bounds.max_y, radius);

    return coverage;
}

SInt HyperRGG2DPolicy::AddWholeCell(const Cell& cell, std::vector<PinRange>& ranges) const {
    const auto stored = TryGetStoredCell(cell);
    if (!stored) {
        return 0;
    }

    if (stored->size > 0) {
        ranges.push_back({.begin = stored->offset, .end = stored->offset + stored->size});
    }

    return stored->size;
}

HyperRGG2DPolicy::CellBounds HyperRGG2DPolicy::GetCellBounds(const Cell& cell) const {
    const double h = 1.0 / static_cast<double>(gen_->SafeTotalCellsPerDim());

    return {
        .min_x = static_cast<double>(cell.global_cell_x) * h,
        .max_x = static_cast<double>((cell.global_cell_x + 1)) * h,
        .min_y = static_cast<double>(cell.global_cell_y) * h,
        .max_y = static_cast<double>((cell.global_cell_y + 1)) * h,
    };
}

#ifdef KAGEN_ENABLE_HIERARCHICAL_CELLS

HyperRGG2DPolicy::CellBounds HyperRGG2DPolicy::GetCellRegionBounds(const CellRegion& region) const {
    SInt chunk_x;
    SInt chunk_y;

    gen_->Decode(region.chunk_id, chunk_x, chunk_y);

    const SInt global_start_x = chunk_x * gen_->cells_per_dim_ + region.start_x;

    const SInt global_start_y = chunk_y * gen_->cells_per_dim_ + region.start_y;

    const double h = 1.0 / static_cast<double>(gen_->SafeTotalCellsPerDim());

    return {
        .min_x = static_cast<double>(global_start_x) * h,

        .max_x = static_cast<double>(global_start_x + region.columns) * h,

        .min_y = static_cast<double>(global_start_y) * h,

        .max_y = static_cast<double>(global_start_y + region.rows) * h,
    };
}

CellBallRelation
HyperRGG2DPolicy::ClassifyCellRegion(const Center& center, const LPFloat radius, const CellRegion& region) const {
    const CellBounds bounds = GetCellRegionBounds(region);

    const LPFloat radius_sq = radius * radius;

    const LPFloat closest_x = std::clamp(center.x, bounds.min_x, bounds.max_x);

    const LPFloat closest_y = std::clamp(center.y, bounds.min_y, bounds.max_y);

    const LPFloat dx_min = center.x - closest_x;

    const LPFloat dy_min = center.y - closest_y;

    if (dx_min * dx_min + dy_min * dy_min > radius_sq) {
        return CellBallRelation::OUTSIDE;
    }

    const LPFloat dx1 = center.x - bounds.min_x;

    const LPFloat dx2 = center.x - bounds.max_x;

    const LPFloat dy1 = center.y - bounds.min_y;

    const LPFloat dy2 = center.y - bounds.max_y;

    const LPFloat max_dist_sq = std::max(dx1 * dx1, dx2 * dx2) + std::max(dy1 * dy1, dy2 * dy2);

    if (max_dist_sq <= radius_sq) {
        return CellBallRelation::INSIDE;
    }

    return CellBallRelation::PARTIAL;
}

SInt HyperRGG2DPolicy::AddWholeCellRegion(const CellRegion& region, std::vector<PinRange>& ranges) const {
    gen_->GenerateCells(region.chunk_id);

    SInt total_added = 0;

    for (SInt local_y = region.start_y; local_y < region.start_y + region.rows; ++local_y) {
        bool have_vertices = false;

        SInt range_begin = 0;
        SInt range_end   = 0;

        for (SInt local_x = region.start_x; local_x < region.start_x + region.columns; ++local_x) {
            const SInt cell_id = gen_->EncodeCell(local_x, local_y);

            const SInt global_cell_id = gen_->ComputeGlobalCellId(region.chunk_id, cell_id);

            const auto it = gen_->cells_.find(global_cell_id);

            //
            // Empty cells are deliberately not stored by
            // GenerateCells().
            //
            if (it == gen_->cells_.end()) {
                continue;
            }

            const SInt size = std::get<0>(it->second);

            const SInt offset = std::get<4>(it->second);

            if (size <= 0) {
                continue;
            }

            if (!have_vertices) {
                range_begin = offset;

                have_vertices = true;
            }

            range_end = offset + size;

            total_added += size;
        }

        if (have_vertices) {
            ranges.push_back({
                .begin = range_begin,
                .end   = range_end,
            });
        }
    }

    return total_added;
}

void HyperRGG2DPolicy::TraverseCellHierarchy(
    const Center& center, const LPFloat radius, const CellRegion& region, std::vector<Cell>& partial_cells,
    std::vector<PinRange>& ranges, bool& has_inside_region) const {
    const CellBallRelation relation = ClassifyCellRegion(center, radius, region);

    if (relation == CellBallRelation::OUTSIDE) {
        return;
    }

    if (relation == CellBallRelation::INSIDE) {
        AddWholeCellRegion(region, ranges);

        has_inside_region = true;
        return;
    }

    //
    // PARTIAL 1×1 region:
    // hand the individual cell back to the existing builder.
    //
    if (region.columns == 1 && region.rows == 1) {
        SInt chunk_x;
        SInt chunk_y;

        gen_->Decode(region.chunk_id, chunk_x, chunk_y);

        const SSInt global_cell_x = chunk_x * gen_->cells_per_dim_ + region.start_x;

        const SSInt global_cell_y = chunk_y * gen_->cells_per_dim_ + region.start_y;

        if (auto cell = TryMakeCell(global_cell_x, global_cell_y)) {
            partial_cells.push_back(*cell);
        }

        return;
    }
    if (region.columns >= region.rows && region.columns > 1) {
        const SInt left_columns = (region.columns + 1) / 2;

        const SInt right_columns = region.columns - left_columns;

        TraverseCellHierarchy(
            center, radius,
            {
                .chunk_id = region.chunk_id,

                .start_x = region.start_x,

                .start_y = region.start_y,

                .columns = left_columns,

                .rows = region.rows,
            },
            partial_cells, ranges, has_inside_region);

        if (right_columns > 0) {
            TraverseCellHierarchy(
                center, radius,
                {
                    .chunk_id = region.chunk_id,

                    .start_x = region.start_x + left_columns,

                    .start_y = region.start_y,

                    .columns = right_columns,

                    .rows = region.rows,
                },
                partial_cells, ranges, has_inside_region);
        }

        return;
    }
    const SInt upper_rows = (region.rows + 1) / 2;

    const SInt lower_rows = region.rows - upper_rows;

    TraverseCellHierarchy(
        center, radius,
        {
            .chunk_id = region.chunk_id,

            .start_x = region.start_x,

            .start_y = region.start_y,

            .columns = region.columns,

            .rows = upper_rows,
        },
        partial_cells, ranges, has_inside_region);

    if (lower_rows > 0) {
        TraverseCellHierarchy(
            center, radius,
            {
                .chunk_id = region.chunk_id,

                .start_x = region.start_x,

                .start_y = region.start_y + upper_rows,

                .columns = region.columns,

                .rows = lower_rows,
            },
            partial_cells, ranges, has_inside_region);
    }
}

bool HyperRGG2DPolicy::HierarchicalCandidateCells(
    const Center& center, const LPFloat radius, std::vector<Cell>& partial_cells, std::vector<PinRange>& ranges) const {
    const SInt total_cells_per_dim = gen_->SafeTotalCellsPerDim();

    const LPFloat cell_size = 1.0 / static_cast<LPFloat>(total_cells_per_dim);

    SSInt min_cell_x = static_cast<SSInt>(std::floor((center.x - radius) / cell_size));

    SSInt max_cell_x = static_cast<SSInt>(std::floor((center.x + radius) / cell_size));

    SSInt min_cell_y = static_cast<SSInt>(std::floor((center.y - radius) / cell_size));

    SSInt max_cell_y = static_cast<SSInt>(std::floor((center.y + radius) / cell_size));

    min_cell_x = std::max<SSInt>(0, min_cell_x);

    min_cell_y = std::max<SSInt>(0, min_cell_y);

    max_cell_x = std::min<SSInt>(total_cells_per_dim - 1, max_cell_x);

    max_cell_y = std::min<SSInt>(total_cells_per_dim - 1, max_cell_y);

    if (min_cell_x > max_cell_x || min_cell_y > max_cell_y) {
        return false;
    }

    const SInt first_chunk_x = static_cast<SInt>(min_cell_x) / gen_->cells_per_dim_;

    const SInt last_chunk_x = static_cast<SInt>(max_cell_x) / gen_->cells_per_dim_;

    const SInt first_chunk_y = static_cast<SInt>(min_cell_y) / gen_->cells_per_dim_;

    const SInt last_chunk_y = static_cast<SInt>(max_cell_y) / gen_->cells_per_dim_;

    bool has_inside_region = false;

    for (SInt chunk_x = first_chunk_x; chunk_x <= last_chunk_x; ++chunk_x) {
        for (SInt chunk_y = first_chunk_y; chunk_y <= last_chunk_y; ++chunk_y) {
            const SInt chunk_id = gen_->Encode(chunk_x, chunk_y);

            TraverseCellHierarchy(
                center, radius,
                {
                    .chunk_id = chunk_id,

                    .start_x = 0,
                    .start_y = 0,

                    .columns = gen_->cells_per_dim_,

                    .rows = gen_->cells_per_dim_,
                },
                partial_cells, ranges, has_inside_region);
        }
    }
    std::sort(partial_cells.begin(), partial_cells.end(), [](const Cell& a, const Cell& b) {
        if (a.global_cell_x != b.global_cell_x) {
            return a.global_cell_x < b.global_cell_x;
        }

        return a.global_cell_y < b.global_cell_y;
    });

    return has_inside_region;
}

#endif

std::optional<HyperRGG2DPolicy::PartialCellSample>
HyperRGG2DPolicy::PreparePartialCellSample(const Cell& cell, const double coverage) const {
    const auto stored = TryGetStoredCell(cell);

    if (!stored || stored->size <= 0 || coverage <= 0.0) {
        return std::nullopt;
    }

    const double clamped_coverage = std::clamp(coverage, 0.0, 1.0);

    const SInt count = static_cast<SInt>(std::round(static_cast<double>(stored->size) * clamped_coverage));

    if (count <= 0) {
        return std::nullopt;
    }

    return PartialCellSample{
        .stored = *stored,
        .count  = count,
    };
}

SInt HyperRGG2DPolicy::AddPartialCellRange(
    const Center& /*center*/, const Cell& cell, const double coverage, std::vector<SInt>& /*pins*/,
    std::vector<PinRange>& ranges) const {
    const auto sample = PreparePartialCellSample(cell, coverage);

    if (!sample) {
        return 0;
    }

    auto sampled = getRandomPinRange(sample->stored.size, sample->count, sample->stored.offset, gen_->mersenne);

    ranges.push_back(sampled);

    return sample->count;
}

SInt HyperRGG2DPolicy::AddPartialCellFloyd(
    const Center& /*center*/, const Cell& cell, const double coverage, std::vector<SInt>& pins,
    std::vector<PinRange>& /*ranges*/) const {
    const auto sample = PreparePartialCellSample(cell, coverage);

    if (!sample) {
        return 0;
    }

    FloydSampleGeometricAppend(sample->stored.offset, sample->stored.size, sample->count, rng_, pins, floyd_scratch_);

    return sample->count;
}

SInt HyperRGG2DPolicy::AddPartialCellExact(
    const Center& center, LPFloat radius, const Cell& cell, std::vector<SInt>& pins) const {
    const auto& vertices = ExactVerticesByX(cell);

    const double cx = static_cast<double>(center.x);
    const double cy = static_cast<double>(center.y);
    const double r  = static_cast<double>(radius);
    const double r2 = r * r;

    const double min_x = cx - r;
    const double max_x = cx + r;

    auto first = std::lower_bound(
        vertices.begin(), vertices.end(), min_x, [](const Vertex& v, double x) { return std::get<0>(v) < x; });

    auto last =
        std::upper_bound(first, vertices.end(), max_x, [](double x, const Vertex& v) { return x < std::get<0>(v); });

    SInt vertex_counter = 0;

    for (auto it = first; it != last; ++it) {
        const auto& vertex = *it;

        const double dx = cx - static_cast<double>(std::get<0>(vertex));
        const double dy = cy - static_cast<double>(std::get<1>(vertex));

        if ((dx * dx) + (dy * dy) <= r2) {
            pins.push_back(std::get<2>(vertex));
            ++vertex_counter;
        }
    }

    gen_->AddExactDebugStats(1, static_cast<SInt>(std::distance(first, last)), vertex_counter);
    return vertex_counter;
}

void HyperRGG2DPolicy::EmitHyperedge(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges) {
    gen_->PushHyperedgeCompressed(pins, ranges);
}

const std::vector<Vertex>& HyperRGG2DPolicy::ExactVerticesByX(const Cell& cell) const {
    if (IsLocalCell(cell) && !gen_->config_.coordinates) {
        auto it = gen_->vertices_.find(cell.global_cell_id);

        if (it != gen_->vertices_.end()) {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
            RecordLocalExactAccess(cell.global_cell_id);
#endif

            auto [_, inserted] = local_vertices_sorted_by_x_.insert(cell.global_cell_id);

            if (inserted) {
                std::sort(it->second.begin(), it->second.end(), [](const Vertex& a, const Vertex& b) {
                    return std::get<0>(a) < std::get<0>(b);
                });
            }

            return it->second;
        }
    }

    // Remote cell, or local cell whose vertices have not been generated yet.
    return ExactRemoteCell(cell).vertices_by_x;
}

const HyperRGG2DPolicy::CachedExactCell& HyperRGG2DPolicy::ExactRemoteCell(const Cell& cell) const {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    RecordRemoteAccess(cell.global_cell_id);
#endif

    // Cache hit.
    auto it = exact_vertices_.find(cell.global_cell_id);
    if (it != exact_vertices_.end()) {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        ++exact_remote_cache_hits_;
#endif

        const auto lru_it = exact_lru_pos_.find(cell.global_cell_id);

        if (lru_it == exact_lru_pos_.end()) {
            throw std::logic_error("LRU cache metadata is inconsistent");
        }

        exact_lru_.splice(exact_lru_.begin(), exact_lru_, lru_it->second);

        lru_it->second = exact_lru_.begin();

        return it->second;
    }

    // Cache miss.
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    ++exact_remote_cache_misses_;
#endif

    if (exact_remote_cache_limit_ == 0) {
        throw std::logic_error("exact remote cache limit must be positive");
    }

    // Evict before inserting, so the cache never temporarily contains
    // exact_remote_cache_limit_ + 1 cells.
    if (exact_vertices_.size() >= exact_remote_cache_limit_) {
        if (exact_lru_.empty()) {
            throw std::logic_error("LRU cache metadata is inconsistent");
        }

        const SInt victim = exact_lru_.back();

        auto victim_it = exact_vertices_.find(victim);
        if (victim_it == exact_vertices_.end()) {
            throw std::logic_error("LRU cache metadata is inconsistent");
        }

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        exact_remote_live_vertices_ -= static_cast<SInt>(victim_it->second.vertices_by_x.capacity());
#endif

        exact_lru_.pop_back();
        exact_lru_pos_.erase(victim);
        exact_vertices_.erase(victim_it);
    }

    // Insert the requested cell.
    auto [inserted_it, inserted] = exact_vertices_.emplace(cell.global_cell_id, CachedExactCell{});

    if (!inserted) {
        throw std::logic_error("failed to insert exact remote cell into cache");
    }

    exact_lru_.push_front(cell.global_cell_id);
    exact_lru_pos_[cell.global_cell_id] = exact_lru_.begin();

    auto& cached = inserted_it->second;

    gen_->GenerateRemoteCellVertices(cell.chunk_id, cell.cell_id, cached.vertices_by_x);

    std::sort(cached.vertices_by_x.begin(), cached.vertices_by_x.end(), [](const Vertex& a, const Vertex& b) {
        return std::get<0>(a) < std::get<0>(b);
    });

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    const auto capacity = static_cast<SInt>(cached.vertices_by_x.capacity());

    exact_remote_live_vertices_ += capacity;
    exact_remote_peak_vertices_ = std::max(exact_remote_peak_vertices_, exact_remote_live_vertices_);

    exact_remote_cached_vertices_ += static_cast<SInt>(cached.vertices_by_x.size());
#endif

    return cached;
}

bool HyperRGG2DPolicy::IsLocalCell(const Cell& cell) const {
    return gen_->IsLocalChunk(cell.chunk_id);
}

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
void HyperRGG2DPolicy::RecordRemoteAccess(SInt global_cell_id) const {
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

void HyperRGG2DPolicy::PrintExactCacheStats() const {}

void HyperRGG2DPolicy::RecordLocalExactAccess(const SInt global_cell_id) const {
    const SInt t = ++local_exact_access_counter_;

    auto it = local_exact_last_access_.find(global_cell_id);

    if (it == local_exact_last_access_.end()) {
        local_exact_last_access_[global_cell_id] = t;
        return;
    }

    const SInt distance = t - it->second;
    it->second          = t;

    ++local_exact_reuse_count_;
    local_exact_reuse_distance_sum_ += distance;
    local_exact_reuse_distance_max_ = std::max(local_exact_reuse_distance_max_, distance);

    if (distance <= 1) {
        ++local_exact_reuse_distance_le_1_;
    } else if (distance <= 4) {
        ++local_exact_reuse_distance_le_4_;
    } else if (distance <= 16) {
        ++local_exact_reuse_distance_le_16_;
    } else {
        ++local_exact_reuse_distance_gt_16_;
    }
}
#endif

bool HyperRGG2DPolicy::ShouldApproximatePartialCell(const Cell& cell) const {
    //
    // Local cells should always use exact processing.
    //
    // If their vertices are already materialized, exact processing is
    // extremely cheap.
    //
    // Even if they are not materialized yet, our measurements show that
    // exact processing is still preferable.
    //
    if (IsLocalCell(cell)) {
        return false;
    }

    //
    // A remote cell already present in the exact-cell cache should also
    // use exact processing. Cache hits are extremely cheap.
    //
    if (exact_vertices_.contains(cell.global_cell_id)) {
        return false;
    }

    //
    // This is the only interesting approximation case:
    //
    //     remote cell
    //     + not currently present in the exact cache
    //
    // Exact processing would cause deterministic vertex regeneration.
    //
    return true;
}

} // namespace kagen
