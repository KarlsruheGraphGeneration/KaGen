#include "kagen/generators/hyper/h_geometric/hyper_rgg2d_policy.h"

#include "kagen/generators/hyper/h_geometric/hyper_rgg2d.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <optional>
#include <vector>

namespace kagen {

HyperRGG2DPolicy::HyperRGG2DPolicy(HyperRGG2D& generator) : gen_(&generator), rng_(RNGWrapper(generator.config_)) {}

void HyperRGG2DPolicy::AddCenter(const Center&, std::vector<SInt>&) const {}

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

double HyperRGG2DPolicy::MinimumRadius(const Center& /*unused*/) {
    constexpr double target = 2.0;

    const double density = static_cast<double>(gen_->config_.n);

    return std::sqrt(target / (M_PI * density));
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
    SInt center_cell_x;
    SInt center_cell_y;
    gen_->DecodeCell(center.cell_id, center_cell_x, center_cell_y);

    SInt center_chunk_x;
    SInt center_chunk_y;
    gen_->Decode(center.chunk_id, center_chunk_x, center_chunk_y);

    const SInt total_cells_per_dim = gen_->SafeTotalCellsPerDim();

    const auto global_center_cell_x = static_cast<SSInt>((center_chunk_x * gen_->cells_per_dim_) + center_cell_x);
    const auto global_center_cell_y = static_cast<SSInt>((center_chunk_y * gen_->cells_per_dim_) + center_cell_y);

    const auto cell_radius = static_cast<SSInt>(std::ceil(radius * static_cast<LPFloat>(total_cells_per_dim)));

    for (SSInt dx = -cell_radius; dx <= cell_radius; ++dx) {
        const auto max_dy =
            static_cast<SSInt>(std::floor(std::sqrt(static_cast<LPFloat>((cell_radius * cell_radius) - dx * dx))));

        for (SSInt dy = -max_dy; dy <= max_dy; ++dy) {
            const SSInt global_cell_x = global_center_cell_x + dx;
            const SSInt global_cell_y = global_center_cell_y + dy;

            if (auto cell = TryMakeCell(global_cell_x, global_cell_y)) {
                cells.push_back(*cell);
            }
        }
    }
}

std::optional<HyperRGG2DPolicy::StoredCell> HyperRGG2DPolicy::TryGetStoredCell(const Cell& cell) const {
    const auto it = gen_->cells_.find(cell.global_cell_id);
    if (it == gen_->cells_.end()) {
        return std::nullopt;
    }

    return StoredCell{
        .size   = std::get<0>(it->second),
        .offset = std::get<4>(it->second),
    };
}

CellBallRelation HyperRGG2DPolicy::ClassifyCell(const Center& center, const LPFloat radius, const Cell& cell) const {
    const CellBounds bounds    = GetCellBounds(cell);
    const LPFloat    radius_sq = radius * radius;

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
    const double r2 = radius * radius;

    auto inside = [&](const double x, const double y) {
        const double dx = x - center_x;
        const double dy = y - center_y;
        return dx * dx + dy * dy <= r2;
    };

    const double dx = (max_x - min_x) / 3.0;
    const double dy = (max_y - min_y) / 3.0;

    int hits = 0;

    for (int ix = 0; ix < 3; ++ix) {
        for (int iy = 0; iy < 3; ++iy) {
            const double x = min_x + ((ix + 0.5) * dx);
            const double y = min_y + ((iy + 0.5) * dy);
            hits += inside(x, y) ? 1 : 0;
        }
    }

    return static_cast<double>(hits) / 9.0;
}

double HyperRGG2DPolicy::CellCoverage(const Center& center, const LPFloat radius, const Cell& cell) const {
    const CellBounds bounds = GetCellBounds(cell);

    return EstimatedCircleRectCoverage(
        center.x, center.y, bounds.min_x, bounds.max_x, bounds.min_y, bounds.max_y, radius);
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

std::optional<HyperRGG2DPolicy::PartialCellSample>
HyperRGG2DPolicy::PreparePartialCellSample(const Center& center, const Cell& cell, double coverage) const {
    PartialCellSample sample;
    const auto        stored = TryGetStoredCell(cell);
    if (!stored) {
        return std::nullopt;
    }
    sample.stored = *stored;
    const SInt k  = std::floor(static_cast<double>(stored->size) * coverage);
    if (k <= 0) {
        return std::nullopt;
    }
    sample.count = k;

    sample.seed = sampling::Spooky::hash(gen_->config_.seed + (131 * center.sampled_id) + (9973 * cell.global_cell_id));

    return sample;
}

SInt HyperRGG2DPolicy::AddPartialCellRange(
    const Center& center, const Cell& cell, const double coverage, std::vector<SInt>& /*pins*/,
    std::vector<PinRange>& ranges) const {
    const auto sample = PreparePartialCellSample(center, cell, coverage);
    if (!sample) {
        return 0;
    }
    auto sampled =
        getRandomPinRange(sample->stored.size, sample->count, sample->stored.offset, sample->seed, gen_->mersenne);
    ranges.insert(ranges.end(), sampled);

    return sample->count;
}

SInt HyperRGG2DPolicy::AddPartialCellFloyd(
    const Center& center, const Cell& cell, const double coverage, std::vector<SInt>& pins,
    std::vector<PinRange>& /*ranges*/) const {
    const auto sample = PreparePartialCellSample(center, cell, coverage);
    if (!sample) {
        return 0;
    }

    SInt seed = sample->seed;

    FloydSampleGeometricAppend(
        sample->stored.offset, sample->stored.size, sample->count, rng_, seed, pins, floyd_scratch_);

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
        if (it == gen_->vertices_.end()) {
            static const std::vector<Vertex> empty;
            return empty; // or assert/throw, because this should not happen
        }

        auto [_, inserted] = local_vertices_sorted_by_x_.insert(cell.global_cell_id);
        if (inserted) {
            std::sort(it->second.begin(), it->second.end(), [](const Vertex& a, const Vertex& b) {
                return std::get<0>(a) < std::get<0>(b);
            });
        }

        return it->second;
    }

    return ExactRemoteCell(cell).vertices_by_x;
}
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

const HyperRGG2DPolicy::CachedExactCell&
HyperRGG2DPolicy::ExactRemoteCell(const Cell& cell) const {
    RecordRemoteAccess(cell.global_cell_id);

    auto it = exact_vertices_.find(cell.global_cell_id);

    if (it != exact_vertices_.end()) {
        ++exact_remote_cache_hits_;

        const auto lru_it = exact_lru_pos_.find(cell.global_cell_id);
        if (lru_it != exact_lru_pos_.end()) {
            exact_lru_.splice(
                exact_lru_.begin(),
                exact_lru_,
                lru_it->second
            );

            lru_it->second = exact_lru_.begin();
        }

        return it->second;
    }

    ++exact_remote_cache_misses_;

    auto [inserted_it, inserted] =
        exact_vertices_.emplace(
            cell.global_cell_id,
            CachedExactCell{}
        );

    if (!inserted) {
        throw std::logic_error(
            "failed to insert exact remote cell into cache"
        );
    }

    exact_lru_.push_front(cell.global_cell_id);
    exact_lru_pos_[cell.global_cell_id] = exact_lru_.begin();

    auto& cached = inserted_it->second;

    gen_->GenerateVertices(
        cell.chunk_id,
        cell.cell_id,
        cached.vertices_by_x
    );

    std::sort(
        cached.vertices_by_x.begin(),
        cached.vertices_by_x.end(),
        [](const Vertex& a, const Vertex& b) {
            return std::get<0>(a) < std::get<0>(b);
        }
    );

    exact_remote_cached_vertices_ +=
        static_cast<SInt>(cached.vertices_by_x.size());

    while (exact_vertices_.size() > exact_remote_cache_limit_) {
        const SInt victim = exact_lru_.back();

        exact_lru_.pop_back();
        exact_lru_pos_.erase(victim);

        // Do not evict the newly inserted cell before returning it.
        if (victim == cell.global_cell_id) {
            throw std::logic_error(
                "exact remote cache limit is too small"
            );
        }

        exact_vertices_.erase(victim);
    }

    return cached;
}

void HyperRGG2DPolicy::PrintExactCacheStats() const {
    std::cerr << " exact_remote_hits=" << exact_remote_cache_hits_
              << " exact_remote_misses=" << exact_remote_cache_misses_
              << " exact_remote_cached_vertices=" << exact_remote_cached_vertices_ << '\n';
    const double avg_reuse = exact_remote_reuse_count_ ? static_cast<double>(exact_remote_reuse_distance_sum_)
                                                             / static_cast<double>(exact_remote_reuse_count_)
                                                       : 0.0;

    std::cerr << " remote_reuse_count=" << exact_remote_reuse_count_ << " remote_reuse_avg=" << avg_reuse
              << " remote_reuse_max=" << exact_remote_reuse_distance_max_
              << " reuse<=1=" << exact_remote_reuse_distance_le_1_ << " reuse<=4=" << exact_remote_reuse_distance_le_4_
              << " reuse<=16=" << exact_remote_reuse_distance_le_16_
              << " reuse>16=" << exact_remote_reuse_distance_gt_16_ << '\n';
}

bool HyperRGG2DPolicy::IsLocalCell(const Cell& cell) const {
    return gen_->IsLocalChunk(cell.chunk_id);
}

} // namespace kagen
