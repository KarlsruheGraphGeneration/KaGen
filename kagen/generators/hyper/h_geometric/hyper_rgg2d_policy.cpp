#include "kagen/generators/hyper/h_geometric/hyper_rgg2d_policy.h"

#include "kagen/generators/hyper/h_geometric/hyper_rgg2d.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>

namespace kagen {

HyperRGG2DPolicy::HyperRGG2DPolicy(HyperRGG2D& generator) : gen_(&generator), rng_(RNGWrapper(generator.config_)) {}

void HyperRGG2DPolicy::AddCenter(const Center&, std::vector<SInt>&) const {}

SInt HyperRGG2DPolicy::GetNumVerticeOfCellCoord(const SSInt global_cell_x, const SSInt global_cell_y) {
    const auto total_cells_per_dim_signed = static_cast<SSInt>(gen_->SafeTotalCellsPerDim());
    if (global_cell_x < 0 || global_cell_y < 0 || global_cell_x >= total_cells_per_dim_signed
        || global_cell_y >= total_cells_per_dim_signed) {
        return SInt{0};
    }

    const SInt chunk_x = static_cast<SInt>(global_cell_x) / gen_->cells_per_dim_;
    const SInt chunk_y = static_cast<SInt>(global_cell_y) / gen_->cells_per_dim_;
    const SInt cell_x  = static_cast<SInt>(global_cell_x) % gen_->cells_per_dim_;
    const SInt cell_y  = static_cast<SInt>(global_cell_y) % gen_->cells_per_dim_;

    const SInt global_cell_id =
        gen_->ComputeGlobalCellId(gen_->Encode(chunk_x, chunk_y), gen_->EncodeCell(cell_x, cell_y));

    const auto cell_it = gen_->cells_.find(global_cell_id);
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
    const double     lower_bound = const_cast<HyperRGG2DPolicy*>(this)->MinimumRadius(center);
    constexpr double upper_bound = 1.0;
    const double     actualRadius =
        getSampledOrConstantRadius(gen_->config_, center.sampled_id, lower_bound, upper_bound, rng_);
    return static_cast<LPFloat>(actualRadius);
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

            const SSInt total_cells_per_dim_signed = static_cast<SSInt>(total_cells_per_dim);

            if (global_cell_x < 0 || global_cell_y < 0 || global_cell_x >= total_cells_per_dim_signed
                || global_cell_y >= total_cells_per_dim_signed) {
                continue;
            }

            const SInt neighbor_chunk_x = static_cast<SInt>(global_cell_x) / gen_->cells_per_dim_;
            const SInt neighbor_chunk_y = static_cast<SInt>(global_cell_y) / gen_->cells_per_dim_;

            const SInt neighbor_cell_x = static_cast<SInt>(global_cell_x) % gen_->cells_per_dim_;
            const SInt neighbor_cell_y = static_cast<SInt>(global_cell_y) % gen_->cells_per_dim_;

            const SInt chunk_id = gen_->Encode(neighbor_chunk_x, neighbor_chunk_y);
            const SInt cell_id  = gen_->EncodeCell(neighbor_cell_x, neighbor_cell_y);

            cells.push_back(
                Cell{
                    .chunk_id      = chunk_id,
                    .cell_id       = cell_id,
                    .global_cell_x = static_cast<SInt>(global_cell_x),
                    .global_cell_y = static_cast<SInt>(global_cell_y),

                    .global_cell_id = gen_->ComputeGlobalCellId(chunk_id, cell_id),
                });
        }
    }
}

CellBallRelation HyperRGG2DPolicy::ClassifyCell(const Center& center, const LPFloat radius, const Cell& cell) const {
    const SInt total_cells_per_dim = gen_->SafeTotalCellsPerDim();

    const LPFloat min_x = static_cast<LPFloat>(cell.global_cell_x) / total_cells_per_dim;
    const LPFloat max_x = static_cast<LPFloat>(cell.global_cell_x + 1) / total_cells_per_dim;

    const LPFloat min_y = static_cast<LPFloat>(cell.global_cell_y) / total_cells_per_dim;
    const LPFloat max_y = static_cast<LPFloat>(cell.global_cell_y + 1) / total_cells_per_dim;

    const LPFloat radius_sq = radius * radius;

    const LPFloat closest_x = std::clamp(center.x, min_x, max_x);
    const LPFloat closest_y = std::clamp(center.y, min_y, max_y);

    const LPFloat dx_min = center.x - closest_x;
    const LPFloat dy_min = center.y - closest_y;

    if (dx_min * dx_min + dy_min * dy_min > radius_sq) {
        return CellBallRelation::OUTSIDE;
    }

    const LPFloat dx1 = center.x - min_x;
    const LPFloat dx2 = center.x - max_x;
    const LPFloat dy1 = center.y - min_y;
    const LPFloat dy2 = center.y - max_y;

    const LPFloat max_dist_sq = std::max(dx1 * dx1, dx2 * dx2) + std::max(dy1 * dy1, dy2 * dy2);

    if (max_dist_sq <= radius_sq) {
        return CellBallRelation::INSIDE;
    }

    return CellBallRelation::PARTIAL;
}
double HyperRGG2DPolicy::EstimatedCellCoverage(const Center& center, const LPFloat radius, const Cell& cell) const {
    const SInt total_cells_per_dim = gen_->SafeTotalCellsPerDim();

    const double h = 1.0 / static_cast<double>(total_cells_per_dim);

    const double min_x = static_cast<double>(cell.global_cell_x) * h;
    const double min_y = static_cast<double>(cell.global_cell_y) * h;

    const double cx = static_cast<double>(center.x);
    const double cy = static_cast<double>(center.y);
    const double r2 = static_cast<double>(radius) * static_cast<double>(radius);

    auto inside = [&](const double x, const double y) {
        const double dx = x - cx;
        const double dy = y - cy;
        return dx * dx + dy * dy <= r2;
    };

    int hits = 0;

    // 3x3 stratified sample points
    for (int ix = 0; ix < 3; ++ix) {
        for (int iy = 0; iy < 3; ++iy) {
            const double x = min_x + (static_cast<double>(ix) + 0.5) * h / 3.0;
            const double y = min_y + (static_cast<double>(iy) + 0.5) * h / 3.0;

            hits += inside(x, y) ? 1 : 0;
        }
    }

    return static_cast<double>(hits) / 9.0;
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
    const SInt   total_cells_per_dim = gen_->SafeTotalCellsPerDim();
    const double min_x = static_cast<double>(cell.global_cell_x) / static_cast<double>(total_cells_per_dim);
    const double max_x = static_cast<double>(cell.global_cell_x + 1) / static_cast<double>(total_cells_per_dim);
    const double min_y = static_cast<double>(cell.global_cell_y) / static_cast<double>(total_cells_per_dim);
    const double max_y = static_cast<double>(cell.global_cell_y + 1) / static_cast<double>(total_cells_per_dim);

    return EstimatedCircleRectCoverage(center.x, center.y, min_x, max_x, min_y, max_y, radius);
    /*
    return ExactCircleRectCoverage(center.x, center.y, min_x, max_x, min_y, max_y, radius); Currently debugged for
     faster estimate.
    */
}

SInt HyperRGG2DPolicy::AddWholeCell(const Cell& cell, std::vector<PinRange>& ranges) const {
    const auto it = gen_->cells_.find(cell.global_cell_id);
    if (it == gen_->cells_.end()) {
        return 0;
    }

    const auto& stored_cell = it->second;

    const SInt size  = std::get<0>(stored_cell);
    const SInt begin = std::get<4>(stored_cell);

    if (size > 0) {
        ranges.push_back({.begin = begin, .end = begin + size});
    }
    return size;
}

SInt HyperRGG2DPolicy::AddPartialCellRange(
    const Center& center, const Cell& cell, const double coverage, std::vector<SInt>& pins,
    std::vector<PinRange>& ranges) const {
    const auto it = gen_->cells_.find(cell.global_cell_id);
    if (it == gen_->cells_.end()) {
        return 0;
    }

    const auto& target_cell = it->second;

    const SInt size   = std::get<0>(target_cell);
    const SInt offset = std::get<4>(target_cell);
    const SInt k      = std::floor(static_cast<double>(size) * coverage);

    if (k <= 0) {
        return 0;
    }

    SInt seed = sampling::Spooky::hash(gen_->config_.seed + (131 * center.sampled_id) + (9973 * cell.global_cell_id));

    // auto sampled = FloydSample(offset, size, k, rng_, seed);
    auto sampled = getRandomPinRange(size, k, offset, seed, gen_->mersenne);
    ranges.insert(ranges.end(), sampled);

    return k;
}

SInt HyperRGG2DPolicy::AddPartialCellFloyd(
    const Center& center, const Cell& cell, const double coverage, std::vector<SInt>& pins,
    std::vector<PinRange>& ranges) const {
    const auto it = gen_->cells_.find(cell.global_cell_id);
    if (it == gen_->cells_.end()) {
        return 0;
    }

    const auto& target_cell = it->second;

    const SInt size   = std::get<0>(target_cell);
    const SInt offset = std::get<4>(target_cell);
    const SInt k      = std::floor(static_cast<double>(size) * coverage);

    if (k <= 0) {
        return 0;
    }

    SInt seed = sampling::Spooky::hash(gen_->config_.seed + (131 * center.sampled_id) + (9973 * cell.global_cell_id));

    FloydSampleGeometricAppend(offset, size, k, rng_, seed, pins, floyd_scratch_);

    return k;
}

SInt HyperRGG2DPolicy::AddPartialCellExact(
    const Center& center, LPFloat radius, const Cell& cell, std::vector<SInt>& pins) const {
    const auto& cached   = ExactCell(cell);
    const auto& vertices = cached.vertices_by_x;

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

} // namespace kagen
