#include "kagen/generators/hyper/h_geometric/hyper_rgg2d_policy.h"

#include "kagen/generators/geometric/geometric_2d.h"
#include "kagen/generators/hyper/h_geometric/hyper_rgg2d.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <tuple>
#include <vector>

namespace kagen {

namespace {

double CircleSegmentIntegral(double x, const double r) {
    x = std::clamp(x, -r, r);
    return 0.5 * ((x * std::sqrt(std::max(0.0, r * r - x * x))) + (r * r * std::asin(x / r)));
}

double IntegrateSqrt(const double a, const double b, const double r) {
    return CircleSegmentIntegral(b, r) - CircleSegmentIntegral(a, r);
}

double ExactCircleRectCoverage(
    const double center_x, const double center_y, const double min_x, const double max_x, const double min_y,
    const double max_y, const double radius) {
    const double cell_area = (max_x - min_x) * (max_y - min_y);

    if (cell_area <= 0.0 || radius <= 0.0) {
        return 0.0;
    }

    const double x0 = min_x - center_x;
    const double x1 = max_x - center_x;
    const double y0 = min_y - center_y;
    const double y1 = max_y - center_y;

    const double xa = std::max(x0, -radius);
    const double xb = std::min(x1, radius);

    if (xa >= xb) {
        return 0.0;
    }

    std::set<double> split_points;
    split_points.insert(xa);
    split_points.insert(xb);

    for (const double y: {y0, y1}) {
        if (std::abs(y) < radius) {
            const double x_cross = std::sqrt(radius * radius - y * y);

            const double left  = -x_cross;
            const double right = x_cross;

            if (left > xa && left < xb) {
                split_points.insert(left);
            }

            if (right > xa && right < xb) {
                split_points.insert(right);
            }
        }
    }

    const std::vector<double> xs(split_points.begin(), split_points.end());

    double area = 0.0;

    for (std::size_t i = 0; i + 1 < xs.size(); ++i) {
        const double a = xs[i];
        const double b = xs[i + 1];
        const double m = 0.5 * (a + b);

        const double circle_top    = std::sqrt(std::max(0.0, radius * radius - m * m));
        const double circle_bottom = -circle_top;

        const double top_inside    = std::min(y1, circle_top);
        const double bottom_inside = std::max(y0, circle_bottom);

        if (top_inside <= bottom_inside) {
            continue;
        }

        const bool top_is_circle    = circle_top <= y1;
        const bool bottom_is_circle = circle_bottom >= y0;

        if (top_is_circle && bottom_is_circle) {
            area += 2.0 * IntegrateSqrt(a, b, radius);
        } else if (top_is_circle && !bottom_is_circle) {
            area += IntegrateSqrt(a, b, radius) - y0 * (b - a);
        } else if (!top_is_circle && bottom_is_circle) {
            area += y1 * (b - a) + IntegrateSqrt(a, b, radius);
        } else {
            area += (y1 - y0) * (b - a);
        }
    }

    return std::clamp(area / cell_area, 0.0, 1.0);
}

} // namespace

HyperRGG2DPolicy::HyperRGG2DPolicy(HyperRGG2D& generator) : gen_(generator) {}

void HyperRGG2DPolicy::AddCenter(const Center&, std::vector<SInt>&) const {}

LPFloat HyperRGG2DPolicy::Radius(const Center& center) const {
    const SInt  global_cell_id = gen_.ComputeGlobalCellId(center.chunk_id, center.cell_id);
    const auto& cell           = gen_.cells_[global_cell_id];
    SInt        size           = std::get<0>(cell);
    double      lower_bound    = 0;
    if (size != 0) {
        SInt          cells_per_dim    = gen_.SafeTotalCellsPerDim();
        LPFloat       single_cell_area = LPFloat{1.0} / static_cast<LPFloat>(cells_per_dim * cells_per_dim);
        const LPFloat density          = single_cell_area > 0.0 ? static_cast<LPFloat>(size) / single_cell_area : 0.0;
        lower_bound                    = std::sqrt(2.0 / (M_PI * density));
    }
    return static_cast<LPFloat>(getSampledOrConstantRadius(gen_.config_, center.sampled_id, lower_bound));
}

std::vector<HyperRGG2DPolicy::Cell> HyperRGG2DPolicy::CandidateCells(const Center& center, const LPFloat radius) const {
    SInt center_cell_x;
    SInt center_cell_y;
    gen_.DecodeCell(center.cell_id, center_cell_x, center_cell_y);

    SInt center_chunk_x;
    SInt center_chunk_y;
    gen_.Decode(center.chunk_id, center_chunk_x, center_chunk_y);

    const SInt total_cells_per_dim = gen_.SafeTotalCellsPerDim();

    const SSInt global_center_cell_x = static_cast<SSInt>(center_chunk_x * gen_.cells_per_dim_ + center_cell_x);
    const SSInt global_center_cell_y = static_cast<SSInt>(center_chunk_y * gen_.cells_per_dim_ + center_cell_y);

    const SSInt cell_radius = static_cast<SSInt>(std::ceil(radius * static_cast<LPFloat>(total_cells_per_dim)));

    std::vector<Cell> cells;

    for (SSInt dx = -cell_radius; dx <= cell_radius; ++dx) {
        const SSInt max_dy =
            static_cast<SSInt>(std::floor(std::sqrt(static_cast<LPFloat>(cell_radius * cell_radius - dx * dx))));

        for (SSInt dy = -max_dy; dy <= max_dy; ++dy) {
            const SSInt global_cell_x = global_center_cell_x + dx;
            const SSInt global_cell_y = global_center_cell_y + dy;

            const SSInt total_cells_per_dim_signed = static_cast<SSInt>(total_cells_per_dim);

            if (global_cell_x < 0 || global_cell_y < 0 || global_cell_x >= total_cells_per_dim_signed
                || global_cell_y >= total_cells_per_dim_signed) {
                continue;
            }

            const SInt neighbor_chunk_x = static_cast<SInt>(global_cell_x) / gen_.cells_per_dim_;
            const SInt neighbor_chunk_y = static_cast<SInt>(global_cell_y) / gen_.cells_per_dim_;

            const SInt neighbor_cell_x = static_cast<SInt>(global_cell_x) % gen_.cells_per_dim_;
            const SInt neighbor_cell_y = static_cast<SInt>(global_cell_y) % gen_.cells_per_dim_;

            cells.push_back(
                Cell{
                    .chunk_id      = gen_.Encode(neighbor_chunk_x, neighbor_chunk_y),
                    .cell_id       = gen_.EncodeCell(neighbor_cell_x, neighbor_cell_y),
                    .global_cell_x = static_cast<SInt>(global_cell_x),
                    .global_cell_y = static_cast<SInt>(global_cell_y),
                });
        }
    }

    return cells;
}

CellBallRelation HyperRGG2DPolicy::ClassifyCell(const Center& center, const LPFloat radius, const Cell& cell) const {
    const SInt total_cells_per_dim = gen_.SafeTotalCellsPerDim();

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

double HyperRGG2DPolicy::CellCoverage(const Center& center, const LPFloat radius, const Cell& cell) const {
    const SInt total_cells_per_dim = gen_.SafeTotalCellsPerDim();

    const double min_x = static_cast<double>(cell.global_cell_x) / static_cast<double>(total_cells_per_dim);
    const double max_x = static_cast<double>(cell.global_cell_x + 1) / static_cast<double>(total_cells_per_dim);
    const double min_y = static_cast<double>(cell.global_cell_y) / static_cast<double>(total_cells_per_dim);
    const double max_y = static_cast<double>(cell.global_cell_y + 1) / static_cast<double>(total_cells_per_dim);

    return ExactCircleRectCoverage(center.x, center.y, min_x, max_x, min_y, max_y, radius);
}

void HyperRGG2DPolicy::AddWholeCell(const Cell& cell, std::vector<PinRange>& ranges) const {
    const SInt global_cell_id = gen_.ComputeGlobalCellId(cell.chunk_id, cell.cell_id);

    if (gen_.cells_.find(global_cell_id) == gen_.cells_.end()) {
        return;
    }

    const auto& stored_cell = gen_.cells_[global_cell_id];

    const SInt size  = std::get<0>(stored_cell);
    const SInt begin = std::get<4>(stored_cell);

    if (size > 0) {
        ranges.push_back({.begin = begin, .end = begin + size});
    }
}

void HyperRGG2DPolicy::AddPartialCell(
    const Center& center, const Cell& cell, const double coverage, std::vector<PinRange>& ranges) const {
    const SInt global_cell_id = gen_.ComputeGlobalCellId(cell.chunk_id, cell.cell_id);

    if (gen_.cells_.find(global_cell_id) == gen_.cells_.end()) {
        return;
    }

    const auto& target_cell = gen_.cells_[global_cell_id];

    const SInt size                     = std::get<0>(target_cell);
    const SInt offset                   = std::get<4>(target_cell);
    const SInt estimated_hyperedge_part = std::floor(size * coverage);
    const SInt seed                     = sampling::Spooky::hash(gen_.config_.seed + (131 * center.sampled_id));

    ranges.push_back(getRandomPinRange(size, estimated_hyperedge_part, offset, seed, gen_.config_));
}

void HyperRGG2DPolicy::AddPartialCellExact(
    const Center& center, const LPFloat radius, const Cell& cell, std::vector<SInt>& pins) const {
    std::vector<Vertex> vertices;
    gen_.GenerateVertices(cell.chunk_id, cell.cell_id, vertices);

    const double radius_sq = static_cast<double>(radius) * static_cast<double>(radius);

    for (const auto& vertex: vertices) {
        const double dx = static_cast<double>(center.x) - static_cast<double>(std::get<0>(vertex));
        const double dy = static_cast<double>(center.y) - static_cast<double>(std::get<1>(vertex));

        if ((dx * dx) + (dy * dy) <= radius_sq) {
            pins.push_back(std::get<2>(vertex));
        }
    }
}

void HyperRGG2DPolicy::EmitHyperedge(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges) {
    gen_.PushHyperedgeCompressed(pins, ranges);
}

} // namespace kagen
