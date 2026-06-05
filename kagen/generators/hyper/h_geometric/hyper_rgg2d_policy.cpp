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

namespace {

double CircleSegmentIntegral(double x_value, const double radius) {
    x_value = std::clamp(x_value, -radius, radius);
    return 0.5
           * ((x_value * std::sqrt(std::max(0.0, (radius * radius) - (x_value * x_value))))
              + (radius * radius * std::asin(x_value / radius)));
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

    double split_points[6];
    int    count = 0;

    auto add_split = [&](double x) {
        if (x <= xa || x >= xb) {
            return;
        }
        for (int i = 0; i < count; ++i) {
            if (std::abs(split_points[i] - x) < 1e-15) {
                return;
            }
        }
        split_points[count++] = x;
    };
    split_points[count++]  = xa;
    split_points[count++]  = xb;
    const double radius_sq = radius * radius;

    for (double y: {y0, y1}) {
        if (std::abs(y) < radius) {
            const double x_cross = std::sqrt(radius_sq - (y * y));

            add_split(-x_cross);
            add_split(x_cross);
        }
    }
    std::sort(split_points, split_points + count);

    double area = 0.0;

    for (int i = 0; i + 1 < count; ++i) {
        const double a = split_points[i];
        const double b = split_points[i + 1];
        const double m = 0.5 * (a + b);

        const double circle_top    = std::sqrt(std::max(0.0, radius_sq - m * m));
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

double HyperRGG2DPolicy::MinimumRadius(const Center& center) {
    constexpr double target_expected_vertices = 2.0;
    constexpr double max_radius               = 1.0;

    const SInt total_cells_per_dim = gen_->SafeTotalCellsPerDim();
    if (total_cells_per_dim <= 0) {
        return 0.0;
    }

    const auto   total_cells_per_dim_signed = static_cast<SSInt>(total_cells_per_dim);
    const double cell_size                  = 1.0 / static_cast<double>(total_cells_per_dim);

    struct CachedCell {
        double min_x;
        double max_x;
        double min_y;
        double max_y;
        SInt   size;
    };

    std::vector<CachedCell> cached_cells;

    const double min_x = std::max(0.0, static_cast<double>(center.x) - max_radius);
    const double max_x = std::min(1.0, static_cast<double>(center.x) + max_radius);
    const double min_y = std::max(0.0, static_cast<double>(center.y) - max_radius);
    const double max_y = std::min(1.0, static_cast<double>(center.y) + max_radius);

    const SSInt first_x = std::max<SSInt>(0, static_cast<SSInt>(std::floor(min_x * total_cells_per_dim)));
    const SSInt last_x  = std::min<SSInt>(
        total_cells_per_dim_signed - 1,
        static_cast<SSInt>(std::floor(std::nextafter(max_x, 0.0) * total_cells_per_dim)));

    const SSInt first_y = std::max<SSInt>(0, static_cast<SSInt>(std::floor(min_y * total_cells_per_dim)));
    const SSInt last_y  = std::min<SSInt>(
        total_cells_per_dim_signed - 1,
        static_cast<SSInt>(std::floor(std::nextafter(max_y, 0.0) * total_cells_per_dim)));

    for (SSInt global_cell_x = first_x; global_cell_x <= last_x; ++global_cell_x) {
        const double rect_min_x = static_cast<double>(global_cell_x) * cell_size;
        const double rect_max_x = static_cast<double>(global_cell_x + 1) * cell_size;

        for (SSInt global_cell_y = first_y; global_cell_y <= last_y; ++global_cell_y) {
            const SInt size = GetNumVerticeOfCellCoord(global_cell_x, global_cell_y);
            if (size == 0) {
                continue;
            }

            cached_cells.push_back({
                .min_x = rect_min_x,
                .max_x = rect_max_x,
                .min_y = static_cast<double>(global_cell_y) * cell_size,
                .max_y = static_cast<double>(global_cell_y + 1) * cell_size,
                .size  = size,
            });
        }
    }

    auto expected_vertices = [&](const double radius) {
        if (radius <= 0.0) {
            return 0.0;
        }

        double expected = 0.0;

        for (const CachedCell& cell: cached_cells) {
            const double coverage = EstimatedCircleRectCoverage(
                static_cast<double>(center.x), static_cast<double>(center.y), cell.min_x, cell.max_x, cell.min_y,
                cell.max_y, radius);
            ;

            expected += static_cast<double>(cell.size) * coverage;

            if (expected >= target_expected_vertices) {
                return expected;
            }
        }

        return expected;
    };

    if (expected_vertices(max_radius) < target_expected_vertices) {
        return max_radius;
    }

    double lower = 0.0;
    double upper = cell_size;

    while (upper < max_radius && expected_vertices(upper) < target_expected_vertices) {
        upper = std::min(max_radius, 2.0 * upper);
    }

    for (int iteration = 0; iteration < 32; ++iteration) {
        const double mid = 0.5 * (lower + upper);

        if (expected_vertices(mid) >= target_expected_vertices) {
            upper = mid;
        } else {
            lower = mid;
        }
    }

    return upper;
}
LPFloat HyperRGG2DPolicy::Radius(const Center& center) const {
    const double     lower_bound = const_cast<HyperRGG2DPolicy*>(this)->MinimumRadius(center);
    constexpr double upper_bound = 1.0;
    const double actualRadius = getSampledOrConstantRadius(gen_->config_, center.sampled_id, lower_bound, upper_bound);
    return static_cast<LPFloat>(actualRadius);
}

std::vector<HyperRGG2DPolicy::Cell> HyperRGG2DPolicy::CandidateCells(const Center& center, const LPFloat radius) const {
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

    std::vector<Cell> cells;

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

            cells.push_back(
                Cell{
                    .chunk_id      = gen_->Encode(neighbor_chunk_x, neighbor_chunk_y),
                    .cell_id       = gen_->EncodeCell(neighbor_cell_x, neighbor_cell_y),
                    .global_cell_x = static_cast<SInt>(global_cell_x),
                    .global_cell_y = static_cast<SInt>(global_cell_y),
                });
        }
    }

    return cells;
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
            const double x = min_x + (ix + 0.5) * dx;
            const double y = min_y + (iy + 0.5) * dy;
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
    const SInt global_cell_id = gen_->ComputeGlobalCellId(cell.chunk_id, cell.cell_id);

    if (gen_->cells_.find(global_cell_id) == gen_->cells_.end()) {
        return 0;
    }

    const auto& stored_cell = gen_->cells_[global_cell_id];

    const SInt size  = std::get<0>(stored_cell);
    const SInt begin = std::get<4>(stored_cell);

    if (size > 0) {
        ranges.push_back({.begin = begin, .end = begin + size});
    }
    return size;
}

SInt HyperRGG2DPolicy::AddPartialCell(
    const Center& center, const Cell& cell, const double coverage, std::vector<PinRange>& ranges) const {
    const SInt global_cell_id = gen_->ComputeGlobalCellId(cell.chunk_id, cell.cell_id);

    if (gen_->cells_.find(global_cell_id) == gen_->cells_.end()) {
        return 0;
    }

    const auto& target_cell = gen_->cells_[global_cell_id];

    const SInt size                     = std::get<0>(target_cell);
    const SInt offset                   = std::get<4>(target_cell);
    const SInt estimated_hyperedge_part = std::floor(size * coverage);
    const SInt seed = sampling::Spooky::hash(gen_->config_.seed + 131 * center.sampled_id + (9973 * global_cell_id));

    ranges.push_back(getRandomPinRange(size, estimated_hyperedge_part, offset, seed, gen_->mersenne));
    return estimated_hyperedge_part;
}

SInt HyperRGG2DPolicy::AddPartialCellExact(
    const Center& center, const LPFloat radius, const Cell& cell, std::vector<SInt>& pins) const {
    std::vector<Vertex> vertices;
    gen_->GenerateVertices(cell.chunk_id, cell.cell_id, vertices);

    const double radius_sq      = static_cast<double>(radius) * static_cast<double>(radius);
    SInt         vertex_counter = 0;
    for (const auto& vertex: vertices) {
        const double dx = static_cast<double>(center.x) - static_cast<double>(std::get<0>(vertex));
        const double dy = static_cast<double>(center.y) - static_cast<double>(std::get<1>(vertex));

        if ((dx * dx) + (dy * dy) <= radius_sq) {
            pins.push_back(std::get<2>(vertex));
            vertex_counter++;
        }
    }
    return vertex_counter;
}

void HyperRGG2DPolicy::EmitHyperedge(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges) {
    gen_->PushHyperedgeCompressed(pins, ranges);
}

} // namespace kagen
