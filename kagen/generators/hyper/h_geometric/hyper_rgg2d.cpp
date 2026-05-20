#include "kagen/generators/hyper/h_geometric/hyper_rgg2d.h"

#include "kagen/context.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/sampling/hash.hpp"
#include "kagen/tools/geometry.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace kagen {

HyperRGG2D::HyperRGG2D(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : SpatialGrid2D(config, rank, size) {}

bool HyperRGG2D::GenerateBallHyperedge(
    const LPFloat center_x, const LPFloat center_y, const SInt sampled_center_id, const SInt center_chunk_id,
    const SInt center_cell_id) {
    std::vector<SInt>     pins;
    std::vector<PinRange> ranges;

    SInt center_cell_x;
    SInt center_cell_y;
    DecodeCell(center_cell_id, center_cell_x, center_cell_y);

    SInt center_chunk_x;
    SInt center_chunk_y;
    Decode(center_chunk_id, center_chunk_x, center_chunk_y);

    const LPFloat radius    = static_cast<LPFloat>(getSampledOrConstantRadius(config_, sampled_center_id));
    const LPFloat radius_sq = radius * radius;

    const SInt total_cells_per_dim = SafeTotalCellsPerDim();

    const SSInt global_center_cell_x = static_cast<SSInt>(center_chunk_x * cells_per_dim_ + center_cell_x);
    const SSInt global_center_cell_y = static_cast<SSInt>(center_chunk_y * cells_per_dim_ + center_cell_y);

    const SSInt cell_radius = static_cast<SSInt>(std::ceil(radius * static_cast<LPFloat>(total_cells_per_dim)));

    for (SSInt dx = -cell_radius; dx <= cell_radius; ++dx) {
        const SSInt max_dy =
            static_cast<SSInt>(std::floor(std::sqrt(static_cast<LPFloat>(cell_radius * cell_radius - dx * dx))));

        for (SSInt dy = -max_dy; dy <= max_dy; ++dy) {
            const SSInt global_cell_x = global_center_cell_x + dx;
            const SSInt global_cell_y = global_center_cell_y + dy;

            if (global_cell_x < 0 || global_cell_y < 0 || global_cell_x >= total_cells_per_dim
                || global_cell_y >= total_cells_per_dim) {
                continue;
            }

            const SInt neighbor_chunk_x = global_cell_x / cells_per_dim_;
            const SInt neighbor_chunk_y = global_cell_y / cells_per_dim_;

            const SInt neighbor_cell_x = global_cell_x % cells_per_dim_;
            const SInt neighbor_cell_y = global_cell_y % cells_per_dim_;

            const SInt neighbor_chunk_id = Encode(neighbor_chunk_x, neighbor_chunk_y);
            const SInt neighbor_cell_id  = EncodeCell(neighbor_cell_x, neighbor_cell_y);

            const auto relation =
                ClassifyCellAgainstBall(center_x, center_y, radius, neighbor_chunk_id, neighbor_cell_id);

            if (relation == CellBallRelation::OUTSIDE) {
                continue;
            }

            if (relation == CellBallRelation::INSIDE) {
                AddWholeCellRange(neighbor_chunk_id, neighbor_cell_id, ranges);
                continue;
            }

            std::vector<Vertex> candidates;
            GenerateVertices(neighbor_chunk_id, neighbor_cell_id, candidates);

            for (const Vertex& candidate: candidates) {
                const LPFloat diff_x           = center_x - std::get<0>(candidate);
                const LPFloat diff_y           = center_y - std::get<1>(candidate);
                const LPFloat squared_distance = diff_x * diff_x + diff_y * diff_y;

                if (squared_distance <= radius_sq) {
                    pins.push_back(std::get<2>(candidate));
                }
            }
        }
    }

    auto normalized = NormalizeHyperedge(std::move(pins), std::move(ranges));

    pins   = std::move(normalized.first);
    ranges = std::move(normalized.second);

    SInt range_pin_count = 0;
    for (const PinRange& range: ranges) {
        range_pin_count += range.end - range.begin;
    }

    if (pins.size() + range_pin_count >= 2) {
        PushHyperedgeCompressed(pins, ranges);
        return true;
    }
    return false;
}
std::pair<std::vector<SInt>, std::vector<PinRange>>
HyperRGG2D::NormalizeHyperedge(std::vector<SInt> pins, std::vector<PinRange> ranges) {
    std::sort(pins.begin(), pins.end());
    pins.erase(std::unique(pins.begin(), pins.end()), pins.end());

    std::sort(ranges.begin(), ranges.end(), [](const PinRange& a, const PinRange& b) { return a.begin < b.begin; });

    std::vector<PinRange> normalized_ranges;

    for (const PinRange& range: ranges) {
        if (range.begin >= range.end) {
            continue;
        }

        if (!normalized_ranges.empty() && normalized_ranges.back().end >= range.begin) {
            normalized_ranges.back().end = std::max(normalized_ranges.back().end, range.end);
        } else {
            normalized_ranges.push_back(range);
        }
    }

    std::vector<SInt> normalized_pins;
    normalized_pins.reserve(pins.size());

    for (const SInt pin: pins) {
        bool covered_by_range = false;

        for (const PinRange& range: normalized_ranges) {
            if (pin < range.begin) {
                break;
            }

            if (pin >= range.begin && pin < range.end) {
                covered_by_range = true;
                break;
            }
        }

        if (!covered_by_range) {
            normalized_pins.push_back(pin);
        }
    }

    return {std::move(normalized_pins), std::move(normalized_ranges)};
}

CellBallRelation HyperRGG2D::ClassifyCellAgainstBall(
    const LPFloat center_x, const LPFloat center_y, const LPFloat radius, const SInt chunk_id, const SInt cell_id) {
    SInt chunk_x;
    SInt chunk_y;
    Decode(chunk_id, chunk_x, chunk_y);

    SInt cell_x;
    SInt cell_y;
    DecodeCell(cell_id, cell_x, cell_y);

    const auto total_cells_per_dim = static_cast<LPFloat>(SafeTotalCellsPerDim());

    const auto min_x = static_cast<LPFloat>((chunk_x * cells_per_dim_) + cell_x) / total_cells_per_dim;
    const auto max_x = static_cast<LPFloat>((chunk_x * cells_per_dim_) + cell_x + 1) / total_cells_per_dim;

    const auto min_y = static_cast<LPFloat>((chunk_y * cells_per_dim_) + cell_y) / total_cells_per_dim;
    const auto max_y = static_cast<LPFloat>((chunk_y * cells_per_dim_) + cell_y + 1) / total_cells_per_dim;

    const LPFloat radius_sq = radius * radius;

    const LPFloat closest_x = std::clamp(center_x, min_x, max_x);
    const LPFloat closest_y = std::clamp(center_y, min_y, max_y);

    const LPFloat dx_min = center_x - closest_x;
    const LPFloat dy_min = center_y - closest_y;

    const LPFloat min_dist_sq = dx_min * dx_min + dy_min * dy_min;

    if (min_dist_sq > radius_sq) {
        return CellBallRelation::OUTSIDE;
    }

    const LPFloat dx1 = center_x - min_x;
    const LPFloat dx2 = center_x - max_x;
    const LPFloat dy1 = center_y - min_y;
    const LPFloat dy2 = center_y - max_y;

    const LPFloat max_dx_sq = std::max(dx1 * dx1, dx2 * dx2);
    const LPFloat max_dy_sq = std::max(dy1 * dy1, dy2 * dy2);

    const LPFloat max_dist_sq = max_dx_sq + max_dy_sq;

    if (max_dist_sq <= radius_sq) {
        return CellBallRelation::INSIDE;
    }

    return CellBallRelation::PARTIAL;
}

void HyperRGG2D::AddWholeCellRange(
    const SInt neighbor_chunk_id, const SInt neighbor_cell_id, std::vector<PinRange>& ranges) {
    const SInt  global_cell_id = ComputeGlobalCellId(neighbor_chunk_id, neighbor_cell_id);
    const auto& cell           = cells_[global_cell_id];

    const SInt size  = std::get<0>(cell);
    const SInt begin = std::get<4>(cell);
    const SInt end   = begin + size;

    if (size > 0) {
        ranges.push_back({.begin = begin, .end = end});
    }
}

void HyperRGG2D::GenerateEdges(SInt chunk_row, SInt chunk_column) {
    const SInt chunk_id = Encode(chunk_column, chunk_row);

    if (!IsLocalChunk(chunk_id)) {
        return;
    }
    if (cells_per_dim_ == 0 || chunks_per_dim_ == 0) {
        throw ConfigurationError("Invalid grid: cells_per_dim_ or chunks_per_dim_ is zero");
    }
    const SInt total_cells_per_dim = SafeTotalCellsPerDim();

    if (total_cells_per_dim <= 0) {
        throw ConfigurationError("Invalid grid: total_cells_per_dim is zero or overflowed");
    }

    const SInt total_cells = total_cells_per_dim * total_cells_per_dim;

    if (total_cells <= 0) {
        throw ConfigurationError("Invalid grid: total_cells is zero or overflowed");
    }

    for (SInt cell_row = 0; cell_row < cells_per_dim_; ++cell_row) {
        for (SInt cell_column = 0; cell_column < cells_per_dim_; ++cell_column) {
            const SInt cell_id = EncodeCell(cell_column, cell_row);

            const SInt global_cell_x = chunk_column * cells_per_dim_ + cell_column;
            const SInt global_cell_y = chunk_row * cells_per_dim_ + cell_row;

            const SInt global_cell_id = global_cell_y * total_cells_per_dim + global_cell_x;

            const SInt base_m    = config_.m / total_cells;
            const SInt remainder = config_.m % total_cells;

            const SInt cell_m = base_m + static_cast<SInt>(global_cell_id < remainder);

            const SInt first_sample_id = global_cell_id * base_m + std::min<SInt>(global_cell_id, remainder);

            const LPFloat cell_width = LPFloat{1.0} / static_cast<LPFloat>(total_cells_per_dim);

            SInt emitted  = 0;
            SInt attempts = 0;

            const SInt max_attempts = std::max<SInt>(100, 100 * cell_m);

            while (emitted < cell_m && attempts < max_attempts) {
                const SInt sampled_center_id = (global_cell_id * max_attempts) + attempts;

                const SInt seed_x = sampling::Spooky::hash(config_.seed + 17 * sampled_center_id);
                const SInt seed_y = sampling::Spooky::hash(config_.seed + 31 * sampled_center_id);

                const LPFloat u = rng_.GenerateUniform<LPFloat>(seed_x);
                const LPFloat v = rng_.GenerateUniform<LPFloat>(seed_y);

                const LPFloat center_x = (static_cast<LPFloat>(global_cell_x) + u) * cell_width;
                const LPFloat center_y = (static_cast<LPFloat>(global_cell_y) + v) * cell_width;

                if (GenerateBallHyperedge(center_x, center_y, sampled_center_id, chunk_id, cell_id)) {
                    ++emitted;
                }

                ++attempts;
            }
        }
    }
}

void HyperRGG2D::PushHyperedgeCompressed(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges) {
    if (graph_.hyperedge_offsets.empty()) {
        graph_.hyperedge_offsets.push_back(0);
    }

    if (graph_.hyperedge_range_offsets.empty()) {
        graph_.hyperedge_range_offsets.push_back(0);
    }

    graph_.hyperedge_pins.insert(graph_.hyperedge_pins.end(), pins.begin(), pins.end());

    graph_.hyperedge_offsets.push_back(graph_.hyperedge_pins.size());

    graph_.hyperedge_ranges.insert(graph_.hyperedge_ranges.end(), ranges.begin(), ranges.end());

    graph_.hyperedge_range_offsets.push_back(graph_.hyperedge_ranges.size());
}

void HyperRGG2D::GenerateCSR() {
    GenerateGeometry();
}

void HyperRGG2D::FinalizeCSR(MPI_Comm) {}

} // namespace kagen