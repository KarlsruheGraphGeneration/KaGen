#pragma once

#include "kagen/context.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"

#include <algorithm>
#include <chrono>
#include <optional>
#include <string>
#include <utility>
#include <vector>
// TODO(clickup)[2026-05-24]: Remove after Debugging
#include "kagen/hypergraph/debug_logger_geometric.h"

#include <CGAL/number_utils_classes.h>

namespace kagen {

template <typename GeometryPolicy>
class HyperedgeBuilder {
public:
    using Center = typename GeometryPolicy::Center;
    using Cell   = typename GeometryPolicy::Cell;
    using Vertex = typename GeometryPolicy::Vertex;

    using Clock = std::chrono::steady_clock;

    explicit HyperedgeBuilder(GeometryPolicy& geometry, const PGeneratorConfig& config)
        : geometry_(geometry),
          partial_cell_mode_(config.partial_cell_mode),
          config_(config) {
        if (config_.debug) {
            logger_.emplace(MakeDebugFilename(), true);
        }
    }

    std::string MakeDebugFilename() {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::string output = config_.output_graph.filename + "_" + std::to_string(config_.n) + "_"
                             + std::to_string(config_.m) + "_" + std::to_string(config_.hyperedge_radius_exponent)
                             + "_debug_rank_" + std::to_string(rank) + ".csv";
        return output;
    }

    void Build(const Center& center) {
        cells_.clear();
        pins_.clear();
        ranges_.clear();
        const auto timer_start_total = config_.debug ? Clock::now() : Clock::time_point{};

        geometry_.AddCenter(center, pins_);
        const auto radius = geometry_.Radius(center);
        geometry_.CandidateCells(center, radius, cells_);

        SInt outside_cells          = 0;
        SInt inside_cells           = 0;
        SInt partial_cells          = 0;
        SInt partial_estimated_size = 0;
        SInt inside_estimated_size  = 0;

        for (const Cell& cell: cells_) {
            const CellBallRelation relation = geometry_.ClassifyCell(center, radius, cell);

            switch (relation) {
                case CellBallRelation::OUTSIDE: {
                    if (config_.debug) {
                        ++outside_cells;
                    }
                    continue;
                }

                case CellBallRelation::INSIDE: {
                    const SInt added = geometry_.AddWholeCell(cell, ranges_);
                    if (config_.debug) {
                        inside_estimated_size += added;
                        ++inside_cells;
                    }
                    continue;
                }

                case CellBallRelation::PARTIAL: {
                    if (partial_cell_mode_ == PartialCellMode::GenerateAndCheck) {
                        const SInt added = geometry_.AddPartialCellExact(center, radius, cell, pins_);
                        if (config_.debug) {
                            partial_estimated_size += added;
                            ++partial_cells;
                        }
                    } else {
                        const double coverage = geometry_.CellCoverage(center, radius, cell);

                        if (coverage <= 0.0) {
                            if (config_.debug) {
                                ++outside_cells;
                            }
                            continue;
                        }

                        if (coverage >= 1.0) {
                            const SInt added = geometry_.AddWholeCell(cell, ranges_);
                            if (config_.debug) {
                                inside_estimated_size += added;
                                ++inside_cells;
                            }
                        } else if (partial_cell_mode_ == PartialCellMode::EstimateByCoverageRange) {
                            const SInt added = geometry_.AddPartialCellRange(center, cell, coverage, pins_, ranges_);
                            if (config_.debug) {
                                partial_estimated_size += added;
                                ++partial_cells;
                            }
                        } else {
                            const SInt added = geometry_.AddPartialCellFloyd(center, cell, coverage, pins_, ranges_);
                            if (config_.debug) {
                                partial_estimated_size += added;
                                ++partial_cells;
                            }
                        }
                    }
                    continue;
                }
            }
        }

        if (partial_cell_mode_ == PartialCellMode::GenerateAndCheck) {
            if (!ranges_.empty()) {
                NormalizeRangesOnly(ranges_);
            }
        } else {
            Normalize(pins_, ranges_);
        }

        if (config_.debug && logger_) {
            const auto      timer_end_total = Clock::now();
            const long long duration_ns =
                std::chrono::duration_cast<std::chrono::nanoseconds>(timer_end_total - timer_start_total).count();
            SInt hyperedge_size = static_cast<SInt>(pins_.size());
            for (const PinRange& range: ranges_) {
                hyperedge_size += range.end - range.begin;
            }

            ++counter_;
            logger_->LogHyperedge(
                counter_, geometry_.CenterToString(center), static_cast<double>(radius),
                static_cast<SInt>(cells_.size()), inside_cells, partial_cells, outside_cells,
                static_cast<SInt>(pins_.size()), static_cast<SInt>(ranges_.size()), hyperedge_size, duration_ns,
                inside_estimated_size, partial_estimated_size);
        }

        geometry_.EmitHyperedge(pins_, ranges_);
    }

private:
    void Normalize(std::vector<SInt>& pins, std::vector<PinRange>& ranges) const {
        std::sort(pins.begin(), pins.end());
        pins.erase(std::unique(pins.begin(), pins.end()), pins.end());

        std::sort(ranges.begin(), ranges.end(), [](const PinRange& firstRange, const PinRange& secondRange) {
            return firstRange.begin < secondRange.begin;
        });

        std::vector<PinRange> merged_ranges;

        for (const PinRange& range: ranges) {
            if (range.begin >= range.end) {
                continue;
            }

            if (!merged_ranges.empty() && merged_ranges.back().end >= range.begin) {
                merged_ranges.back().end = std::max(merged_ranges.back().end, range.end);
            } else {
                merged_ranges.push_back(range);
            }
        }

        std::vector<SInt> filtered_pins;
        filtered_pins.reserve(pins.size());

        std::size_t range_index = 0;

        for (const SInt pin: pins) {
            while (range_index < merged_ranges.size() && merged_ranges[range_index].end <= pin) {
                ++range_index;
            }

            const bool covered = range_index < merged_ranges.size() && merged_ranges[range_index].begin <= pin
                                 && pin < merged_ranges[range_index].end;

            if (!covered) {
                filtered_pins.push_back(pin);
            }
        }

        pins   = std::move(filtered_pins);
        ranges = std::move(merged_ranges);
    }

    void NormalizeRangesOnly(std::vector<PinRange>& ranges) const {
        std::sort(ranges.begin(), ranges.end(), [](const PinRange& a, const PinRange& b) { return a.begin < b.begin; });

        std::vector<PinRange> merged;
        merged.reserve(ranges.size());

        for (const auto& range: ranges) {
            if (range.begin >= range.end) {
                continue;
            }

            if (!merged.empty() && merged.back().end >= range.begin) {
                merged.back().end = std::max(merged.back().end, range.end);
            } else {
                merged.push_back(range);
            }
        }

        ranges = std::move(merged);
    }

    GeometryPolicy&                               geometry_;
    PartialCellMode                               partial_cell_mode_;
    std::vector<Cell>                             cells_;
    const PGeneratorConfig&                       config_;
    std::optional<GeometricHypergraphDebugLogger> logger_;
    std::vector<SInt>                             pins_;
    std::vector<PinRange>                         ranges_;
    SInt                                          counter_ = 0;
    int                                           rank_    = 0;
};

} // namespace kagen
