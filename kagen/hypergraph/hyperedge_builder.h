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
#include "kagen/hypergraph/debug_logger.h"

#include <CGAL/number_utils_classes.h>

namespace kagen {

template <typename GeometryPolicy>
class HyperedgeBuilder {
public:
    using Center = typename GeometryPolicy::Center;
    using Cell   = typename GeometryPolicy::Cell;

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
        const auto timer_start_total = Clock::now();

        std::vector<SInt>     pins;
        std::vector<PinRange> ranges;

        geometry_.AddCenter(center, pins);
        const auto radius = geometry_.Radius(center);
        const auto cells  = geometry_.CandidateCells(center, radius);

        SInt outside_cells          = 0;
        SInt inside_cells           = 0;
        SInt partial_cells          = 0;
        SInt partial_estimated_size = 0;
        SInt inside_estimated_size  = 0;

        for (const Cell& cell: cells) {
            const CellBallRelation relation = geometry_.ClassifyCell(center, radius, cell);

            switch (relation) {
                case CellBallRelation::OUTSIDE:
                    ++outside_cells;
                    continue;

                case CellBallRelation::INSIDE:
                    inside_estimated_size += geometry_.AddWholeCell(cell, ranges);
                    ++inside_cells;
                    continue;

                case CellBallRelation::PARTIAL:
                    if (partial_cell_mode_ == PartialCellMode::GenerateAndCheck) {
                        partial_estimated_size += geometry_.AddPartialCellExact(center, radius, cell, pins);
                        ++partial_cells;
                    } else {
                        const double coverage = geometry_.CellCoverage(center, radius, cell);

                        if (coverage <= 0.0) {
                            ++outside_cells;
                            continue;
                        }

                        if (coverage >= 1.0) {
                            inside_estimated_size += geometry_.AddWholeCell(cell, ranges);
                            ++inside_cells;
                        } else {
                            partial_estimated_size += geometry_.AddPartialCell(center, cell, coverage, ranges);
                            ++partial_cells;
                        }
                    }
                    continue;
            }
        }

        Normalize(pins, ranges);
        SInt hyperedge_size = static_cast<SInt>(pins.size());
        for (const PinRange& range: ranges) {
            hyperedge_size += range.end - range.begin;
        }

        const auto timer_end_total = Clock::now();

        const long long duration_ns =
            std::chrono::duration_cast<std::chrono::nanoseconds>(timer_end_total - timer_start_total).count();

        ++counter_;
        std::string center_string = geometry_.CenterToString(center);
        if (logger_) {
            logger_->LogHyperedge(
                counter_, center_string, static_cast<double>(radius), static_cast<SInt>(cells.size()), inside_cells,
                partial_cells, outside_cells, static_cast<SInt>(pins.size()), static_cast<SInt>(ranges.size()),
                hyperedge_size, duration_ns, inside_estimated_size, partial_estimated_size);
        }

        geometry_.EmitHyperedge(pins, ranges);
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

    GeometryPolicy                       geometry_;
    PartialCellMode                      partial_cell_mode_;
    PGeneratorConfig                     config_;
    std::optional<HypergraphDebugLogger> logger_;
    SInt                                 counter_ = 0;
    int                                  rank_    = 0;
};

} // namespace kagen
