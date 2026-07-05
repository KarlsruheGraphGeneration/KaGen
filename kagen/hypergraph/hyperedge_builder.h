#pragma once

#include "kagen/context.h"
#include "kagen/hypergraph/debug_logger_geometric.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"

#include <algorithm>
#include <chrono>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace kagen {

template <typename GeometryPolicy>
class HyperedgeBuilder {
public:
    using Center = typename GeometryPolicy::Center;
    using Cell   = typename GeometryPolicy::Cell;
    using Vertex = typename GeometryPolicy::Vertex;

    using Clock  = std::chrono::steady_clock;
    using Double = decltype(std::declval<GeometryPolicy&>().Radius(std::declval<const Center&>()));

    explicit HyperedgeBuilder(GeometryPolicy& geometry, const PGeneratorConfig& config)
        : geometry_(geometry),
          partial_cell_mode_(config.partial_cell_mode),
          config_(config) {
        if (config_.debug) {
            logger_.emplace(MakeDebugFilename(), true);
        }
    }

    std::string MakeDebugFilename() const {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::string output = config_.output_graph.filename + "_" + std::to_string(config_.n) + "_"
                             + std::to_string(config_.m) + "_" + std::to_string(config_.hyperedge_radius_exponent)
                             + "_debug_rank_" + std::to_string(rank) + ".csv";
        return output;
    }

    void Build(const Center& center) {
        ResetBuildState();
        const auto timer_start_total = config_.debug ? Clock::now() : Clock::time_point{};

        const auto radius = CollectCandidateCells(center);

        BuildStats stats;

        for (const Cell& cell: cells_) {
            ProcessCell(center, radius, cell, stats);
        }

        FinalizePinsAndRanges();

        if (config_.debug && logger_) {
            LogBuildStats(center, radius, stats, timer_start_total);
        }

        geometry_.EmitHyperedge(pins_, ranges_);
    }

private:
    struct BuildStats {
        SInt outside_cells          = 0;
        SInt inside_cells           = 0;
        SInt partial_cells          = 0;
        SInt partial_estimated_size = 0;
        SInt inside_estimated_size  = 0;
    };

    void ResetBuildState() {
        cells_.clear();
        pins_.clear();
        ranges_.clear();
    }

    void FinalizePinsAndRanges() {
        if (partial_cell_mode_ == PartialCellMode::GenerateAndCheck) {
            if (!ranges_.empty()) {
                NormalizeRangesOnly(ranges_);
            }
        } else {
            Normalize(pins_, ranges_);
        }
    }

    std::vector<PinRange> MergeRanges(std::vector<PinRange>& ranges) const {
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
        return merged_ranges;
    }

    std::vector<SInt>
    RemovePinsCoveredByRanges(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges) const {
        std::vector<SInt> filtered_pins;
        filtered_pins.reserve(pins.size());

        std::size_t range_index = 0;

        for (const SInt pin: pins) {
            while (range_index < ranges.size() && ranges[range_index].end <= pin) {
                ++range_index;
            }

            const bool covered =
                range_index < ranges.size() && ranges[range_index].begin <= pin && pin < ranges[range_index].end;

            if (!covered) {
                filtered_pins.push_back(pin);
            }
        }

        return filtered_pins;
    }

    void Normalize(std::vector<SInt>& pins, std::vector<PinRange>& ranges) const {
        std::sort(pins.begin(), pins.end());
        pins.erase(std::unique(pins.begin(), pins.end()), pins.end());

        std::sort(ranges.begin(), ranges.end(), [](const PinRange& firstRange, const PinRange& secondRange) {
            return firstRange.begin < secondRange.begin;
        });

        std::vector<PinRange> merged_ranges;
        merged_ranges.reserve(ranges.size());

        merged_ranges = MergeRanges(ranges);

        auto filtered_pins = RemovePinsCoveredByRanges(pins, merged_ranges);

        pins   = std::move(filtered_pins);
        ranges = std::move(merged_ranges);
    }

    void NormalizeRangesOnly(std::vector<PinRange>& ranges) const {
        std::sort(ranges.begin(), ranges.end(), [](const PinRange& a, const PinRange& b) { return a.begin < b.begin; });

        auto merged = MergeRanges(ranges);

        ranges = std::move(merged);
    }

    void CountOutside(BuildStats& stats) const {
        if (config_.debug) {
            ++stats.outside_cells;
        }
    }

    void CountInside(BuildStats& stats, const SInt added) const {
        if (config_.debug) {
            stats.inside_estimated_size += added;
            ++stats.inside_cells;
        }
    }

    void CountPartial(BuildStats& stats, const SInt added) const {
        if (config_.debug) {
            stats.partial_estimated_size += added;
            ++stats.partial_cells;
        }
    }

    void ProcessPartialCell(const Center& center, Double radius, const Cell& cell, BuildStats& stats) {
        if (partial_cell_mode_ == PartialCellMode::GenerateAndCheck) {
            CountPartial(stats, geometry_.AddPartialCellExact(center, radius, cell, pins_));
            return;
        }
        const Double coverage = geometry_.CellCoverage(center, radius, cell);

        if (coverage <= 0.0) {
            CountOutside(stats);
            return;
        }

        if (coverage >= 1.0) {
            CountInside(stats, geometry_.AddWholeCell(cell, ranges_));
            return;
        }
        if (partial_cell_mode_ == PartialCellMode::EstimateByCoverageRange) {
            CountPartial(stats, geometry_.AddPartialCellRange(center, cell, coverage, pins_, ranges_));
            return;
        }
        CountPartial(stats, geometry_.AddPartialCellFloyd(center, cell, coverage, pins_, ranges_));
    }

    void ProcessCell(const Center& center, Double radius, const Cell& cell, BuildStats& stats) {
        const CellBallRelation relation = geometry_.ClassifyCell(center, radius, cell);

        switch (relation) {
            case CellBallRelation::OUTSIDE: {
                CountOutside(stats);
                return;
            }

            case CellBallRelation::INSIDE: {
                CountInside(stats, geometry_.AddWholeCell(cell, ranges_));
                return;
            }

            case CellBallRelation::PARTIAL: {
                ProcessPartialCell(center, radius, cell, stats);
                return;
            }
        }
    }

    void LogBuildStats(const Center& center, Double radius, const BuildStats& stats, Clock::time_point start) {
        const auto      timer_end_total = Clock::now();
        const long long duration_ns =
            std::chrono::duration_cast<std::chrono::nanoseconds>(timer_end_total - start).count();
        SInt hyperedge_size = static_cast<SInt>(pins_.size());
        for (const PinRange& range: ranges_) {
            hyperedge_size += range.end - range.begin;
        }

        ++counter_;
        logger_->LogHyperedge(
            counter_, geometry_.CenterToString(center), radius, static_cast<SInt>(cells_.size()), stats.inside_cells,
            stats.partial_cells, stats.outside_cells, static_cast<SInt>(pins_.size()),
            static_cast<SInt>(ranges_.size()), hyperedge_size, duration_ns, stats.inside_estimated_size,
            stats.partial_estimated_size);
    }

    Double CollectCandidateCells(const Center& center) {
        geometry_.AddCenter(center, pins_);
        const auto radius = geometry_.Radius(center);
        geometry_.CandidateCells(center, radius, cells_);
        return radius;
    }

    GeometryPolicy&                               geometry_;
    PartialCellMode                               partial_cell_mode_;
    std::vector<Cell>                             cells_;
    const PGeneratorConfig&                       config_;
    std::optional<GeometricHypergraphDebugLogger> logger_;
    std::vector<SInt>                             pins_;
    std::vector<PinRange>                         ranges_;
    SInt                                          counter_ = 0;
};

} // namespace kagen
