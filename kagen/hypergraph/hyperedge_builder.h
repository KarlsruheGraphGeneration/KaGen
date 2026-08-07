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

    explicit HyperedgeBuilder(
        GeometryPolicy& geometry, const PGeneratorConfig& config, GeometricHypergraphDebugLogger* logger)
        : geometry_(geometry),
          partial_cell_mode_(config.partial_cell_mode),
          config_(config),
          logger_(logger) {}

    void Build(const Center& center) {
        ResetBuildState();
        const auto timer_start_total = config_.debug ? Clock::now() : Clock::time_point{};

        const auto radius = CollectCandidateCells(center);

        BuildStats stats;

        std::vector<CellBallRelation> relations;
        relations.reserve(cells_.size());

        bool has_inside_cell = false;

        for (const Cell& cell: cells_) {
            const CellBallRelation relation = geometry_.ClassifyCell(center, radius, cell);

            relations.push_back(relation);

            if (relation == CellBallRelation::INSIDE) {
                has_inside_cell = true;
            }
        }

        const bool approximation_allowed = has_inside_cell;

        for (std::size_t i = 0; i < cells_.size(); ++i) {
            ProcessCell(center, radius, cells_[i], relations[i], approximation_allowed, stats);
        }

        FinalizePinsAndRanges();

        if (config_.debug && logger_) {
            LogBuildStats(center, radius, stats, timer_start_total);
        }

        geometry_.EmitHyperedge(pins_, ranges_);
    }

private:
    struct BuildStats {
        SInt outside_cells            = 0;
        SInt inside_cells             = 0;
        SInt partial_cells            = 0;
        SInt partial_estimated_size   = 0;
        SInt inside_estimated_size    = 0;
        SInt num_vertice_remote_chunk = 0;
    };

    template <typename T>
    void ClearMaybeRelease(std::vector<T>& v, const std::size_t max_keep_capacity = 4096) {
        if (v.capacity() > max_keep_capacity) {
            std::vector<T>().swap(v);
        } else {
            v.clear();
        }
    }

    void ResetBuildState() {
        ClearMaybeRelease(cells_, 1 << 15);
        ClearMaybeRelease(pins_, 1 << 20);
        ClearMaybeRelease(ranges_, 1 << 18);
    }

    void FinalizePinsAndRanges() {
        if (pins_.empty()) {
            NormalizeRangesOnly(ranges_);
        } else {
            Normalize(pins_, ranges_);
        }
    }
    void RemovePinsCoveredByRanges(std::vector<SInt>& pins, const std::vector<PinRange>& ranges) const {
        std::size_t range_index = 0;
        std::size_t write       = 0;

        for (const SInt pin: pins) {
            while (range_index < ranges.size() && ranges[range_index].end <= pin) {
                ++range_index;
            }

            const bool covered =
                range_index < ranges.size() && ranges[range_index].begin <= pin && pin < ranges[range_index].end;

            if (!covered) {
                pins[write++] = pin;
            }
        }

        pins.resize(write);
    }

    void ExtractRunsAsRanges(std::vector<SInt>& pins, std::vector<PinRange>& ranges) const {
        constexpr SInt min_run_length = 3;

        std::size_t read  = 0;
        std::size_t write = 0;

        while (read < pins.size()) {
            std::size_t end = read + 1;

            while (end < pins.size() && pins[end] == pins[end - 1] + 1) {
                ++end;
            }

            const std::size_t run_length = end - read;

            if (run_length >= static_cast<std::size_t>(min_run_length)) {
                ranges.push_back({
                    .begin = pins[read],
                    .end   = pins[end - 1] + 1,
                });
            } else {
                for (std::size_t i = read; i < end; ++i) {
                    pins[write++] = pins[i];
                }
            }

            read = end;
        }

        pins.resize(write);
    }

    void Normalize(std::vector<SInt>& pins, std::vector<PinRange>& ranges) const {
        std::sort(pins.begin(), pins.end());
        pins.erase(std::unique(pins.begin(), pins.end()), pins.end());

        ExtractRunsAsRanges(pins, ranges);

        NormalizeRangesOnly(ranges);

        RemovePinsCoveredByRanges(pins, ranges);
    }

    void NormalizeRangesOnly(std::vector<PinRange>& ranges) const {
        if (ranges.empty()) {
            return;
        }

        std::sort(ranges.begin(), ranges.end(), [](const PinRange& a, const PinRange& b) { return a.begin < b.begin; });

        std::size_t write = 0;

        for (const PinRange& range: ranges) {
            if (range.begin >= range.end) {
                continue;
            }

            if (write > 0 && ranges[write - 1].end >= range.begin) {
                ranges[write - 1].end = std::max(ranges[write - 1].end, range.end);
            } else {
                ranges[write++] = range;
            }
        }

        ranges.resize(write);
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

    void ProcessPartialCell(
        const Center& center, const Double radius, const Cell& cell, const bool approximation_allowed,
        BuildStats& stats) {
        const bool use_exact = partial_cell_mode_ == PartialCellMode::GenerateAndCheck || !approximation_allowed;

        if (use_exact) {
            CountPartial(stats, geometry_.AddPartialCellExact(center, radius, cell, pins_));
            return;
        }

        const Double coverage = std::clamp(geometry_.CellCoverage(center, radius, cell), Double{0.0}, Double{1.0});

        if (partial_cell_mode_ == PartialCellMode::EstimateByCoverageRange) {
            CountPartial(stats, geometry_.AddPartialCellRange(center, cell, coverage, pins_, ranges_));
            return;
        }

        CountPartial(stats, geometry_.AddPartialCellFloyd(center, cell, coverage, pins_, ranges_));
    }

    void ProcessCell(
        const Center& center, const Double radius, const Cell& cell, const CellBallRelation relation,
        const bool approximation_allowed, BuildStats& stats) {
        switch (relation) {
            case CellBallRelation::OUTSIDE:
                CountOutside(stats);
                return;

            case CellBallRelation::INSIDE: {
                const SInt added = geometry_.AddWholeCell(cell, ranges_);

                CountInside(stats, added);
                return;
            }

            case CellBallRelation::PARTIAL:
                ProcessPartialCell(center, radius, cell, approximation_allowed, stats);
                return;
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

        logger_->LogHyperedge(
            geometry_.CenterToString(center), radius, static_cast<SInt>(cells_.size()), stats.inside_cells,
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

    GeometryPolicy&                 geometry_;
    PartialCellMode                 partial_cell_mode_;
    std::vector<Cell>               cells_;
    const PGeneratorConfig&         config_;
    GeometricHypergraphDebugLogger* logger_ = nullptr;
    std::vector<SInt>               pins_;
    std::vector<PinRange>           ranges_;
};

} // namespace kagen
