#pragma once

#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"

#include <algorithm>
#include <chrono>
#include <utility>
#include <vector>
// TODO(clickup)[2026-05-24]: Remove after Debugging
#include "kagen/hypergraph/debug_logger.h"

#include <CGAL/number_utils_classes.h>

#include <fstream>

namespace kagen {

template <typename GeometryPolicy>
class HyperedgeBuilder {
public:
    using Center = typename GeometryPolicy::Center;
    using Cell   = typename GeometryPolicy::Cell;

    using Clock = std::chrono::steady_clock;

    HypergraphDebugLogger logger_;
    SInt                  counter_ = 0;
    int                   rank_    = 0;

    explicit HyperedgeBuilder(GeometryPolicy& geometry, const PartialCellMode partial_cell_mode)
        : geometry_(geometry),
          partial_cell_mode_(partial_cell_mode),
          logger_(MakeDebugFilename(), true) {}

    static std::string MakeDebugFilename() {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return "output/hyper_debug_rank_" + std::to_string(rank) + ".csv";
    }

    void Build(const Center& center) {
        const auto t1 = Clock::now();

        std::vector<SInt>     pins;
        std::vector<PinRange> ranges;

        geometry_.AddCenter(center, pins);

        const auto radius = geometry_.Radius(center);
        const auto cells  = geometry_.CandidateCells(center, radius);

        SInt outside_cells = 0;
        SInt inside_cells  = 0;
        SInt partial_cells = 0;

        for (const Cell& cell: cells) {
            const CellBallRelation relation = geometry_.ClassifyCell(center, radius, cell);

            switch (relation) {
                case CellBallRelation::OUTSIDE:
                    ++outside_cells;
                    continue;

                case CellBallRelation::INSIDE:
                    geometry_.AddWholeCell(cell, ranges);
                    ++inside_cells;
                    continue;

                case CellBallRelation::PARTIAL:
                    if (partial_cell_mode_ == PartialCellMode::GenerateAndCheck) {
                        geometry_.AddPartialCellExact(center, radius, cell, pins);
                        ++partial_cells;
                    } else {
                        const double coverage = geometry_.CellCoverage(center, radius, cell);

                        if (coverage <= 0.0) {
                            ++outside_cells;
                            continue;
                        }

                        if (coverage >= 1.0) {
                            geometry_.AddWholeCell(cell, ranges);
                            ++inside_cells;
                        } else {
                            geometry_.AddPartialCell(center, cell, coverage, ranges);
                            ++partial_cells;
                        }
                    }
                    continue;
            }
        }

        Normalize(pins, ranges);
        static SInt counter = 0;
        if (counter % 1000 == 0) {
            std::cout << counter << '\n';
        }
        ++counter;
        SInt hyperedge_size = static_cast<SInt>(pins.size());
        for (const PinRange& range: ranges) {
            hyperedge_size += range.end - range.begin;
        }

        const auto t2 = Clock::now();

        const long long duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();

        ++counter_;
        std::string center_string = "TODO";

        logger_.LogHyperedge(
            counter_, center_string, static_cast<double>(radius), static_cast<SInt>(cells.size()), inside_cells, partial_cells,
            outside_cells, static_cast<SInt>(pins.size()), static_cast<SInt>(ranges.size()), hyperedge_size,
            duration_ns);

        geometry_.EmitHyperedge(pins, ranges);
    }

private:
    void Normalize(std::vector<SInt>& pins, std::vector<PinRange>& ranges) const {
        std::sort(pins.begin(), pins.end());
        pins.erase(std::unique(pins.begin(), pins.end()), pins.end());

        std::sort(ranges.begin(), ranges.end(), [](const PinRange& a, const PinRange& b) { return a.begin < b.begin; });

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

        for (const SInt pin: pins) {
            bool covered = false;

            for (const PinRange& range: merged_ranges) {
                if (pin < range.begin) {
                    break;
                }

                if (pin >= range.begin && pin < range.end) {
                    covered = true;
                    break;
                }
            }

            if (!covered) {
                filtered_pins.push_back(pin);
            }
        }

        pins   = std::move(filtered_pins);
        ranges = std::move(merged_ranges);
    }

    GeometryPolicy& geometry_;
    PartialCellMode partial_cell_mode_;
};

} // namespace kagen
