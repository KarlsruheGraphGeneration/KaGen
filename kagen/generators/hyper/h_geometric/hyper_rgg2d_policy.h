#pragma once

#include "kagen/generators/hyper/h_geometric/hyper_rgg2d.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"

#include <string>
#include <vector>

namespace kagen {

class HyperRGG2DPolicy {
public:
    struct Center {
        LPFloat x;
        LPFloat y;
        SInt    sampled_id;
        SInt    chunk_id;
        SInt    cell_id;
    };

    struct Cell {
        SInt chunk_id;
        SInt cell_id;
        SInt global_cell_x;
        SInt global_cell_y;
    };

    explicit HyperRGG2DPolicy(HyperRGG2D& generator);

    void AddCenter(const Center& center, std::vector<SInt>& pins) const;

    LPFloat Radius(const Center& center) const;

    SInt GetNumVerticeOfCellCoord(SSInt global_cell_x, SSInt global_cell_y);

    std::vector<Cell> CandidateCells(const Center& center, LPFloat radius) const;

    CellBallRelation ClassifyCell(const Center& center, LPFloat radius, const Cell& cell) const;

    double CellCoverage(const Center& center, LPFloat radius, const Cell& cell) const;

    SInt AddWholeCell(const Cell& cell, std::vector<PinRange>& ranges) const;

    SInt AddPartialCell(const Center& center, const Cell& cell, double coverage, std::vector<PinRange>& ranges) const;

    SInt AddPartialCellExact(const Center& center, LPFloat radius, const Cell& cell, std::vector<SInt>& pins) const;

    void EmitHyperedge(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges);

    static std::string CenterToString(const Center center) {
        return std::to_string(center.x) + ";" + std::to_string(center.y);
    }

private:
    double MinimumRadius(const Center& center);

    double EstimatedCellCoverage(const Center& center, LPFloat radius, const Cell& cell) const;

    double EstimatedCircleRectCoverage(
        double center_x, double center_y, double min_x, double max_x, double min_y, double max_y, double radius) const;

    HyperRGG2D* gen_;

    mutable RNGWrapper<> rng_;
};

} // namespace kagen
