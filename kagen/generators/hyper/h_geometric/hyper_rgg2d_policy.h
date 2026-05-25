#pragma once

#include "kagen/context.h"
#include "kagen/generators/hyper/h_geometric/hyper_rgg2d.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"

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

    std::vector<Cell> CandidateCells(const Center& center, LPFloat radius) const;

    CellBallRelation ClassifyCell(const Center& center, LPFloat radius, const Cell& cell) const;

    double CellCoverage(const Center& center, LPFloat radius, const Cell& cell) const;

    void AddWholeCell(const Cell& cell, std::vector<PinRange>& ranges) const;

    void AddPartialCell(const Center& center, const Cell& cell, double coverage, std::vector<PinRange>& ranges) const;

    void AddPartialCellExact(const Center& center, LPFloat radius, const Cell& cell, std::vector<SInt>& pins) const;

    void EmitHyperedge(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges);

private:
    HyperRGG2D& gen_;
};

} // namespace kagen
