#pragma once

#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"
#include "kagen/generators/hyper/h_geometric/hyper_rgg2d.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"

#include <string>
#include <vector>

namespace kagen {

class HyperRGG2DPolicy {
public:
    using Vertex = kagen::Vertex;
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

        SInt global_cell_id;
    };

    explicit HyperRGG2DPolicy(HyperRGG2D& generator);

    void AddCenter(const Center& center, std::vector<SInt>& pins) const;

    LPFloat Radius(const Center& center) const;

    SInt GetNumVerticeOfCellCoord(SSInt global_cell_x, SSInt global_cell_y);

    void CandidateCells(const Center& center, LPFloat radius, std::vector<Cell>& cells) const;

    CellBallRelation ClassifyCell(const Center& center, LPFloat radius, const Cell& cell) const;

    double CellCoverage(const Center& center, LPFloat radius, const Cell& cell) const;

    SInt AddWholeCell(const Cell& cell, std::vector<PinRange>& ranges) const;

    SInt AddPartialCellRange(
        const Center& center, const Cell& cell, double coverage, std::vector<SInt>& pins,
        std::vector<PinRange>& ranges) const;

    SInt AddPartialCellFloyd(
        const Center& center, const Cell& cell, double coverage, std::vector<SInt>& pins,
        std::vector<PinRange>& ranges) const;

    SInt AddPartialCellExact(const Center& center, LPFloat radius, const Cell& cell, std::vector<SInt>& pins) const;

    void EmitHyperedge(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges);

    static std::string CenterToString(const Center center) {
        return std::to_string(center.x) + ";" + std::to_string(center.y);
    }

    void PrintExactStatistics() const;

private:
    double MinimumRadius(const Center& center);

    double EstimatedCellCoverage(const Center& center, LPFloat radius, const Cell& cell) const;

    double EstimatedCircleRectCoverage(
        double center_x, double center_y, double min_x, double max_x, double min_y, double max_y, double radius) const;

    HyperRGG2D* gen_;

    // Cache
    struct CachedExactCell {
        std::vector<Vertex> vertices_by_x;
    };

    mutable std::unordered_map<SInt, CachedExactCell> exact_vertices_;
    mutable std::unordered_set<SInt>                  floyd_scratch_;

    mutable RNGWrapper<> rng_;

    const CachedExactCell& ExactCell(const Cell& cell) const {
        auto it = exact_vertices_.find(cell.global_cell_id);
        if (it != exact_vertices_.end()) {
            return it->second;
        }

        auto [inserted_it, inserted] = exact_vertices_.emplace(cell.global_cell_id, CachedExactCell{});

        auto& cached = inserted_it->second;
        gen_->GenerateVertices(cell.chunk_id, cell.cell_id, cached.vertices_by_x);

        std::sort(cached.vertices_by_x.begin(), cached.vertices_by_x.end(), [](const Vertex& a, const Vertex& b) {
            return std::get<0>(a) < std::get<0>(b);
        });

        return cached;
    }
};

} // namespace kagen
