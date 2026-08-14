#pragma once

#include "kagen/generators/hyper/h_geometric/hyper_rgg2d.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"

#include <list>
#include <string>
#include <vector>

namespace kagen {

class HyperRGG2DPolicy {
public:
    using Vertex = kagen::Vertex;
    using Center = HyperRGG2D::Center;

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

    double MinimumRadius(const Center& center);

private:
    double EstimatedCircleRectCoverage(
        double center_x, double center_y, double min_x, double max_x, double min_y, double max_y, double radius) const;
    double EstimateCoverageRecursive(
        double center_x, double center_y, double radius_sq, double min_x, double max_x, double min_y, double max_y,
        int depth) const;

    static constexpr int coverage_max_depth_ = 8;

    std::optional<Cell> TryMakeCell(SSInt global_cell_x, SSInt global_cell_y) const;

    struct StoredCell {
        SInt size;
        SInt offset;
    };

    std::optional<StoredCell> TryGetStoredCell(const Cell& cell) const;

    struct PartialCellSample {
        StoredCell stored;
        SInt       count;
    };

    std::optional<PartialCellSample> PreparePartialCellSample(const Cell& cell, double coverage) const;

    struct CellBounds {
        double min_x;
        double max_x;
        double min_y;
        double max_y;
    };

    CellBounds GetCellBounds(const Cell& cell) const;

    HyperRGG2D* gen_;

    // Cache
    struct CachedExactCell {
        std::vector<Vertex> vertices_by_x;
    };

    mutable std::unordered_map<SInt, CachedExactCell>           exact_vertices_;
    mutable std::list<SInt>                                     exact_lru_;
    mutable std::unordered_map<SInt, std::list<SInt>::iterator> exact_lru_pos_;
    std::size_t                                                 exact_remote_cache_limit_ = 8; // configurable later
    mutable std::unordered_set<SInt>                            floyd_scratch_;

    mutable SInt exact_remote_cache_hits_      = 0;
    mutable SInt exact_remote_cache_misses_    = 0;
    mutable SInt exact_remote_cached_vertices_ = 0;
    mutable SInt exact_remote_live_vertices_   = 0;
    mutable SInt exact_remote_peak_vertices_   = 0;

    mutable SInt exact_remote_access_counter_ = 0;

    mutable std::unordered_map<SInt, SInt> exact_remote_last_access_;

    mutable SInt exact_remote_reuse_count_        = 0;
    mutable SInt exact_remote_reuse_distance_sum_ = 0;
    mutable SInt exact_remote_reuse_distance_max_ = 0;

    mutable SInt exact_remote_reuse_distance_le_1_  = 0;
    mutable SInt exact_remote_reuse_distance_le_4_  = 0;
    mutable SInt exact_remote_reuse_distance_le_16_ = 0;
    mutable SInt exact_remote_reuse_distance_gt_16_ = 0;

    mutable SInt local_exact_access_counter_     = 0;
    mutable SInt local_exact_reuse_count_        = 0;
    mutable SInt local_exact_reuse_distance_sum_ = 0;
    mutable SInt local_exact_reuse_distance_max_ = 0;

    mutable SInt local_exact_reuse_distance_le_1_  = 0;
    mutable SInt local_exact_reuse_distance_le_4_  = 0;
    mutable SInt local_exact_reuse_distance_le_16_ = 0;
    mutable SInt local_exact_reuse_distance_gt_16_ = 0;

    mutable HashMap<SInt, SInt> local_exact_last_access_;

    mutable RNGWrapper<> rng_;
    mutable Mersenne     rounding_rng_;

    const std::vector<Vertex>& ExactVerticesByX(const Cell& cell) const;
    const CachedExactCell&     ExactRemoteCell(const Cell& cell) const;

    mutable std::unordered_set<SInt> local_vertices_sorted_by_x_;

    bool IsLocalCell(const Cell& cell) const;
    void RecordRemoteAccess(SInt global_cell_id) const;

    void RecordLocalExactAccess(SInt global_cell_id) const;

public:
    void PrintExactCacheStats() const;

    bool ShouldApproximatePartialCell(const Cell& cell) const;
};

} // namespace kagen
