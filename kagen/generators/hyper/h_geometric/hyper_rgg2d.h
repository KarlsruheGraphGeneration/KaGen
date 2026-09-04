#pragma once

#include "kagen/context.h"
#include "kagen/generators/geometric/spatial_grid_2d.h"
#include "kagen/generators/hyper/h_geometric/h_rgg.h"
#include "kagen/hypergraph/debug_logger_geometric.h"

#include <vector>

namespace kagen {

class HyperRGG2DPolicy;

class HyperRGG2D : public HRGG, protected SpatialGrid2D, private CSROnlyGenerator {
    friend class HyperRGG2DPolicy;

public:
    HyperRGG2D(const PGeneratorConfig& config, PEID rank, PEID size);

    struct Center {
        LPFloat x;
        LPFloat y;
        LPFloat radius;
        SInt    sampled_id;
        SInt    chunk_id;
        SInt    cell_id;
    };

    ~HyperRGG2D() override;

protected:
    void GenerateCSR() override;
    void FinalizeCSR(MPI_Comm comm) override;

    void AddExactDebugStats(SInt partial_cells, SInt vertices_tested, SInt vertices_accepted) const {
        if (!config_.debug) {
            return;
        }

        exact_debug_stats_.partial_cells += partial_cells;
        exact_debug_stats_.vertices_tested += vertices_tested;
        exact_debug_stats_.vertices_accepted += vertices_accepted;
    }

private:
    std::unique_ptr<HyperRGG2DPolicy> policy_;
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    struct WeakScalingInstrumentation {
        std::uint64_t generate_cells_requests = 0;
        std::uint64_t chunks_materialized     = 0;
        std::uint64_t cell_slots_scanned      = 0;
        std::uint64_t cell_records_inserted   = 0;

        double chunk_materialization_seconds = 0.0;
    };

    WeakScalingInstrumentation weak_scaling_instrumentation_;

    void PrintWeakScalingInstrumentation(MPI_Comm comm) const;
#endif

    struct CellMetadata {
        SInt    size;
        SInt    offset;
        LPFloat start_x;
        LPFloat start_y;
    };

    CellMetadata ReconstructCellMetadata(SInt chunk_id, SInt cell_id);

    void GenerateCellMetadataSubtree(
        SInt chunk_id, SInt begin_cell, SInt end_cell, SInt vertices, SInt offset, SInt tree_node);

    SInt SampleLeftCellCount(SInt chunk_id, SInt tree_node, SInt vertices, SInt left_cells, SInt total_cells);

    void GenerateRemoteCellVertices(SInt chunk_id, SInt cell_id, std::vector<Vertex>& buffer);
    void GenerateCells(SInt chunk_id) override;

    void GenerateEdges(SInt chunk_row, SInt chunk_column) override;

    SInt SafeTotalCellsPerDim() const;

    std::string MakeDebugFilename() const;

    struct LocalMemoryStats {
        SInt max_pins_per_hyperedge   = 0;
        SInt max_ranges_per_hyperedge = 0;
    };

    LocalMemoryStats                              local_memory_stats_;
    std::optional<GeometricHypergraphDebugLogger> debug_logger_;

    void ObserveHyperedgeAndMaybeReserve(std::size_t pins, std::size_t ranges);

    void PushHyperedgeCompressed(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges);

    struct CellPosition {
        SInt cell_id;
        SInt global_cell_x;
        SInt global_cell_y;
        SInt global_cell_id;
    };

    CellPosition
    GetCellPosition(SInt chunk_row, SInt chunk_column, SInt cell_row, SInt cell_column, SInt total_cells_per_dim) const;

    struct ExactDebugStats {
        SInt partial_cells     = 0;
        SInt vertices_tested   = 0;
        SInt vertices_accepted = 0;
    };

    mutable ExactDebugStats exact_debug_stats_;

    class CenterSampler {
    public:
        CenterSampler(HyperRGG2D& gen, LPFloat cell_width, LPFloat lower_bound, LPFloat upper_bound);

        Center Sample(SInt chunk_id, const CellPosition& cell, SInt emitted, SInt base_m, SInt remainder);

    private:
        HyperRGG2D& gen_;
        LPFloat     cell_width_;
        LPFloat     lower_bound_;
        LPFloat     upper_bound_;
    };
};

} // namespace kagen
