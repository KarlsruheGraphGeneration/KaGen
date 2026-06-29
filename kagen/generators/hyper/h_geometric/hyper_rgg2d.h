#pragma once

#include "kagen/context.h"
#include "kagen/generators/geometric/spatial_grid_2d.h"
#include "kagen/generators/hyper/h_geometric/h_rgg.h"

#include <vector>

namespace kagen {

class HyperRGG2DPolicy;

class HyperRGG2D : public HRGG, protected SpatialGrid2D, private CSROnlyGenerator {
    friend class HyperRGG2DPolicy;

public:
    HyperRGG2D(const PGeneratorConfig& config, PEID rank, PEID size);

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
    void GenerateEdges(SInt chunk_row, SInt chunk_column) override;

    SInt SafeTotalCellsPerDim() const;

    void PushHyperedgeCompressed(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges);

    struct ExactDebugStats {
        SInt partial_cells     = 0;
        SInt vertices_tested   = 0;
        SInt vertices_accepted = 0;
    };

    mutable ExactDebugStats exact_debug_stats_;
};

} // namespace kagen
