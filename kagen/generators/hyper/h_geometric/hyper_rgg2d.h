#pragma once

#include "kagen/context.h"
#include "kagen/generators/geometric/geometric_2d.h"
#include "kagen/generators/geometric/spatial_grid_2d.h"
#include "kagen/generators/hyper/h_geometric/h_rgg.h"
#include "kagen/hypergraph/hypergraph_utils.h"

#include <vector>

namespace kagen {

class HyperRGG2D : public HRGG, protected SpatialGrid2D, private CSROnlyGenerator {
public:
    HyperRGG2D(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void GenerateCSR() override;
    void FinalizeCSR(MPI_Comm comm) override;

private:
    void GenerateEdges(SInt chunk_row, SInt chunk_column);

    bool GenerateBallHyperedge(
        LPFloat center_x, LPFloat center_y, SInt sampled_center_id, SInt center_chunk_id, SInt center_cell_id);

    CellBallRelation
    ClassifyCellAgainstBall(LPFloat center_x, LPFloat center_y, LPFloat radius, SInt chunk_id, SInt cell_id);

    void AddWholeCellRange(SInt neighbor_chunk_id, SInt neighbor_cell_id, std::vector<PinRange>& ranges);

    void PushHyperedgeCompressed(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges);

    std::pair<std::vector<SInt>, std::vector<PinRange>>
    NormalizeHyperedge(std::vector<SInt> pins, std::vector<PinRange> ranges);
};

} // namespace kagen