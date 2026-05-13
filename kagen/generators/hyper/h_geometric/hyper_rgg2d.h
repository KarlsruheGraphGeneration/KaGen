#pragma once

#include "kagen/context.h"
#include "kagen/generators/geometric/spatial_grid_2d.h"
#include "kagen/generators/hyper/h_geometric/h_rgg.h"

namespace kagen {

class HyperRGG2D : public HRGG, protected SpatialGrid2D, private CSROnlyGenerator {
public:
    HyperRGG2D(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void GenerateCSR() override;
    void FinalizeCSR(MPI_Comm comm) override;

private:
    void GenerateEdges(SInt chunk_row, SInt chunk_column) override;
    void GenerateBallHyperedge(const Vertex& center, SInt center_chunk_id, SInt center_cell_id);
};

} // namespace kagen