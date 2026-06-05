#include "kagen/generators/hyper/h_geometric/hyper_rgg2d.h"

#include "kagen/generators/hyper/h_geometric/hyper_rgg2d_policy.h"
#include "kagen/hypergraph/hyperedge_builder.h"
#include "kagen/sampling/hash.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace kagen {

HyperRGG2D::HyperRGG2D(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : SpatialGrid2D(config, rank, size) {}

SInt HyperRGG2D::SafeTotalCellsPerDim() const {
    if (cells_per_dim_ <= 0 || chunks_per_dim_ <= 0) {
        throw ConfigurationError("Invalid grid: cells_per_dim_ or chunks_per_dim_ is zero");
    }

    if (chunks_per_dim_ > std::numeric_limits<SInt>::max() / cells_per_dim_) {
        throw ConfigurationError("Invalid grid: total_cells_per_dim overflow");
    }

    return chunks_per_dim_ * cells_per_dim_;
}

void HyperRGG2D::GenerateEdges(const SInt chunk_row, const SInt chunk_column) {
    const SInt chunk_id = Encode(chunk_column, chunk_row);

    if (!IsLocalChunk(chunk_id)) {
        return;
    }

    const SInt total_cells_per_dim = SafeTotalCellsPerDim();

    if (total_cells_per_dim > std::numeric_limits<SInt>::max() / total_cells_per_dim) {
        throw ConfigurationError("Invalid grid: total_cells overflow");
    }

    const SInt total_cells = total_cells_per_dim * total_cells_per_dim;

    HyperRGG2DPolicy                   policy(*this);
    HyperedgeBuilder<HyperRGG2DPolicy> builder(policy, config_);

    for (SInt cell_row = 0; cell_row < cells_per_dim_; ++cell_row) {
        for (SInt cell_column = 0; cell_column < cells_per_dim_; ++cell_column) {
            const SInt cell_id = EncodeCell(cell_column, cell_row);

            const SInt global_cell_x  = (chunk_column * cells_per_dim_) + cell_column;
            const SInt global_cell_y  = (chunk_row * cells_per_dim_) + cell_row;
            const SInt global_cell_id = (global_cell_y * total_cells_per_dim) + global_cell_x;

            const SInt    base_m     = config_.m / total_cells;
            const SInt    remainder  = config_.m % total_cells;
            const SInt    cell_m     = base_m + static_cast<SInt>(global_cell_id < remainder);
            const LPFloat cell_width = LPFloat{1.0} / static_cast<LPFloat>(total_cells_per_dim);

            for (SInt emitted = 0; emitted < cell_m; emitted++) {
                const SInt sampled_center_id =
                    (global_cell_id * base_m) + std::min(global_cell_id, remainder) + emitted;

                const SInt seed_x = sampling::Spooky::hash(config_.seed + (17 * sampled_center_id));
                const SInt seed_y = sampling::Spooky::hash(config_.seed + (31 * sampled_center_id));

                const LPFloat random_x = rng_.GenerateUniform<LPFloat>(seed_x);
                const LPFloat random_y = rng_.GenerateUniform<LPFloat>(seed_y);

                const LPFloat center_x = (static_cast<LPFloat>(global_cell_x) + random_x) * cell_width;
                const LPFloat center_y = (static_cast<LPFloat>(global_cell_y) + random_y) * cell_width;

                const HyperRGG2DPolicy::Center center{
                    .x          = center_x,
                    .y          = center_y,
                    .sampled_id = sampled_center_id,
                    .chunk_id   = chunk_id,
                    .cell_id    = cell_id,
                };

                builder.Build(center);
            }
        }
    }
}

void HyperRGG2D::PushHyperedgeCompressed(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges) {
    if (graph_.hyperedge_offsets.empty()) {
        graph_.hyperedge_offsets.push_back(0);
    }

    if (graph_.hyperedge_range_offsets.empty()) {
        graph_.hyperedge_range_offsets.push_back(0);
    }

    graph_.hyperedge_pins.insert(graph_.hyperedge_pins.end(), pins.begin(), pins.end());
    graph_.hyperedge_offsets.push_back(graph_.hyperedge_pins.size());

    graph_.hyperedge_ranges.insert(graph_.hyperedge_ranges.end(), ranges.begin(), ranges.end());
    graph_.hyperedge_range_offsets.push_back(graph_.hyperedge_ranges.size());
}

void HyperRGG2D::GenerateCSR() {
    GenerateGeometry();
}

void HyperRGG2D::FinalizeCSR(MPI_Comm /*comm*/) {}

} // namespace kagen
