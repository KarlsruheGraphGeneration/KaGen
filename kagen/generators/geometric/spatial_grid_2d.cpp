#include "kagen/generators/geometric/spatial_grid_2d.h"

#include "kagen/hypergraph/hypergraph_utils.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace kagen {

SpatialGrid2D::SpatialGrid2D(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : Geometric2D(config, rank, size) {
    InitSpatialGrid2D(config);
    InitDatastructures();
}

void SpatialGrid2D::InitSpatialGrid2D(const PGeneratorConfig& config) {
    const LPFloat representative_radius = static_cast<LPFloat>(QuantileOrConstantHyperedgeRadius(config));

    if (representative_radius <= 0.0 || !std::isfinite(representative_radius)) {
        throw ConfigurationError("SpatialGrid2D requires a positive finite representative radius");
    }

    total_chunks_ = config.k;

    chunks_per_dim_ = static_cast<SInt>(std::sqrt(total_chunks_));

    chunk_size_ = 1.0 / static_cast<LPFloat>(chunks_per_dim_);

    const LPFloat desired_cell_size = 2.0 * representative_radius;

    cells_per_dim_ = std::max<SInt>(1, static_cast<SInt>(std::floor(chunk_size_ / desired_cell_size)));

    cells_per_chunk_ = cells_per_dim_ * cells_per_dim_;

    cell_size_ = chunk_size_ / static_cast<LPFloat>(cells_per_dim_);

    if (!config.random_radius) {
        target_r_ = config.r * config.r;
    } else {
        target_r_ = representative_radius * representative_radius;
    }
}

SInt SpatialGrid2D::EncodeCell(SInt x, SInt y) const {
    // x = column, y = row
    return (y * cells_per_dim_) + x;
}

void SpatialGrid2D::DecodeCell(SInt id, SInt& x, SInt& y) const {
    // x = column, y = row
    x = id % cells_per_dim_;
    y = id / cells_per_dim_;
}

SInt SpatialGrid2D::SafeTotalCellsPerDim() const {
    if (chunks_per_dim_ <= 0 || cells_per_dim_ <= 0) {
        throw ConfigurationError("Invalid grid: chunks_per_dim_ or cells_per_dim_ not initialized");
    }

    if (chunks_per_dim_ > std::numeric_limits<SInt>::max() / cells_per_dim_) {
        throw ConfigurationError("Invalid grid: total_cells_per_dim overflow");
    }

    return chunks_per_dim_ * cells_per_dim_;
}

} // namespace kagen