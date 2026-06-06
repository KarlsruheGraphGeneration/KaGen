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
    const LPFloat grid_radius = config.random_radius ? static_cast<LPFloat>(QuantileOrConstantHyperedgeRadius(config))
                                                     : static_cast<LPFloat>(config.r);
    if (grid_radius <= 0.0 || !std::isfinite(grid_radius)) {
        throw ConfigurationError("SpatialGrid2D requires positive finite grid radius");
    }

    if (!config.random_radius) {
        target_r_ = config.r * config.r;
    } else {
        target_r_ = grid_radius * grid_radius; // harmless if unused
    }

    total_chunks_ = config.k;

    chunks_per_dim_ = static_cast<SInt>(std::sqrt(total_chunks_));
    chunk_size_     = 1.0 / static_cast<LPFloat>(chunks_per_dim_);

    cells_per_dim_   = static_cast<SInt>(std::floor(chunk_size_ / grid_radius));
    cells_per_dim_   = std::max<SInt>(1, cells_per_dim_);
    cells_per_chunk_ = cells_per_dim_ * cells_per_dim_;
    cell_size_       = chunk_size_ / static_cast<LPFloat>(cells_per_dim_);
}

SInt SpatialGrid2D::EncodeCell(const SInt x, const SInt y) const {
    return x * cells_per_dim_ + y;
}

void SpatialGrid2D::DecodeCell(const SInt id, SInt& x, SInt& y) const {
    x = id / cells_per_dim_;
    y = id % cells_per_dim_;
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