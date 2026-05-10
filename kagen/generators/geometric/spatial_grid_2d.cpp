#include "kagen/generators/geometric/spatial_grid_2d.h"

#include <cmath>

namespace kagen {
SpatialGrid2D::SpatialGrid2D(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : Geometric2D(config, rank, size) {
    InitSpatialGrid2D(config);
    InitDatastructures();
}

void SpatialGrid2D::InitSpatialGrid2D(const PGeneratorConfig& config) {
    target_r_ = config.r * config.r;

    total_chunks_ = config.k;

    chunks_per_dim_ = static_cast<SInt>(std::sqrt(total_chunks_));
    chunk_size_     = 1.0 / static_cast<LPFloat>(chunks_per_dim_);

    cells_per_dim_   = static_cast<SInt>(std::floor(chunk_size_ / config.r));
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

} // namespace kagen