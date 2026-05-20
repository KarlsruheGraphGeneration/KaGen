#pragma once

#include "kagen/generators/geometric/geometric_2d.h"

namespace kagen {

class SpatialGrid2D : public Geometric2D {
public:
    SpatialGrid2D(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void InitSpatialGrid2D(const PGeneratorConfig& config);

    SInt EncodeCell(SInt x, SInt y) const;
    void DecodeCell(SInt id, SInt& x, SInt& y) const;

    SInt SafeTotalCellsPerDim() const;

protected:
    LPFloat target_r_;
};

} // namespace kagen