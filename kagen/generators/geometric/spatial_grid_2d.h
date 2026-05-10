#pragma once

#include "kagen/context.h"
#include "kagen/generators/geometric/geometric_2d.h"

namespace kagen {

class SpatialGrid2D : protected virtual Geometric2D {
protected:
    SpatialGrid2D(const PGeneratorConfig& config, PEID rank, PEID size);

    void InitSpatialGrid2D(const PGeneratorConfig& config);

    SInt EncodeCell(SInt x, SInt y) const;
    void DecodeCell(SInt id, SInt& x, SInt& y) const;

protected:
    LPFloat target_r_;
};

} // namespace kagen