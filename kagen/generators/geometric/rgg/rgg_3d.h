/*******************************************************************************
 * include/generators/geometric/rgg/rgg_3d.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Daniel Funke <funke@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include "kagen/generators/geometric/geometric_3d.h"

namespace kagen {
class RGG3D : public Geometric3D {
public:
    RGG3D(const PGeneratorConfig& config, const PEID rank, const PEID size);

protected:
    LPFloat target_r_;

    void GenerateEdges(SInt chunk_row, SInt chunk_column, SInt chunk_depth) override;

    void GenerateGridEdges(SInt first_chunk_id, SInt first_cell_id, SInt second_chunk_id, SInt second_cell_id);

    void GenerateCells(SInt chunk_id) override;

    bool IsAdjacentCell(SInt chunk_id, SInt cell_id);

    SInt EncodeCell(SInt x, SInt y, SInt z) const;

    void DecodeCell(SInt id, SInt& x, SInt& y, SInt& z) const;
};
} // namespace kagen
