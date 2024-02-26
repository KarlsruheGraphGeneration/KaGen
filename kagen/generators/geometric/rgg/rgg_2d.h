/*******************************************************************************
 * include/generators/geometric/rgg/rgg_2d.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Daniel Funke <funke@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include "kagen/generators/geometric/geometric_2d.h"

namespace kagen {
class RGG2D : public Geometric2D {
public:
    RGG2D(const PGeneratorConfig& config, PEID rank, PEID size);

private:
    LPFloat target_r_;

    void GenerateEdges(SInt chunk_row, SInt chunk_column) override;

    void GenerateGridEdges(SInt first_chunk_id, SInt first_cell_id, SInt second_chunk_id, SInt second_cell_id);

    void GenerateCells(SInt chunk_id) override;

    bool IsAdjacentCell(SInt chunk_id, SInt cell_id);

    SInt EncodeCell(SInt x, SInt y) const;

    void DecodeCell(SInt id, SInt& x, SInt& y) const;
};
} // namespace kagen
