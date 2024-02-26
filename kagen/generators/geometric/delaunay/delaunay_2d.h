/*******************************************************************************
 * include/generators/geometric/delaunay/delaunay_2d.h
 *
 * Copyright (C) 2017 Daniel Funke <funke@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include "kagen/generators/geometric/geometric_2d.h"
#include "kagen/kagen.h"

#include <climits>

namespace kagen {
class Delaunay2D : public Geometric2D {
public:
    Delaunay2D(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    static constexpr SInt COPY_FLAG = SInt(1) << (sizeof(SInt) * CHAR_BIT - 1);

    SInt max_radius_;

    void GenerateEdges(SInt chunk_row, SInt chunk_column) override;

private:
    void SortCellVertices(std::vector<Vertex>& vertices);
};
} // namespace kagen
