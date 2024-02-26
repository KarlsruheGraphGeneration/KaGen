/*******************************************************************************
 * include/generators/geometric/delaunay/delaunay_3d.h
 *
 * Copyright (C) 2017 Daniel Funke <funke@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include "kagen/generators/geometric/geometric_3d.h"

#include <climits>

namespace kagen {
class Delaunay3D : public Geometric3D {
public:
    Delaunay3D(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    static constexpr SInt COPY_FLAG = SInt(1) << (sizeof(SInt) * CHAR_BIT - 1);

    SInt max_radius_;

    void GenerateEdges(SInt chunk_row, SInt chunk_column, SInt chunk_depth) override;

private:
    void SortCellVertices(std::vector<Vertex>& vertices) const;
};
} // namespace kagen
