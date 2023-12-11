/*******************************************************************************
 * include/generators/grid/grid_3d.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/generator.h"
#include "kagen/tools/rng_wrapper.h"

namespace kagen {
class Grid3DFactory : public GeneratorFactory {
public:
    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, PEID rank, PEID size, bool output) const final;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const final;
};

class Grid3D : public virtual Generator, private EdgeListOnlyGenerator {
public:
    Grid3D(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void GenerateEdgeList() final;

private:
    // Config
    const PGeneratorConfig& config_;

    // Variates
    RNGWrapper<> rng_;

    // Constants and variables
    SInt    start_node_, end_node_, num_nodes_;
    LPFloat edge_probability_;
    SInt    total_x_, total_y_, total_z_;
    SInt    total_chunks_, chunks_per_dim_;
    SInt    x_per_chunk_, y_per_chunk_, z_per_chunk_;
    SInt    remaining_x_, remaining_y_, remaining_z_;
    SInt    vertices_per_chunk_;

    PEID rank_;
    PEID size_;

    void GenerateChunk(SInt chunk);

    void GenerateEdges(SInt chunk, SInt vertex);

    void QueryInDirection(SInt chunk, SInt vertex, Direction direction);

    bool IsLocalVertex(SSInt local_x, SSInt local_y, SSInt local_z, SInt xs, SInt ys, SInt zs) const;

    bool IsValidChunk(SSInt chunk_x, SSInt chunk_y, SSInt chunk_z) const;

    SInt LocateVertexInChunk(SInt chunk, SInt local_x, SInt local_y, SInt local_z, Direction direction) const;

    void GenerateEdge(SInt source, SInt target);

    SSInt DirectionX(Direction direction) const;

    SSInt DirectionY(Direction direction) const;

    SSInt DirectionZ(Direction direction) const;

    SInt OffsetForChunk(SInt chunk) const;

    void Decode(SInt id, SInt& x, SInt& y, SInt& z) const;

    SInt Encode(SInt x, SInt y, SInt z) const;
};
} // namespace kagen
