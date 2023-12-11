/*******************************************************************************
 * include/generators/grid/grid_2d.h
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
class Grid2DFactory : public GeneratorFactory {
public:
    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, PEID rank, PEID size, bool output) const final;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const final;
};

class Grid2D : public virtual Generator, private EdgeListOnlyGenerator {
public:
    Grid2D(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void GenerateEdgeList() override;

private:
    // Config
    const PGeneratorConfig& config_;

    PEID rank_;
    PEID size_;

    // Variates
    RNGWrapper<> rng_;

    // Constants and variables
    SInt    start_node_, end_node_, num_nodes_;
    LPFloat edge_probability_;
    SInt    total_rows_, total_cols_;
    SInt    total_chunks_, chunks_per_dim_;
    SInt    rows_per_chunk_, cols_per_chunk_;
    SInt    remaining_rows_, remaining_cols_;
    SInt    vertices_per_chunk_;

    void GenerateChunk(SInt chunk);

    void GenerateEdges(SInt chunk, SInt vertex);

    void QueryInDirection(SInt chunk, SInt vertex, Direction direction);

    bool IsLocalVertex(SInt local_row, SInt local_col, SInt rows, SInt cols) const;

    bool IsValidChunk(SInt chunk_row, SInt chunk_col) const;

    SInt LocateVertexInChunk(SInt chunk, SInt local_row, SInt local_col, Direction direction) const;

    void GenerateEdge(SInt source, SInt target);

    SSInt DirectionRow(Direction direction) const;

    SSInt DirectionColumn(Direction direction) const;

    SInt OffsetForChunk(SInt chunk) const;

    void Decode(SInt id, SInt& x, SInt& y) const;

    SInt Encode(SInt x, SInt y) const;
};
} // namespace kagen
