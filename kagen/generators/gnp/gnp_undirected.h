/*******************************************************************************
 * include/generators/gnp/gnp_undirected.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/kagen.h"
#include "kagen/tools/rng_wrapper.h"

namespace kagen {
class GNPUndirectedFactory : public GeneratorFactory {
public:
    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const override;
};

class GNPUndirected : public virtual Generator, private EdgeListOnlyGenerator {
public:
    GNPUndirected(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void GenerateEdgeList() final;

private:
    // Config
    const PGeneratorConfig& config_;

    PEID rank_;
    PEID size_;

    // Variates
    RNGWrapper<> rng_;

    // Constants and variables
    SInt nodes_per_chunk;
    SInt start_node_, end_node_, num_nodes_;

    void GenerateTriangleChunk(
        SInt row_id, SInt column_id, SInt row_node_id, SInt column_node_id, SInt row_n, SInt column_n);

    void GenerateRectangleChunk(
        SInt row_id, SInt column_id, SInt row_node_id, SInt column_node_id, SInt row_n, SInt column_n);

    void GenerateTriangularEdges(
        SInt row_n, SInt column_n, double p, SInt row_id, SInt column_id, SInt offset_row, SInt offset_column);

    void GenerateRectangleEdges(
        SInt row_n, SInt column_n, double p, SInt row_id, SInt column_id, SInt offset_row, SInt offset_column);
};
} // namespace kagen
