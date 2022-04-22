/*******************************************************************************
 * include/generators/gnm/gnm_undirected.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include "hash.hpp"
#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/generator.h"
#include "kagen/tools/rng_wrapper.h"

namespace kagen {
class GNMUndirected : public Generator {
public:
    GNMUndirected(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void GenerateImpl() final;

private:
    // Config
    const PGeneratorConfig& config_;

    PEID rank_;
    PEID size_;

    // Globals
    SInt leftover_chunks_, nodes_per_chunk_, remaining_nodes_;
    SInt start_node_, end_node_, num_nodes_;

    // Variates
    RNGWrapper rng_;

    void GenerateChunks(SInt row);

    void QueryTriangular(
        SInt m, SInt num_rows, SInt num_columns, SInt row_id, SInt column_id, SInt offset_row, SInt offset_column,
        SInt level);

    void QueryRowRectangle(
        SInt m, SInt num_rows, SInt num_columns, SInt row_id, SInt offset_row, SInt offset_column, SInt level);

    void QueryColumnRectangle(
        SInt m, SInt num_rows, SInt num_columns, SInt column_id, SInt offset_row, SInt offset_column, SInt level);

    void GenerateTriangularEdges(SInt m, SInt row_id, SInt column_id);

    void GenerateRectangleEdges(SInt m, SInt row_id, SInt column_id);

    inline SInt NodesInRows(SInt rows, SInt offset) const;

    inline SInt NodesInColumns(SInt columns, SInt offset) const;

    inline SInt NodesInRow(SInt row) const;

    inline SInt NodesInColumn(SInt column) const;

    inline SInt OffsetInRow(SInt row) const;

    inline SInt OffsetInColumn(SInt column) const;

    inline SInt ChunkStart(SInt row, SInt column) const;

    inline HPFloat NumTriangleEdges(HPFloat row, HPFloat column, bool loops = false) const;

    inline HPFloat NumRectangleEdges(HPFloat row, HPFloat column) const;
};
} // namespace kagen