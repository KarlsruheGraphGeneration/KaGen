/*******************************************************************************
 * include/generators/gnm/gnm_directed.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include "hash.hpp"
#include "kagen/definitions.h"
#include "kagen/generator_config.h"
#include "kagen/generators/generator.h"
#include "kagen/tools/rng_wrapper.h"

namespace kagen {
class GNMDirected : public Generator {
public:
    GNMDirected(PGeneratorConfig& config, PEID rank, PEID size);

    GeneratorRequirement Requirements() const final;

    GeneratorFeature Features() const final;

protected:
    void GenerateImpl() final;

private:
    // Config
    PGeneratorConfig& config_;

    PEID rank_;
    PEID size_;

    // Variates
    RNGWrapper rng_;

    // Constants and variables
    SInt edges_per_node_;
    SInt start_node_, end_node_, num_nodes_;

    void GenerateChunk(SInt chunk_id);

    void GenerateChunk(SInt n, SInt m, SInt k, SInt chunk_id, SInt chunk_start, SInt node_start, SInt level);

    void GenerateEdges(SInt n, SInt m, SInt chunk_id, SInt offset);
};
} // namespace kagen
