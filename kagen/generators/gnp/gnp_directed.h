/*******************************************************************************
 * include/generators/gnp/gnp_directed.h
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
class GNPDirected : public Generator {
public:
    GNPDirected(const PGeneratorConfig& config, PEID rank, PEID size);

    GeneratorRequirement Requirements() const final;

    GeneratorFeature Features() const final;

protected:
    void GenerateImpl() final;

private:
    // Config
    const PGeneratorConfig& config_;

    PEID rank_;
    PEID size_;

    // Variates
    RNGWrapper rng_;

    // Constants and variables
    SInt edges_per_node;
    SInt start_node_, end_node_, num_nodes_;

    void GenerateChunk(SInt chunk_id, SInt node_id, SInt n);

    void GenerateEdges(SInt n, double p, SInt chunk_id, SInt offset);
};
} // namespace kagen
