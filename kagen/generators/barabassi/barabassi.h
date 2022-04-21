/*******************************************************************************
 * include/generators/barabassi/barabassi.h
 *
 * Copyright (C) 2016-2017 Christian Schulz <christian.schulz@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include "kagen/definitions.h"
#include "kagen/generator_config.h"
#include "kagen/generators/generator.h"

namespace kagen {
class Barabassi : public Generator {
public:
    Barabassi(PGeneratorConfig& config, PEID rank, PEID size);

    GeneratorRequirement Requirements() const final;

    GeneratorFeature Features() const final;

protected:
    void GenerateImpl() final;

private:
    // Config
    PGeneratorConfig& config_;

    // Constants and variables
    SInt min_degree_;
    SInt total_degree_;
    SInt from_, to_;

    void GenerateEdges();
};
} // namespace kagen
