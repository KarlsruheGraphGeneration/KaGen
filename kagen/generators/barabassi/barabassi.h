/*******************************************************************************
 * include/generators/barabassi/barabassi.h
 *
 * Copyright (C) 2016-2017 Christian Schulz <christian.schulz@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/generator.h"

namespace kagen {
class Barabassi : public Generator {
public:
    Barabassi(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void GenerateImpl() final;

private:
    // Config
    const PGeneratorConfig& config_;

    // Constants and variables
    SInt min_degree_;
    SInt total_degree_;
    SInt from_, to_;

    void GenerateEdges();
};
} // namespace kagen
