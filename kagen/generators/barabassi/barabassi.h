/*******************************************************************************
 * include/generators/barabassi/barabassi.h
 *
 * Copyright (C) 2016-2017 Christian Schulz <christian.schulz@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/kagen.h"

#include <mpi.h>

namespace kagen {
class BarabassiFactory : public GeneratorFactory {
public:
    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, PEID rank, PEID size, bool output) const final;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const final;
};

class Barabassi : public virtual Generator, private EdgeListOnlyGenerator {
public:
    Barabassi(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void GenerateEdgeList() final;

    void FinalizeEdgeList(MPI_Comm comm) final;

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
