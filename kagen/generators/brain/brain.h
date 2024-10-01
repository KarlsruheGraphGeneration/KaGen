#pragma once

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "sim/Simulation.h"

#include <Types.h>
#include <mpi.h>

namespace kagen {
class BrainFactory : public GeneratorFactory {
public:
    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, PEID rank, PEID size, bool output) const final;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const final;
};

  class BrainGenerator : public virtual Generator, private EdgeListOnlyGenerator {
public:
    BrainGenerator(const PGeneratorConfig& config, const PEID rank, const PEID size);

    void GenerateEdgeList() final;

  private:

    Simulation _sim;
    BrainConfig _config;
    PEID _rank;

};

} // namespace kagen
