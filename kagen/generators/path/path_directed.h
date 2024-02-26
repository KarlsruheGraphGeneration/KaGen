#pragma once

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/kagen.h"

namespace kagen {
class PathDirectedFactory : public GeneratorFactory {
public:
    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, PEID rank, PEID size, bool output) const final;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const final;
};

class PathDirected : public virtual Generator, private EdgeListOnlyGenerator {
public:
    PathDirected(const PGeneratorConfig& config, const PEID rank, const PEID size);

protected:
    void GenerateEdgeList() final;

private:
    const PGeneratorConfig& config_;

    PEID rank_;
    PEID size_;
};
} // namespace kagen
