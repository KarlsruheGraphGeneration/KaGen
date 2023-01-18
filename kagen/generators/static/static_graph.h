#pragma once

#include "kagen/generators/generator.h"

namespace kagen {
class StaticGraphFactory : public GeneratorFactory {
public:
    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, PEID size, bool output) const override;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const override;
};

class StaticGraph : public Generator {
public:
    StaticGraph(const PGeneratorConfig& config, const PEID rank, const PEID size);

protected:
    void GenerateEdgeList() final;

    void GenerateCSR() final;

private:
    void GenerateImpl(GraphRepresentation representation);

    const PGeneratorConfig& config_;

    PEID rank_;
    PEID size_;
};
} // namespace kagen
