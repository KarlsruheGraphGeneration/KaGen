#pragma once

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/generators/graph500_generator.h"

namespace kagen {
class RMATFactory : public GeneratorFactory {
public:
    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, PEID rank, PEID size, bool output) const final;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const final;
};

class RMAT : public Graph500Generator {
public:
    RMAT(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void GenerateEdgeList() final;

private:
    const PGeneratorConfig& config_;
    PEID                    rank_;
    SInt                    num_edges_;
};
} // namespace kagen
