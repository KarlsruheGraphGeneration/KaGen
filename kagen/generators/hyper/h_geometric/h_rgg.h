#pragma once

#include "kagen/generators/generator.h"

#include <cmath>

namespace kagen {
class HyperRGG2DFactory : public GeneratorFactory {
public:
    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, PEID rank, PEID size, bool output) const final;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const final;
};

// TODO(clickup)[2026-05-08]: Implementation of Geometric Factories

class HRGG : public virtual Generator {
protected:
    void
    PushWeightIfRequested(const EdgeWeightConfig& config, const double& squared_distance, const double& squared_radius);
};
} // namespace kagen
