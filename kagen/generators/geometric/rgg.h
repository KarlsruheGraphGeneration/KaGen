#pragma once

#include "kagen/generators/generator.h"

#include <cmath>

namespace kagen {
class RGG2DFactory : public GeneratorFactory {
public:
    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, PEID rank, PEID size, bool output) const final;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const final;
};

class RGG3DFactory : public GeneratorFactory {
public:
    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, PEID rank, PEID size, bool output) const final;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const final;
};

class RGG : public virtual Generator {
protected:
    void
    PushWeightIfRequested(const EdgeWeightConfig& config, const double& squared_distance, const double& squared_radius);
};
} // namespace kagen
