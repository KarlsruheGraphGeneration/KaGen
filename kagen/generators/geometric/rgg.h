#pragma once

#include "kagen/generators/generator.h"
#include "kagen/generators/geometric/rgg/rgg_2d.h"
#include "kagen/generators/geometric/rgg/rgg_3d.h"

namespace kagen {
class RGG2DFactory : public GeneratorFactory {
public:
    int Requirements() const override;

    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, bool output) const override;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const override;
};

class RGG3DFactory : public GeneratorFactory {
public:
    int Requirements() const override;

    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, bool output) const override;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const override;
};
} // namespace kagen
