#pragma once

#include "kagen/generators/generator.h"
#include "kagen/generators/geometric/delaunay/delaunay_2d.h"
#include "kagen/generators/geometric/delaunay/delaunay_3d.h"

namespace kagen {
class Delaunay2DFactory : public GeneratorFactory {
public:
    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, PEID rank, PEID size, bool output) const final;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const final;
};

class Delaunay3DFactory : public GeneratorFactory {
public:
    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, PEID rank, PEID size, bool output) const final;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const final;
};
} // namespace kagen
