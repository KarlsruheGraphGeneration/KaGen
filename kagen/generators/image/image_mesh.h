#pragma once

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/generator.h"

namespace kagen {
class ImageMeshFactory : public GeneratorFactory {
public:
    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, PEID size, bool output) const override;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const override;
};

class ImageMesh : public Generator {
public:
    ImageMesh(const PGeneratorConfig& config, const PEID rank, const PEID size);

protected:
    void GenerateImpl() final;

private:
    const PGeneratorConfig& config_;

    PEID rank_;
    PEID size_;
};
} // namespace kagen
