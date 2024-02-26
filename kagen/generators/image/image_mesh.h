#pragma once

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/kagen.h"

namespace kagen {
class ImageMeshFactory : public GeneratorFactory {
public:
    PGeneratorConfig NormalizeParameters(PGeneratorConfig config, PEID rank, PEID size, bool output) const final;

    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const final;
};

class ImageMesh : public virtual Generator, private EdgeListOnlyGenerator {
public:
    ImageMesh(const PGeneratorConfig& config, const PEID rank, const PEID size);

protected:
    void GenerateEdgeList() final;

private:
    const PGeneratorConfig& config_;

    PEID rank_;
    PEID size_;
};
} // namespace kagen
