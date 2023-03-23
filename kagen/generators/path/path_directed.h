#pragma once

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/generator.h"
#include "kagen/tools/random_permutation.h"

namespace kagen {
class PathDirectedFactory : public GeneratorFactory {
public:
    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const final;
};

class PathDirected : public virtual Generator, private EdgeListOnlyGenerator {
public:
    PathDirected(const PGeneratorConfig& config, const PEID rank, const PEID size);

protected:
    void GenerateEdgeList() final;

private:
    std::pair<SInt, bool>
    GetTargetVertex(SInt j, const random_permutation::FeistelPseudoRandomPermutation& permutator_) const;

    const PGeneratorConfig& config_;

    PEID rank_;
    PEID size_;
};
} // namespace kagen
