#pragma once

#include "kagen/context.h"
#include "kagen/edgeweight_generators/edge_weight_generator.h"
#include "kagen/edgeweight_generators/per_edge_weight_generator.h"
#include "kagen/kagen.h"

namespace kagen {
class HashingBasedEdgeWeightGenerator : public PerEdgeWeightGenerator<HashingBasedEdgeWeightGenerator> {
public:
    HashingBasedEdgeWeightGenerator(EdgeWeightConfig config);

    SSInt GenerateEdgeWeight(SInt u, SInt v);

private:
    const EdgeWeightConfig config_;
};

class HashingBasedEdgeWeightGeneratorFactory : public EdgeWeightGeneratorFactory {
public:
    std::unique_ptr<EdgeWeightGenerator> Create(EdgeWeightConfig config) const final;
};
} // namespace kagen
