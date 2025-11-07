#pragma once

#include "kagen/context.h"
#include "kagen/edgeweight_generators/per_edge_weight_generator.h"
#include "kagen/kagen.h"

namespace kagen {
class HashingBasedEdgeWeightGenerator : public PerEdgeWeightGenerator<HashingBasedEdgeWeightGenerator> {
public:
    HashingBasedEdgeWeightGenerator(EdgeWeightConfig config, VertexRange vertex_range);

    SSInt GenerateEdgeWeight(SInt u, SInt v);

private:
    const EdgeWeightConfig config_;
};
} // namespace kagen
