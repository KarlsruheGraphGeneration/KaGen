#pragma once

#include "kagen/context.h"
#include "kagen/edgeweight_generators/edge_weight_generator.h"
#include "kagen/kagen.h"

namespace kagen {
class VoidingEdgeWeightGenerator : public EdgeWeightGenerator {
public:
    VoidingEdgeWeightGenerator(EdgeWeightConfig) {}

    void GenerateEdgeWeights(const Edgelist&, EdgeWeights& weights) final {
        weights.clear();
    }

    void GenerateEdgeWeights(const XadjArray&, const AdjncyArray&, EdgeWeights& weights) final {
        weights.clear();
    }
};
} // namespace kagen
