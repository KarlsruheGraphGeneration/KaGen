#pragma once

#include "kagen/context.h"
#include "kagen/edgeweight_generators/edge_weight_generator.h"
#include "kagen/kagen.h"

namespace kagen {
class NoneEdgeWeightGenerator : public EdgeWeightGenerator {
public:
    NoneEdgeWeightGenerator(EdgeWeightConfig) {}

    void GenerateEdgeWeights(const Edgelist&, EdgeWeights&) final {}

    void GenerateEdgeWeights(const XadjArray&, const AdjncyArray&, EdgeWeights&) final {}
};
} // namespace kagen
