#pragma once

#include "kagen/context.h"
#include "kagen/edgeweight_generators/edge_weight_generator.h"
#include "kagen/kagen.h"

namespace kagen {
class EuclideanDistanceEdgeWeightGenerator : public EdgeWeightGenerator {
public:
    EuclideanDistanceEdgeWeightGenerator(EdgeWeightConfig) {}

    void GenerateEdgeWeights(const Edgelist&, EdgeWeights&) final {}

    void GenerateEdgeWeights(const XadjArray&, const AdjncyArray&, EdgeWeights&) final {}
};
} // namespace kagen
