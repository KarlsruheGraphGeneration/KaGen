#pragma once

#include "kagen/context.h"
#include "kagen/edgeweight_generators/edge_weight_generator.h"
#include "kagen/edgeweight_generators/per_edge_weight_generator.h"
#include "kagen/kagen.h"

namespace kagen {
class NoneEdgeWeightGenerator : public EdgeWeightGenerator {
public:
    NoneEdgeWeightGenerator(EdgeWeightConfig) {}
    EdgeWeights GenerateEdgeWeights(const Edgelist&) final {
        return {};
    }

    EdgeWeights GenerateEdgeWeights(const XadjArray&, const AdjncyArray&) final {
        return {};
    }
};
} // namespace kagen
