#pragma once

#include "kagen/context.h"
#include "kagen/kagen.h"
#include "kagen/vertexweight_generators/vertex_weight_generator.h"

namespace kagen {
class VoidingVertexWeightGenerator : public VertexWeightGenerator {
public:
    VoidingVertexWeightGenerator(VertexWeightConfig) {}
    void GenerateVertexWeights(const VertexRange&, const Edgelist&, VertexWeights& weights) final {
        weights.clear();
    }
    void GenerateVertexWeights(const VertexRange&, const XadjArray&, const AdjncyArray&, VertexWeights& weights) final {
        weights.clear();
    }
};
} // namespace kagen
