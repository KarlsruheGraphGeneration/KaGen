#pragma once

#include "kagen/context.h"
#include "kagen/kagen.h"
#include "kagen/vertexweight_generators/vertex_weight_generator.h"

namespace kagen {
class DefaultVertexWeightGenerator : public VertexWeightGenerator {
public:
    DefaultVertexWeightGenerator(VertexWeightConfig) {}
    void GenerateVertexWeights(const VertexRange&, const Edgelist&, VertexWeights&) final {}

    void GenerateVertexWeights(const VertexRange&, const XadjArray&, const AdjncyArray&, VertexWeights&) final {}
};
} // namespace kagen
