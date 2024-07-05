#pragma once

#include "kagen/context.h"
#include "kagen/kagen.h"
#include "kagen/vertexweight_generators/vertex_weight_generator.h"

namespace kagen {
class NoneVertexWeightGenerator : public VertexWeightGenerator {
public:
    NoneVertexWeightGenerator(VertexWeightConfig) {}
    VertexWeights GenerateVertexWeights(const VertexRange&, const Edgelist&) final {
        return {};
    }
    VertexWeights GenerateVertexWeights(const VertexRange&, const XadjArray&, const AdjncyArray&) final {
        return {};
    }
};
} // namespace kagen
