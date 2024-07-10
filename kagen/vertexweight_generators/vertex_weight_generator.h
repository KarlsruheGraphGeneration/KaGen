#pragma once

#include "kagen/context.h"

namespace kagen {
class VertexWeightGenerator {
public:
    virtual ~VertexWeightGenerator() = default;

    virtual void
    GenerateVertexWeights(const VertexRange& vertex_range, const Edgelist& edgelist, VertexWeights& weights) = 0;
    virtual void GenerateVertexWeights(
        const VertexRange& vertex_range, const XadjArray& xadj, const AdjncyArray& adjncy, VertexWeights& weights) = 0;
};
} // namespace kagen
