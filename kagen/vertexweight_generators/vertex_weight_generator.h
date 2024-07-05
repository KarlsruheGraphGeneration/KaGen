#pragma once

#include "kagen/context.h"

namespace kagen {
class VertexWeightGenerator {
public:
    virtual ~VertexWeightGenerator() = default;

    virtual VertexWeights GenerateVertexWeights(const VertexRange& vertex_range, const Edgelist& edgelist) = 0;
    virtual VertexWeights GenerateVertexWeights(const VertexRange& vertex_range, const XadjArray& xadj, const AdjncyArray& adjncy) = 0;
};
} // namespace kagen
