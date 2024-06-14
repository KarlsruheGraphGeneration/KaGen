#pragma once

#include "kagen/context.h"

namespace kagen {
class EdgeWeightGenerator {
public:
    virtual ~EdgeWeightGenerator() = default;

    virtual EdgeWeights GenerateEdgeWeights(const Edgelist& edgelist) = 0;

    virtual EdgeWeights
    GenerateEdgeWeights(const XadjArray& xadj, const AdjncyArray& adjncy) = 0;
};
} // namespace kagen
