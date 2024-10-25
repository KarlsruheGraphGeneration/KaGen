#pragma once

#include "kagen/kagen.h"

namespace kagen {
class EdgeWeightGenerator {
public:
    virtual ~EdgeWeightGenerator() = default;

    virtual void GenerateEdgeWeights(const Edgelist& edgelist, EdgeWeights& weights) = 0;

    virtual void GenerateEdgeWeights(const XadjArray& xadj, const AdjncyArray& adjncy, EdgeWeights& weights) = 0;
};
} // namespace kagen
