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

class EdgeWeightGeneratorFactory {
public:
    virtual ~EdgeWeightGeneratorFactory() = default;

    virtual std::unique_ptr<EdgeWeightGenerator> Create(EdgeWeightConfig config) const = 0;
};
} // namespace kagen
