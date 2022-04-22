#pragma once

#include <vector>

#include "kagen/definitions.h"

namespace kagen {
SInt FindNumberOfGlobalNodes(VertexRange vertex_range);

SInt FindNumberOfGlobalEdges(const EdgeList& edges);

std::vector<SInt> GatherNumberOfEdges(const EdgeList& edges);

SInt    ReduceSum(SInt value);
SInt    ReduceMin(SInt value);
LPFloat ReduceMean(SInt value);
SInt    ReduceMax(SInt value);
LPFloat ReduceSD(SInt value);
} // namespace kagen
