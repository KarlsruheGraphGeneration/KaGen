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

struct DegreeStatistics {
    SInt    min;
    LPFloat mean;
    SInt    max;
};

DegreeStatistics ReduceDegreeStatistics(const EdgeList& edges, SInt global_num_nodes);

std::vector<SInt> ComputeDegreeBins(const EdgeList& edges, VertexRange vertex_range);

double ComputeEdgeLocalicty(const EdgeList& edges, VertexRange vertex_range);

SInt ComputeNumberOfGhostNodes(const EdgeList& edges, VertexRange vertex_range);

void PrintBasicStatistics(const EdgeList& edges, VertexRange vertex_range, bool root);

void PrintAdvancedStatistics(EdgeList& edges, VertexRange vertex_range, bool root);
} // namespace kagen
