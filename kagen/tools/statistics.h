#pragma once

#include "kagen/comm/comm.h"
#include "kagen/kagen.h"

#include <vector>

namespace kagen {
SInt FindNumberOfGlobalNodes(VertexRange vertex_range, Comm& comm);

SInt FindNumberOfGlobalEdges(const Edgelist& edges, Comm& comm);

std::vector<SInt> GatherNumberOfEdges(const Edgelist& edges, Comm& comm);

SInt    ReduceSum(SInt value, Comm& comm);
SInt    ReduceMin(SInt value, Comm& comm);
LPFloat ReduceMean(SInt value, Comm& comm);
SInt    ReduceMax(SInt value, Comm& comm);
LPFloat ReduceSD(SInt value, Comm& comm);

struct DegreeStatistics {
    SInt    min;
    LPFloat mean;
    SInt    max;
};

DegreeStatistics ReduceDegreeStatistics(const Edgelist& edges, SInt global_num_nodes, Comm& comm);

std::vector<SInt> ComputeDegreeBins(const Edgelist& edges, VertexRange vertex_range, Comm& comm);

double ComputeEdgeLocality(const Edgelist& edges, VertexRange vertex_range, Comm& comm);

SInt ComputeNumberOfGhostNodes(const Edgelist& edges, VertexRange vertex_range, Comm& comm);

void PrintBasicStatistics(
    const XadjArray& xadj, const AdjncyArray& adjncy, VertexRange vertex_range, bool root, Comm& comm);

void PrintBasicStatistics(const Edgelist& edges, VertexRange vertex_range, bool root, Comm& comm);

void PrintAdvancedStatistics(Edgelist& edges, VertexRange vertex_range, bool root, Comm& comm);
} // namespace kagen
