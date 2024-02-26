#pragma once

#include "kagen/kagen.h"

#include <mpi.h>

#include <vector>

namespace kagen {
SInt FindNumberOfGlobalNodes(VertexRange vertex_range, MPI_Comm comm);

SInt FindNumberOfGlobalEdges(const Edgelist& edges, MPI_Comm comm);

std::vector<SInt> GatherNumberOfEdges(const Edgelist& edges, MPI_Comm comm);

SInt    ReduceSum(SInt value, MPI_Comm comm);
SInt    ReduceMin(SInt value, MPI_Comm comm);
LPFloat ReduceMean(SInt value, MPI_Comm comm);
SInt    ReduceMax(SInt value, MPI_Comm comm);
LPFloat ReduceSD(SInt value, MPI_Comm comm);

struct DegreeStatistics {
    SInt    min;
    LPFloat mean;
    SInt    max;
};

DegreeStatistics ReduceDegreeStatistics(const Edgelist& edges, SInt global_num_nodes, MPI_Comm comm);

std::vector<SInt> ComputeDegreeBins(const Edgelist& edges, VertexRange vertex_range, MPI_Comm comm);

double ComputeEdgeLocalicty(const Edgelist& edges, VertexRange vertex_range, MPI_Comm comm);

SInt ComputeNumberOfGhostNodes(const Edgelist& edges, VertexRange vertex_range, MPI_Comm comm);

void PrintBasicStatistics(
    const XadjArray& xadj, const AdjncyArray& adjncy, VertexRange vertex_range, bool root, MPI_Comm comm);

void PrintBasicStatistics(const Edgelist& edges, VertexRange vertex_range, bool root, MPI_Comm comm);

void PrintAdvancedStatistics(Edgelist& edges, VertexRange vertex_range, bool root, MPI_Comm comm);
} // namespace kagen
