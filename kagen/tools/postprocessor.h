#pragma once

#include "kagen/comm/comm.h"
#include "kagen/kagen.h"

namespace kagen {
void AddNonlocalReverseEdges(Edgelist& edge_list, EdgeWeights& edge_weights, VertexRange vertex_range, Comm& comm);

void RedistributeEdgesByVertexRange(
    Edgelist& edge_list, VertexRange vertex_range, Comm& comm, bool use_binary_search = false);

VertexRange RedistributeEdges(Edgelist& source, Edgelist& destination, SInt n, bool remap_round_robin, Comm& comm);

std::vector<SInt> ComputeBalancedVertexDistribution(SInt n, Comm& comm);

std::vector<SInt> RoundRobinRemapping(Edgelist& edges, SInt n, Comm& comm);

VertexRange
RedistributeEdgesBalanced(Edgelist& source, Edgelist& destination, SInt n, bool remap_round_robin, Comm& comm);
} // namespace kagen
