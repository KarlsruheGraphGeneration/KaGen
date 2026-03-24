#pragma once

#include "kagen/comm/comm.h"
#include "kagen/kagen.h"

namespace kagen {
bool ValidateVertexRanges(const Edgelist& edge_list, VertexRange vertex_range, Comm& comm);

bool ValidateGraph(
    Graph& graph, bool allow_self_loops, bool allow_directed_graphs, bool allow_multi_edges, Comm& comm);

bool ValidateGraphInplace(
    Graph& graph, bool allow_self_loops, bool allow_directed_graphs, bool allow_multi_edges, Comm& comm);
} // namespace kagen
