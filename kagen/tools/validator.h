#pragma once

#include <mpi.h>

#include "kagen/definitions.h"
#include "kagen/kagen.h"

namespace kagen {
bool ValidateVertexRanges(const EdgeList& edge_list, VertexRange vertex_range, MPI_Comm comm);

bool ValidateGraph(
    Graph& graph, bool allow_self_loops, bool allow_directed_graphs, bool allow_multi_edges, MPI_Comm comm);
} // namespace kagen
