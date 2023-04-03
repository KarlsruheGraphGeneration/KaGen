#pragma once

#include <mpi.h>

#include "kagen/definitions.h"

namespace kagen {
void AddReverseEdges(EdgeList& edge_list, VertexRange vertex_range, MPI_Comm comm);

void AddReverseEdgesAndRedistribute(EdgeList& edge_list, VertexRange vertex_range, MPI_Comm comm);
} // namespace kagen
