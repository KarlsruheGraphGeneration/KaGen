#pragma once

#include <mpi.h>

#include "kagen/definitions.h"

namespace kagen {
bool ValidateVertexRanges(const EdgeList& edge_list, VertexRange vertex_range, MPI_Comm comm);

bool ValidateSimpleGraph(
    const EdgeList& edge_list, VertexRange vertex_range, const VertexWeights& vertex_weights,
    const EdgeWeights& edge_weights, MPI_Comm comm);
} // namespace kagen
