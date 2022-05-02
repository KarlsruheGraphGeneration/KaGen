#pragma once

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/definitions.h"

namespace kagen {
void WriteGraph(const PGeneratorConfig& config, EdgeList& edges, VertexRange vertex_range, Coordinates &coordinates, MPI_Comm comm);

void WriteMetis(const std::string& filename, EdgeList& edges, VertexRange vertex_range, MPI_Comm comm);

void WriteHMetis(const std::string& filename, EdgeList& edges, VertexRange vertex_range, MPI_Comm comm);
} // namespace kagen
