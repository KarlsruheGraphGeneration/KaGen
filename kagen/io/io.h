#pragma once

#include <string>

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/definitions.h"

namespace kagen {
void WriteGraph(const PGeneratorConfig& config, EdgeList& edges, VertexRange vertex_range, MPI_Comm comm);

void WriteEdgeList(
    const std::string& filename, bool omit_header, bool single_file, const EdgeList& edges, VertexRange vertex_range,
    MPI_Comm comm);

void WriteBinaryEdgeList(
    const std::string& filename, bool omit_header, bool single_file, const EdgeList& edges, VertexRange vertex_range,
    MPI_Comm comm);

void WriteMetis(const std::string& filename, EdgeList& edges, VertexRange vertex_range, MPI_Comm comm);

void WriteHMetis(const std::string& filename, EdgeList& edges, VertexRange vertex_range, MPI_Comm comm);
} // namespace kagen
