#pragma once

#include <mpi.h>

#include "kagen/definitions.h"
#include "kagen/io/graph_writer.h"

namespace kagen {
class EdgeListWriter : public SequentialGraphWriter {
public:
    EdgeListWriter(EdgeList& edges, VertexRange vertex_range, Coordinates& coordinates, MPI_Comm comm);

    std::string DefaultExtension() const final;

protected:
    void AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) final;

    void AppendTo(const std::string& filename);
};

class BinaryEdgeListWriter : public SequentialGraphWriter {
public:
    BinaryEdgeListWriter(EdgeList& edges, VertexRange vertex_range, Coordinates& coordinates, MPI_Comm comm);

    std::string DefaultExtension() const final;

protected:
    void AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) final;

    void AppendTo(const std::string& filename);
};
} // namespace kagen
