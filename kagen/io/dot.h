#pragma once

#include <mpi.h>

#include "kagen/io/graph_writer.h"

namespace kagen {
class DotWriter : public SequentialGraphWriter {
public:
    DotWriter(EdgeList& edges, VertexRange vertex_range, Coordinates& coordinates, MPI_Comm comm);

    std::string DefaultExtension() const final;

protected:
    void AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) final;

    void AppendTo(const std::string& filename) final;

    void AppendFooterTo(const std::string& filename) final;
};
} // namespace kagen
