#pragma once

#include <mpi.h>

#include "kagen/definitions.h"
#include "kagen/io/graph_writer.h"

namespace kagen {
class EdgeListWriter : public SequentialGraphWriter {
public:
    EdgeListWriter(Graph& graph, MPI_Comm comm);

    std::string DefaultExtension() const final;

protected:
    void AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) final;

    void AppendTo(const std::string& filename) final;
};

class BinaryEdgeListWriter : public SequentialGraphWriter {
public:
    BinaryEdgeListWriter(Graph& graph, MPI_Comm comm);

    std::string DefaultExtension() const final;

protected:
    void AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) final;

    void AppendTo(const std::string& filename) final;
};
} // namespace kagen
