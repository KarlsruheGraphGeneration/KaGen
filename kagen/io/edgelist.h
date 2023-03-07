#pragma once

#include <mpi.h>

#include "kagen/definitions.h"
#include "kagen/io/graph_writer.h"

namespace kagen {
class EdgeListWriter : public SequentialGraphWriter {
public:
    EdgeListWriter(Graph& graph, MPI_Comm comm, bool header, bool undirected);

    std::string DefaultExtension() const final;

protected:
    int Requirements() const final;

    void AppendHeaderTo(const std::string& filename, SInt n, SInt m) final;

    void AppendTo(const std::string& filename) final;

private:
    bool header_;
    bool undirected_;
};

class BinaryEdgeListWriter : public SequentialGraphWriter {
public:
    BinaryEdgeListWriter(Graph& graph, MPI_Comm comm, int width, bool header, bool undirected);

    std::string DefaultExtension() const final;

protected:
    int Requirements() const final;

    void AppendHeaderTo(const std::string& filename, SInt n, SInt m) final;

    void AppendTo(const std::string& filename) final;

private:
    int  width_;
    bool header_;
    bool undirected_;
};
} // namespace kagen
