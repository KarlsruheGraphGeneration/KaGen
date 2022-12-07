#pragma once

#include <mpi.h>

#include "kagen/io/graph_writer.h"

namespace kagen {
class MetisWriter : public SequentialGraphWriter {
public:
    MetisWriter(Graph& graph, MPI_Comm comm);

    std::string DefaultExtension() const final;

protected:
    void AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) final;

    void AppendTo(const std::string& filename) final;

    Requirement Requirements() const final;
};
} // namespace kagen
