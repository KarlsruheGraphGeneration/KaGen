#pragma once

#include <mpi.h>

#include "kagen/io/graph_writer.h"

namespace kagen {
class DotWriter : public SequentialGraphWriter {
public:
    DotWriter(Graph& graph, MPI_Comm comm, bool directed);

    std::string DefaultExtension() const final;

protected:
    int Requirements() const final;

    void AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) final;

    void AppendTo(const std::string& filename) final;

    void AppendFooterTo(const std::string& filename) final;

private:
    bool directed_;
};
} // namespace kagen
