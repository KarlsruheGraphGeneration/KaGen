#pragma once

#include <mpi.h>

#include "kagen/io/graph_writer.h"

namespace kagen {
class DotWriter : public SequentialGraphWriter {
public:
    DotWriter(Graph& graph, bool directed_output, MPI_Comm comm);

    std::string DefaultExtension() const final;

protected:
    int Requirements() const final;

    void AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) final;

    void AppendTo(const std::string& filename) final;

    void AppendFooterTo(const std::string& filename) final;

private:
    bool directed_output_;
};
} // namespace kagen
