#pragma once

#include <mpi.h>

#include "kagen/io/graph_writer.h"

namespace kagen {
class HMetisWriter : public SequentialGraphWriter {
public:
    HMetisWriter(Graph& graph, MPI_Comm comm, bool directed);

    std::string DefaultExtension() const final;

protected:
    void AppendHeaderTo(const std::string& filename, SInt n, SInt m) final;

    void AppendTo(const std::string& filename) final;

    void AppendFooterTo(const std::string& filename) final;

private:
    bool directed_;
};
} // namespace kagen
