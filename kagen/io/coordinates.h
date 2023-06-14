#pragma once

#include <string>

#include <mpi.h>

#include "kagen/io/graph_format.h"
#include "kagen/io/seq_graph_writer.h"

namespace kagen {
class CoordinatesWriter : public SequentialGraphWriter {
public:
    CoordinatesWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm);

protected:
    int Requirements() const final;

    void AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) final;

    void AppendTo(const std::string& filename) final;
};

class CoordinatesFactory : public FileFormatFactory {
public:
    std::string DefaultExtension() const final {
        return "xyz";
    }

    std::unique_ptr<GraphWriter> CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const final;
};
} // namespace kagen
