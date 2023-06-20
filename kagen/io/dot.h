#pragma once

#include <string>

#include <mpi.h>

#include "kagen/io/graph_format.h"

namespace kagen {
class DotWriter : public SequentialGraphWriter {
public:
    DotWriter(const bool directed, const OutputGraphConfig& config, Graph& graph, MPI_Comm comm);

protected:
    int Requirements() const final;

    void AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) final;

    void AppendTo(const std::string& filename) final;

    void AppendFooterTo(const std::string& filename) final;

private:
    bool directed_;
};

class DotFactory : public FileFormatFactory {
public:
    std::string DefaultExtension() const final {
        return "dot";
    }

    std::unique_ptr<GraphWriter> CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const final;
};

class DirectedDotFactory : public FileFormatFactory {
public:
    std::string DefaultExtension() const final {
        return "dot";
    }

    std::unique_ptr<GraphWriter> CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const final;
};
} // namespace kagen
