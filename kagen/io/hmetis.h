#pragma once

#include <string>

#include <mpi.h>

#include "kagen/io/graph_format.h"
#include "kagen/io/seq_graph_writer.h"

namespace kagen {
class HmetisWriter : public SequentialGraphWriter {
public:
    HmetisWriter(bool directed, const OutputGraphConfig& config, Graph& graph, MPI_Comm comm);

protected:
    void AppendHeaderTo(const std::string& filename, SInt n, SInt m) final;

    void AppendTo(const std::string& filename) final;

    void AppendFooterTo(const std::string& filename) final;

private:
    bool directed_;
};

class HmetisFactory : public FileFormatFactory {
public:
    std::string DefaultExtension() const final {
        return "hmetis";
    }

    std::unique_ptr<GraphReader> CreateReader(const InputGraphConfig& config) const final;

    std::unique_ptr<GraphWriter> CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const final;
};

class DirectedHmetisFactory : public FileFormatFactory {
public:
    std::string DefaultExtension() const final {
        return "hmetis";
    }

    std::unique_ptr<GraphReader> CreateReader(const InputGraphConfig& config) const final;

    std::unique_ptr<GraphWriter> CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const final;
};
} // namespace kagen
