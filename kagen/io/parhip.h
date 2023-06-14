#pragma once

#include <fstream>
#include <string>

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/io/graph_format.h"

namespace kagen {
class ParhipWriter : public GraphWriter {
public:
    ParhipWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm);

    void Write(bool report_progress) final;
};

class ParhipReader : public GraphReader {
public:
    ParhipReader(const std::string& filename);

    GraphSize ReadSize() final;

    Graph Read(SInt from, SInt to_node, SInt to_edge, GraphRepresentation representation) final;

    SInt FindNodeByEdge(SInt edge) final;

private:
    std::ifstream in_;

    SInt n_       = 0;
    SInt m_       = 0;
    SInt version_ = 0;
};

class ParhipFactory : public FileFormatFactory {
public:
    std::string DefaultExtension() const final {
        return "parhip";
    }

    std::unique_ptr<GraphReader> CreateReader(const InputGraphConfig& config, PEID rank, PEID size) const final;
    std::unique_ptr<GraphWriter> CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const final;
};
} // namespace kagen
