#pragma once

#include <fstream>
#include <string>

#include <mpi.h>

#include "kagen/io/graph_reader.h"
#include "kagen/io/graph_writer.h"

namespace kagen {
class BinaryParHipWriter : public GraphWriter {
public:
    BinaryParHipWriter(Graph& graph, MPI_Comm comm);

    std::string DefaultExtension() const final;

    void Write(const PGeneratorConfig& config) final;
};

class BinaryParhipReader : public GraphReader {
public:
    BinaryParhipReader(const std::string& filename);

    GraphSize ReadSize() final;

    Graph Read(SInt from, SInt to_node, SInt to_edge, GraphRepresentation representation) final;

    SInt FindNodeByEdge(SInt edge) final;

private:
    std::ifstream in_;

    SInt n_       = 0;
    SInt m_       = 0;
    SInt version_ = 0;
};
} // namespace kagen
