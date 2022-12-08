#pragma once

#include <string>

#include <mpi.h>

#include "kagen/io/graph_writer.h"

namespace kagen {
class BinaryParHipWriter : public GraphWriter {
public:
    BinaryParHipWriter(Graph& graph, MPI_Comm comm);

    std::string DefaultExtension() const final;

    void Write(const PGeneratorConfig& config) final;
};
} // namespace kagen
