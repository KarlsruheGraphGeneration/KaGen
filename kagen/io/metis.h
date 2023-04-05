#pragma once

#include <mpi.h>

#include "kagen/io/graph_format.h"
#include "kagen/io/mmap_toker.h"
#include "kagen/io/seq_graph_writer.h"

namespace kagen {
class MetisWriter : public SequentialGraphWriter {
public:
    MetisWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm);

protected:
    void AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) final;

    void AppendTo(const std::string& filename) final;

    int Requirements() const final;
};

class MetisReader : public GraphReader {
public:
    MetisReader(const std::string& filename);

    GraphSize ReadSize() final;

    Graph Read(SInt from_vertex, SInt to_vertex, SInt to_edge, GraphRepresentation representation) final;

    SInt FindNodeByEdge(SInt edge) final;

private:
    MappedFileToker toker_;

    SInt cached_first_vertex_   = 0;
    SInt cached_first_edge_     = 0;
    SInt cached_first_node_pos_ = 0;
};

class MetisFactory : public FileFormatFactory {
public:
    std::string DefaultExtension() const final {
        return "metis";
    }

    std::unique_ptr<GraphReader> CreateReader(const InputGraphConfig& config) const final;

    std::unique_ptr<GraphWriter> CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const final;
};
} // namespace kagen
