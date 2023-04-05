#pragma once

#include <mpi.h>

#include "kagen/io/graph_reader.h"
#include "kagen/io/graph_writer.h"
#include "kagen/io/mmap_toker.h"

namespace kagen {
class MetisWriter : public SequentialGraphWriter {
public:
    MetisWriter(Graph& graph, MPI_Comm comm);

    std::string DefaultExtension() const final;

protected:
    void AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) final;

    void AppendTo(const std::string& filename) final;

    int Requirements() const final;
};

class MetisReader : public GraphReader {
public:
    MetisReader(const std::string& filename);

    GraphSize ReadSize() final;

    Graph Read(SInt from, SInt to, SInt num_edges, GraphRepresentation representation) final;

    SInt FindNodeByEdge(SInt edge) final;

private:
    MappedFileToker toker_;

    SInt cached_first_node_     = 0;
    SInt cached_first_edge_     = 0;
    SInt cached_first_node_pos_ = 0;
};
} // namespace kagen
