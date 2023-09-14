#pragma once

#include "kagen/io/graph_format.h"
#include "kagen/io/mmap_toker.h"

#include <mpi.h>

namespace kagen {
class MetisWriter : public StandardGraphWriter {
public:
    MetisWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size);

protected:
    void WriteHeader(const std::string& filename, const SInt n, const SInt m) final;

    bool WriteBody(const std::string& filename) final;
};

class MetisReader : public GraphReader {
public:
    MetisReader(const std::string& filename);

    GraphSize ReadSize() final;

    Graph Read(SInt from_vertex, SInt to_vertex, SInt to_edge, GraphRepresentation representation) final;

    SInt FindNodeByEdge(SInt edge) final;

private:
    MappedFileToker toker_;

    SInt cached_first_vertex_     = 0;
    SInt cached_first_edge_       = 0;
    SInt cached_first_vertex_pos_ = 0;
};

class MetisFactory : public FileFormatFactory {
public:
    std::vector<std::string> DefaultExtensions() const final {
        return {"metis", "graph"}; // Keep *.graph for legacy reasons
    }

    std::unique_ptr<GraphReader> CreateReader(const InputGraphConfig& config, PEID rank, PEID size) const final;

    std::unique_ptr<GraphWriter>
    CreateWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size) const final;
};
} // namespace kagen
