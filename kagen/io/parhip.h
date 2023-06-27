#pragma once

#include <fstream>
#include <string>

#include "kagen/context.h"
#include "kagen/io/graph_format.h"

namespace kagen {
namespace parhip {
using ParhipID     = unsigned long long;
using ParhipWeight = SSInt;

ParhipID BuildVersion(
    bool has_vertex_weights, bool has_edge_weights, bool has_32bit_edge_ids, bool has_32bit_vertex_ids,
    bool has_32bit_vertex_weights, bool has_32bit_edge_weights);

bool HasVertexWeights(SInt version);

bool HasEdgeWeights(SInt version);

bool Has32BitEdgeIDs(SInt version);

bool Has32BitVertexIDs(SInt version);

bool Has32BitVertexWeights(SInt version);

bool Has32BitEdgeWeights(SInt version);
} // namespace parhip

class ParhipWriter : public GraphWriter {
public:
    ParhipWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size);

    bool Write(const int pass, const std::string& filename) final;

private:
    void WriteHeader(const std::string& filename);

    void WriteOffsets(const std::string& filename);

    void WriteEdges(const std::string& filename);

    void WriteVertexWeights(const std::string& filename);

    void WriteEdgeWeights(const std::string& filename);
};

class ParhipReader : public GraphReader {
public:
    ParhipReader(const InputGraphConfig& config);

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

    std::unique_ptr<GraphWriter>
    CreateWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size) const final;
};
} // namespace kagen
