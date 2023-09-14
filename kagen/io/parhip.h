/*******************************************************************************
 * Implements an extension to the binary graph format used by the ParHiP graph
 * partitioner.
 *
 * The format is structured as follows:
 *
 * +------------------------------------+
 * | Header (24 bytes)                  |
 * +------------------------------------+
 * | Offsets ([n + 1] * EID bytes)      |
 * +------------------------------------+
 * | Neighborhood lists (m * VID bytes) |
 * +------------------------------------+
 * | Vertex weights (n * VWGT bytes)    |
 * +------------------------------------+
 * | Edge weights (m * EWGT bytes)      |
 * +------------------------------------+
 *
 * Header:
 * - 8 bytes: VERSION
 * - 8 bytes: n (number of vertices)
 * - 8 bytes: m (number of directed edges = 2x the number of undirected edges)
 *
 * VERSION is a bitfield that encodes the following information:
 * LSB -> [!ewgt|!vwgt|32bit_eid|32bit_vid|32bit_vwgt|32bit_ewgt] -> MSB
 *
 * - !ewgt: 1 = no edge weights, 0 = edge weights
 * - !vwgt: 1 = no vertex weights, 0 = vertex weights
 * - 32bit_eid: 1 = EID is 32 bit, 0 = EID is 64 bit
 * - 32bit_vid: 1 = VID is 32 bit, 0 = VID is 64 bit
 * - 32bit_vwgt: 1 = VWGT is 32 bit, 0 = VWGT is 64 bit
 * - 32bit_ewgt: 1 = EWGT is 32 bit, 0 = EWGT is 64 bit
 *
 * The entries of "Offsets" should be addresses relative to the start of the
 * file, such that Offsets[v] points to the first neighbor of vertex v in
 * "Neighborhood lists".
 * Note that this offset changes when using 32 vs 64 bit edge- or vertex IDs.
 *
 * ParHiP only supports graphs with 64 bit IDs and weights, no edge weights and
 * no vertex weights (i.e., VERSION must be 3).
 *
 * @file   parhip.h
 * @author Daniel Seemaier
 ******************************************************************************/
#pragma once

#include "kagen/context.h"
#include "kagen/io/graph_format.h"

#include <fstream>
#include <string>

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
    std::vector<std::string> DefaultExtensions() const final {
        return {"parhip", "bgf"}; // Keep *.bgf for legacy reasons
    }

    std::unique_ptr<GraphReader> CreateReader(const InputGraphConfig& config, PEID rank, PEID size) const final;

    std::unique_ptr<GraphWriter>
    CreateWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size) const final;
};
} // namespace kagen
