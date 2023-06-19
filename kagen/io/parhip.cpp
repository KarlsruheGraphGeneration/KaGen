#include "kagen/io/parhip.h"

#include <algorithm>
#include <array>
#include <fstream>

#include <mpi.h>

#include "kagen/io/buffered_writer.h"
#include "kagen/tools/statistics.h"

namespace kagen {
using ParhipID     = unsigned long long;
using ParhipWeight = SSInt;

ParhipWriter::ParhipWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size)
    : GraphWriter(config, graph, info, rank, size) {}

std::string ParhipWriter::GetFilename() const {
    return config_.extension ? config_.filename + ".parhip" : config_.filename;
}

void ParhipWriter::WriteHeader() {
    // Edges must be sorted in order to convert them to the CSR format
    if (!std::is_sorted(graph_.edges.begin(), graph_.edges.end())) {
        std::sort(graph_.edges.begin(), graph_.edges.end());
    }

    // Header
    if (rank_ == ROOT && config_.header != OutputHeader::NEVER) {
        // To be compatible with the original format used by ParHiP, we use the following versions:
        // 3 = no vertex or edge weights = compatible with ParHiP
        // 2 = no vertex weights, but edge weights
        // 1 = vertex weights, but no edge weights
        // 0 = vertex weights and edge weights
        // I.e., the negated "format" code used by the Metis format
        const ParhipID vertex_weights_bit = (static_cast<SInt>(info_.has_vertex_weights) ^ 1) << 1;
        const ParhipID edge_weights_bit   = static_cast<SInt>(info_.has_edge_weights) ^ 1;
        const ParhipID version            = vertex_weights_bit | edge_weights_bit;

        std::ofstream  out(GetFilename(), std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
        const ParhipID global_n = info_.global_n;
        const ParhipID global_m = info_.global_m;
        out.write(reinterpret_cast<const char*>(&version), sizeof(ParhipID));
        out.write(reinterpret_cast<const char*>(&global_n), sizeof(ParhipID));
        out.write(reinterpret_cast<const char*>(&global_m), sizeof(ParhipID));
    }
}

void ParhipWriter::WriteOffsets() {
    std::vector<ParhipID> offset(graph_.NumberOfLocalVertices() + 1);

    SInt cur_offset = (3 + info_.global_n + 1) * sizeof(ParhipID); // 3 = header size
    SInt cur_edge   = 0;
    SInt cur_vertex = 0;

    for (SInt from = graph_.vertex_range.first; from < graph_.vertex_range.second; ++from) {
        offset[cur_vertex] = cur_offset;

        const SInt cur_edge_before = cur_edge;
        while (cur_edge < graph_.edges.size() && std::get<0>(graph_.edges[cur_edge]) == from) {
            ++cur_edge;
        }
        const SInt degree = cur_edge - cur_edge_before;

        cur_offset += degree * sizeof(ParhipID);
        ++cur_vertex;
    }
    offset[cur_vertex] = cur_offset;

    std::ofstream out(GetFilename(), std::ios_base::out | std::ios_base::binary | std::ios_base::app);
    out.write(reinterpret_cast<const char*>(offset.data()), offset.size() * sizeof(ParhipID));
}

void ParhipWriter::WriteEdges() {
    std::vector<ParhipID> edges(graph_.NumberOfLocalEdges());
    std::transform(graph_.edges.begin(), graph_.edges.end(), edges.begin(), [&](const auto& edge) {
        return static_cast<ParhipID>(std::get<1>(edge));
    });

    std::ofstream out(GetFilename(), std::ios_base::out | std::ios_base::binary | std::ios_base::app);
    out.write(reinterpret_cast<const char*>(edges.data()), edges.size() * sizeof(ParhipID));
}

void ParhipWriter::WriteVertexWeights() {
    std::ofstream out(GetFilename(), std::ios_base::out | std::ios_base::binary | std::ios_base::app);
    out.write(
        reinterpret_cast<const char*>(graph_.vertex_weights.data()),
        graph_.vertex_weights.size() * sizeof(ParhipWeight));
}

void ParhipWriter::WriteEdgeWeights() {
    std::ofstream out(GetFilename(), std::ios_base::out | std::ios_base::binary | std::ios_base::app);
    out.write(
        reinterpret_cast<const char*>(graph_.edge_weights.data()), graph_.edge_weights.size() * sizeof(ParhipWeight));
}

bool ParhipWriter::Write(const int pass) {
    if (config_.distributed) {
        throw IOError("ParHiP format does not support distributed output");
    }

    // @todo create copies if data types mismatch
    static_assert(std::is_same_v<SInt, ParhipID>);
    static_assert(std::is_same_v<SSInt, ParhipWeight>);

    switch (pass) {
        case 0:
            WriteHeader();
            WriteOffsets();
            return true;

        case 1:
            WriteEdges();
            return info_.has_vertex_weights || info_.has_edge_weights;

        case 2:
            if (info_.has_vertex_weights) {
                WriteVertexWeights();
            } else {
                WriteEdgeWeights();
            }
            return info_.has_vertex_weights && info_.has_edge_weights;

        case 3:
            WriteEdgeWeights();
            return false;
    }

    return false;
}

namespace {
SInt OffsetToEdge(const SInt n, const SInt offset) {
    return (offset / sizeof(ParhipID)) - 3 - (n + 1);
}

SInt ReadFirstEdge(std::ifstream& in, const SInt n, const SInt u) {
    const SInt offset = (3 + u) * sizeof(ParhipID);
    in.seekg(static_cast<std::streamsize>(offset));

    ParhipID entry = 0;
    in.read(reinterpret_cast<char*>(&entry), sizeof(ParhipID));

    return OffsetToEdge(n, entry);
}

bool HasVertexWeights(const SInt version) {
    return (version & 2) == 0;
}

bool HasEdgeWeights(const SInt version) {
    return (version & 1) == 0;
}

template <typename T>
std::vector<T> ReadVector(std::ifstream& in, const SInt length) {
    std::vector<T> ans(length);
    in.read(reinterpret_cast<char*>(ans.data()), length * sizeof(T));
    return ans;
}
} // namespace

ParhipReader::ParhipReader(const std::string& filename) : in_(filename) {
    if (!in_) {
        throw IOError("input file cannot be read");
    }

    std::array<ParhipID, 3> size;
    in_.read(reinterpret_cast<char*>(size.data()), 3 * sizeof(ParhipID));
    version_ = size[0];
    n_       = size[1];
    m_       = size[2];
}

std::pair<SInt, SInt> ParhipReader::ReadSize() {
    return {n_, m_};
}

Graph ParhipReader::Read(
    const SInt from_vertex, SInt to_vertex, SInt to_edge, const GraphRepresentation representation) {
    if (to_vertex > n_) {
        to_vertex = FindNodeByEdge(to_edge);
    }

    // Read xadj array of the CSR representation
    const SInt num_local_nodes = to_vertex - from_vertex;
    in_.seekg((3 + from_vertex) * sizeof(ParhipID));
    auto xadj = ReadVector<ParhipID>(in_, num_local_nodes + 1);

    const SInt first_edge_offset = xadj.front();
    for (auto& entry: xadj) { // Transform file offsets to edge offsets
        entry = (entry - first_edge_offset) / sizeof(ParhipID);
    }

    // Read adjncy array
    const SInt num_local_edges = xadj.back();
    in_.seekg(first_edge_offset);
    auto adjncy = ReadVector<ParhipID>(in_, num_local_edges);

    Graph ans;
    ans.vertex_range   = {from_vertex, from_vertex + num_local_nodes};
    ans.representation = representation;

    // Build edge list array from xadj and adjncy
    if (representation == GraphRepresentation::EDGE_LIST) {
        ans.edges.reserve(num_local_edges);
        for (SInt u = 0; u < num_local_nodes; ++u) {
            const SInt first_edge         = xadj[u];
            const SInt first_invalid_edge = xadj[u + 1];
            for (SInt e = first_edge; e < first_invalid_edge; ++e) {
                const SInt v = adjncy[e];
                ans.edges.emplace_back(from_vertex + u, v);
            }
        }
    } else {
        ans.xadj   = std::move(xadj);
        ans.adjncy = std::move(adjncy);
    }

    // Load weights if the graph is weighted
    if (HasVertexWeights(version_)) {
        ans.vertex_weights.resize(num_local_nodes);

        const SInt offset = (3 + n_ + m_ + 1 + from_vertex) * sizeof(ParhipID);
        in_.seekg(offset);
        in_.read(reinterpret_cast<char*>(ans.vertex_weights.data()), num_local_nodes * sizeof(ParhipWeight));
    }
    if (HasEdgeWeights(version_)) {
        ans.edge_weights.resize(num_local_edges);

        SInt offset = first_edge_offset + m_ * sizeof(ParhipID);
        if (HasVertexWeights(version_)) {
            offset += n_ * sizeof(ParhipID);
        }
        in_.seekg(offset);
        in_.read(reinterpret_cast<char*>(ans.edge_weights.data()), num_local_edges * sizeof(ParhipWeight));
    }

    return ans;
}

SInt ParhipReader::FindNodeByEdge(const SInt edge) {
    if (edge == 0) {
        return 0;
    }

    std::pair<SInt, SInt> low  = {0, 0};
    std::pair<SInt, SInt> high = {n_, m_ - 1};
    while (high.first - low.first > 1) {
        std::pair<SInt, SInt> mid;
        mid.first  = (low.first + high.first) / 2;
        mid.second = ReadFirstEdge(in_, n_, mid.first);

        if (mid.second < edge) {
            low = mid;
        } else {
            high = mid;
        }
    }

    return high.first;
}

std::unique_ptr<GraphReader> ParhipFactory::CreateReader(const InputGraphConfig& config, PEID, PEID) const {
    return std::make_unique<ParhipReader>(config.filename);
}

std::unique_ptr<GraphWriter> ParhipFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<ParhipWriter>(config, graph, info, rank, size);
}
} // namespace kagen
