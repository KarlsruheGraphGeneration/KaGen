#include "kagen/io/parhip.h"

#include "kagen/definitions.h"
#include "kagen/io.h"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <fstream>

namespace kagen {
using namespace parhip;

constexpr static bool kDebug = false;

namespace parhip {
ParhipID BuildVersion(
    const bool has_vertex_weights, const bool has_edge_weights, const bool has_32bit_edge_ids,
    const bool has_32bit_vertex_ids, const bool has_32bit_vertex_weights, const bool has_32bit_edge_weights) {
    // To be compatible with the original format used by ParHiP, we use the following versions:
    // 3 = no vertex or edge weights = compatible with ParHiP
    // 2 = no vertex weights, but edge weights
    // 1 = vertex weights, but no edge weights
    // 0 = vertex weights and edge weights
    // I.e., the negated "format" code used by the Metis format
    // Higher bits are used to encode the data types used to store the graph; if set to 0, we use 64 bits
    const ParhipID edge_weights_bit        = static_cast<SInt>(has_edge_weights) ^ 1;
    const ParhipID vertex_weights_bit      = (static_cast<SInt>(has_vertex_weights) ^ 1) << 1;
    const ParhipID edge_id_width_bit       = static_cast<SInt>(has_32bit_edge_ids) << 2;
    const ParhipID vertex_id_width_bit     = static_cast<SInt>(has_32bit_vertex_ids) << 3;
    const ParhipID vertex_weight_width_bit = static_cast<SInt>(has_32bit_vertex_weights) << 3;
    const ParhipID edge_weight_width_bit   = static_cast<SInt>(has_32bit_edge_weights) << 4;

    return vertex_weights_bit | edge_weights_bit | edge_id_width_bit | vertex_id_width_bit | vertex_weight_width_bit
           | edge_weight_width_bit;
}

bool HasVertexWeights(const SInt version) {
    return (version & 2) == 0;
}

bool HasEdgeWeights(const SInt version) {
    return (version & 1) == 0;
}

bool Has32BitEdgeIDs(const SInt version) {
    return version & 4;
}

bool Has32BitVertexIDs(const SInt version) {
    return version & 8;
}

bool Has32BitVertexWeights(const SInt version) {
    return version & 16;
}

bool Has32BitEdgeWeights(const SInt version) {
    return version & 32;
}
} // namespace parhip

ParhipWriter::ParhipWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size)
    : GraphWriter(config, graph, info, rank, size) {}

void ParhipWriter::WriteHeader(const std::string& filename) {
    // Edges must be sorted in order to convert them to the CSR format
    if (!std::is_sorted(graph_.edges.begin(), graph_.edges.end())) {
        std::sort(graph_.edges.begin(), graph_.edges.end());
    }

    // Header
    if (rank_ == ROOT && config_.header != OutputHeader::NEVER) {
        std::ofstream  out(filename, std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
        const ParhipID version =
            BuildVersion(info_.has_vertex_weights, info_.has_edge_weights, false, config_.width == 32, false, false);
        const ParhipID global_n = info_.global_n;
        const ParhipID global_m = info_.global_m;
        out.write(reinterpret_cast<const char*>(&version), sizeof(ParhipID));
        out.write(reinterpret_cast<const char*>(&global_n), sizeof(ParhipID));
        out.write(reinterpret_cast<const char*>(&global_m), sizeof(ParhipID));
    }
}

void ParhipWriter::WriteOffsets(const std::string& filename) {
    std::vector<ParhipID> offset(graph_.NumberOfLocalVertices() + 1);

    const int vertex_id_width = config_.width == 32 ? 4 : 8;
    SInt cur_offset = 3 * sizeof(ParhipID) + (info_.global_n + 1) * sizeof(ParhipID) + info_.offset_m * vertex_id_width;
    SInt cur_edge   = 0;
    SInt cur_vertex = 0;

    for (SInt from = graph_.vertex_range.first; from < graph_.vertex_range.second; ++from) {
        offset[cur_vertex] = cur_offset;

        const SInt cur_edge_before = cur_edge;
        while (cur_edge < graph_.edges.size() && std::get<0>(graph_.edges[cur_edge]) == from) {
            ++cur_edge;
        }
        const SInt degree = cur_edge - cur_edge_before;

        cur_offset += degree * vertex_id_width;
        ++cur_vertex;
    }
    offset[cur_vertex] = cur_offset;

    std::ofstream out(filename, std::ios_base::out | std::ios_base::binary | std::ios_base::app);

    // Case distinction: the last offset acts as a guardian and should only be written by the last PE
    if (rank_ + 1 == size_) {
        out.write(reinterpret_cast<const char*>(offset.data()), offset.size() * sizeof(ParhipID));
    } else {
        out.write(reinterpret_cast<const char*>(offset.data()), (offset.size() - 1) * sizeof(ParhipID));
    }
}

namespace {
template <typename T, SInt buf_size = 1024 * 1024>
void WriteEdgesImpl(std::ofstream& out, const Edgelist& edges) {
    std::vector<T> buf(buf_size);

    for (SInt pos = 0; pos < edges.size(); pos += buf_size) {
        const SInt count = std::min(edges.size() - pos, buf_size);
        std::transform(edges.begin() + pos, edges.begin() + pos + count, buf.begin(), [&](const auto& edge) {
            return static_cast<T>(std::get<1>(edge));
        });
        out.write(reinterpret_cast<const char*>(buf.data()), count * sizeof(T));
    }
}
} // namespace

void ParhipWriter::WriteEdges(const std::string& filename) {
    std::ofstream out(filename, std::ios_base::out | std::ios_base::binary | std::ios_base::app);
    if (config_.width == 32) {
        WriteEdgesImpl<std::uint32_t>(out, graph_.edges);
    } else {
        WriteEdgesImpl<ParhipID>(out, graph_.edges);
    }
}

void ParhipWriter::WriteVertexWeights(const std::string& filename) {
    std::ofstream out(filename, std::ios_base::out | std::ios_base::binary | std::ios_base::app);
    out.write(
        reinterpret_cast<const char*>(graph_.vertex_weights.data()),
        graph_.vertex_weights.size() * sizeof(ParhipWeight));
}

void ParhipWriter::WriteEdgeWeights(const std::string& filename) {
    std::ofstream out(filename, std::ios_base::out | std::ios_base::binary | std::ios_base::app);
    out.write(
        reinterpret_cast<const char*>(graph_.edge_weights.data()), graph_.edge_weights.size() * sizeof(ParhipWeight));
}

bool ParhipWriter::Write(const int pass, const std::string& filename) {
    if (config_.distributed) {
        throw IOError("ParHiP format does not support distributed output");
    }

    graph_.SortEdgelist();

    // @todo create copies if data types mismatch
    static_assert(std::is_same_v<SInt, ParhipID>);
    static_assert(std::is_same_v<SSInt, ParhipWeight>);

    switch (pass) {
        case 0:
            WriteHeader(filename);
            WriteOffsets(filename);
            return true;

        case 1:
            WriteEdges(filename);
            return info_.has_vertex_weights || info_.has_edge_weights;

        case 2:
            if (info_.has_vertex_weights) {
                WriteVertexWeights(filename);
            } else {
                WriteEdgeWeights(filename);
            }
            return info_.has_vertex_weights && info_.has_edge_weights;

        case 3:
            WriteEdgeWeights(filename);
            return false;
    }

    return false;
}

namespace {
SInt OffsetToEdge(const SInt version, const SInt n, const SInt offset) {
    const int edge_id_width   = Has32BitEdgeIDs(version) ? 4 : 8;
    const int vertex_id_width = Has32BitVertexIDs(version) ? 4 : 8;
    return (offset - 3 * sizeof(ParhipID) - (n + 1) * edge_id_width) / vertex_id_width;
}

SInt ReadFirstEdge(std::ifstream& in, const SInt version, const SInt n, const SInt u) {
    const int  edge_id_width = Has32BitVertexIDs(version) ? 4 : 8;
    const SInt offset        = 3 * sizeof(ParhipID) + u * edge_id_width;
    in.seekg(static_cast<std::streamsize>(offset));
    if (in.rdstate()) {
        throw IOError("seeking to offset " + std::to_string(offset) + " failed");
    }

    ParhipID entry = 0;
    in.read(reinterpret_cast<char*>(&entry), sizeof(ParhipID));
    if (in.rdstate()) {
        throw IOError(
            "reading " + std::to_string(sizeof(ParhipID)) + " bytes at offset " + std::to_string(offset) + " failed");
    }

    return OffsetToEdge(version, n, entry);
}

template <typename T, typename F, SInt buf_size = 1024 * 1024>
std::vector<T> ReadVector(std::ifstream& in, const SInt length) {
    std::vector<T> ans(length);
    if constexpr (std::is_same_v<T, F>) {
        const SInt size = length * sizeof(T);
        in.read(reinterpret_cast<char*>(ans.data()), size);
        if (in.rdstate()) {
            throw IOError(
                "reading " + std::to_string(size) + " bytes failed: only read " + std::to_string(in.gcount())
                + " bytes");
        }
    } else {
        std::vector<F> buf(buf_size);
        for (SInt pos = 0; pos < length; pos += buf_size) {
            const SInt count = std::min(length - pos, buf_size);
            const SInt size  = count * sizeof(F);
            const SInt at    = in.tellg();
            in.read(reinterpret_cast<char*>(buf.data()), size);
            if (in.rdstate()) {
                throw IOError(
                    "reading " + std::to_string(size) + " bytes failed: only read " + std::to_string(in.gcount())
                    + " bytes, at " + std::to_string(at));
            }
            std::copy(buf.begin(), buf.begin() + count, ans.begin() + pos);
        }
    }
    return ans;
}
} // namespace

ParhipReader::ParhipReader(const InputGraphConfig& config) : in_(config.filename) {
    if (!in_) {
        throw IOError("input file cannot be read");
    }

    std::array<ParhipID, 3> size;
    in_.read(reinterpret_cast<char*>(size.data()), 3 * sizeof(ParhipID));
    version_ = size[0];
    n_       = size[1];
    m_       = size[2];

    // Check file length
    in_.seekg(0, std::ios_base::end);
    const SInt actual_length = in_.tellg();

    SInt expected_length = 3 * sizeof(ParhipID);
    expected_length += (n_ + 1) * (Has32BitEdgeIDs(version_) ? 4 : 8);
    expected_length += m_ * (Has32BitVertexIDs(version_) ? 4 : 8);
    if (HasVertexWeights(version_)) {
        expected_length += n_ * (Has32BitVertexWeights(version_) ? 4 : 8);
    }
    if (HasEdgeWeights(version_)) {
        expected_length += m_ * (Has32BitEdgeWeights(version_) ? 4 : 8);
    }

    if constexpr (kDebug) {
        std::cout << "[Dbg] ParHiP file version: " << version_ << ", number of vertices: " << n_
                  << ", number of edges: " << m_ << "\n";
        std::cout << "[Dbg] Expected file length: " << expected_length << "\n";
        std::cout << "[Dbg] Actual file length: " << actual_length << "\n";
    }

    if (actual_length < expected_length) {
        throw IOError(
            "file is too short: expected " + std::to_string(expected_length) + " bytes, but file only contains "
            + std::to_string(actual_length) + " bytes");
    } else if (actual_length > expected_length) {
        std::cerr << "Warning: parhip input file is longer than expected: expected " << expected_length
                  << " bytes, but file contains " << actual_length << " bytes\n";
    }
}

std::pair<SInt, SInt> ParhipReader::ReadSize() {
    return {n_, m_};
}

Graph ParhipReader::Read(
    const SInt from_vertex, SInt to_vertex, SInt to_edge, const GraphRepresentation representation) {
    if (to_vertex > n_) {
        to_vertex = FindNodeByEdge(to_edge);
    }

    if constexpr (kDebug) {
        std::cout << "[Dbg] Reading from vertex " << from_vertex << " to vertex " << to_vertex << "; or to edge "
                  << to_edge << ", whichever is reached first\n";
    }

    const int edge_id_width       = Has32BitEdgeIDs(version_) ? 4 : 8;
    const int vertex_id_width     = Has32BitVertexIDs(version_) ? 4 : 8;
    const int vertex_weight_width = Has32BitVertexWeights(version_) ? 4 : 8;
    const int edge_weight_width   = Has32BitEdgeWeights(version_) ? 4 : 8;

    // Read xadj array of the CSR representation
    const SInt num_local_nodes = to_vertex - from_vertex;
    const SInt xadj_offset     = 3 * sizeof(ParhipID) + from_vertex * edge_id_width;
    in_.seekg(xadj_offset);
    if (in_.rdstate()) {
        throw IOError("seeking to offset " + std::to_string(xadj_offset) + " failed");
    }
    auto xadj = Has32BitEdgeIDs(version_) ? ReadVector<ParhipID, std::uint32_t>(in_, num_local_nodes + 1)
                                          : ReadVector<ParhipID, ParhipID>(in_, num_local_nodes + 1);

    const SInt first_global_edge = OffsetToEdge(version_, n_, xadj.front());
    if constexpr (kDebug) {
        std::cout << "[Dbg] First global edge to be read: " << first_global_edge
                  << ", computed from xadj.front()=" << xadj.front() << ", n=" << n_ << ", version_=" << version_
                  << std::endl;
    }

    const SInt offset_by = xadj.front();
    for (auto& entry: xadj) { // Transform file offsets to edge offsets
        entry = (entry - offset_by) / vertex_id_width;
    }

    // Read adjncy array
    const SInt num_local_edges = xadj.back() - xadj.front();
    const SInt adjncy_offset   = 3 * sizeof(ParhipID) + (n_ + 1) * edge_id_width + first_global_edge * vertex_id_width;
    in_.seekg(adjncy_offset);
    if (in_.rdstate()) {
        throw IOError("seeking to offset " + std::to_string(adjncy_offset) + " failed");
    }
    if constexpr (kDebug) {
        std::cout << "[Dbg] Reading " << num_local_edges << " edges from file position " << adjncy_offset << std::endl;
    }

    auto adjncy = Has32BitVertexIDs(version_) ? ReadVector<ParhipID, std::uint32_t>(in_, num_local_edges)
                                              : ReadVector<ParhipID, ParhipID>(in_, num_local_edges);

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
        const SInt vwgt_offset =
            3 * sizeof(ParhipID) + (n_ + 1) * edge_id_width + m_ * vertex_id_width + from_vertex * vertex_weight_width;
        in_.seekg(vwgt_offset);
        if (in_.rdstate()) {
            throw IOError("seeking to offset " + std::to_string(vwgt_offset) + " failed");
        }
        ans.vertex_weights = Has32BitVertexWeights(version_)
                                 ? ReadVector<ParhipWeight, std::int32_t>(in_, num_local_nodes)
                                 : ReadVector<ParhipWeight, ParhipWeight>(in_, num_local_nodes);
    }
    if (HasEdgeWeights(version_)) {
        SInt adjwgt_offset = 3 * sizeof(ParhipID) + (n_ + 1) * edge_id_width + m_ * vertex_id_width
                             + first_global_edge * edge_weight_width;
        if (HasVertexWeights(version_)) {
            adjwgt_offset += n_ * vertex_weight_width;
        }
        in_.seekg(adjwgt_offset);
        if (in_.rdstate()) {
            throw IOError("seeking to offset " + std::to_string(adjwgt_offset) + " failed");
        }
        ans.edge_weights = Has32BitEdgeWeights(version_) ? ReadVector<ParhipWeight, std::int32_t>(in_, num_local_edges)
                                                         : ReadVector<ParhipWeight, ParhipWeight>(in_, num_local_edges);
    }

    if (to_vertex >= n_ && to_edge >= m_) {
        const SInt end_pos = in_.tellg();
        in_.seekg(0, in_.end);
        const SInt real_end_pos = in_.tellg();
        if (end_pos != real_end_pos) {
            std::cerr << "Warning: " << (real_end_pos - end_pos) << " unexpected bytes at the end of the file"
                      << std::endl;
        }
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
        mid.second = ReadFirstEdge(in_, version_, n_, mid.first);

        if (mid.second < edge) {
            low = mid;
        } else {
            high = mid;
        }
    }

    return high.first;
}

std::unique_ptr<GraphReader> ParhipFactory::CreateReader(const InputGraphConfig& config, PEID, PEID) const {
    return std::make_unique<ParhipReader>(config);
}

std::unique_ptr<GraphWriter> ParhipFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<ParhipWriter>(config, graph, info, rank, size);
}
} // namespace kagen
