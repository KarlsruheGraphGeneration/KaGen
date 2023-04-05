#include "kagen/io/binary_parhip.h"

#include <algorithm>
#include <array>
#include <fstream>

#include <mpi.h>

#include "kagen/io/buffered_writer.h"
#include "kagen/io/graph_writer.h"
#include "kagen/tools/statistics.h"

namespace kagen {
using ParHipID     = unsigned long long;
using ParHipWeight = SSInt;

BinaryParHipWriter::BinaryParHipWriter(Graph& graph, MPI_Comm comm) : GraphWriter(graph, comm) {}

std::string BinaryParHipWriter::DefaultExtension() const {
    return "bgf";
}

void BinaryParHipWriter::Write(const PGeneratorConfig& config) {
    PEID rank = 0, size = 0;
    MPI_Comm_rank(comm_, &rank);
    MPI_Comm_size(comm_, &size);

    const bool        output = !config.quiet && rank == ROOT;
    const std::string filename =
        config.output.extension ? config.output.filename + "." + DefaultExtension() : config.output.filename;

    if (config.output.distributed && output) {
        std::cout
            << "Warning: this file format does not support distributed output; writing the graph to a single file."
            << std::endl;
    }

    // Edges must be sorted in order to convert them to the CSR format
    if (!std::is_sorted(edges_.begin(), edges_.end())) {
        std::sort(edges_.begin(), edges_.end());
    }

    const ParHipID num_global_vertices = FindNumberOfGlobalNodes(vertex_range_, comm_);
    const ParHipID num_global_edges    = FindNumberOfGlobalEdges(edges_, comm_);

    // Header
    if (rank == ROOT && config.output.header != OutputHeader::NEVER) {
        // To be compatible with the original format used by ParHiP, we use the following versions:
        // 3 = no vertex or edge weights = compatible with ParHiP
        // 2 = no vertex weights, but edge weights
        // 1 = vertex weights, but no edge weights
        // 0 = vertex weights and edge weights
        // I.e., the negated "format" code used by the Metis format
        const ParHipID vertex_weights_bit = (static_cast<SInt>(HasVertexWeights()) ^ 1) << 1;
        const ParHipID edge_weights_bit   = static_cast<SInt>(HasEdgeWeights()) ^ 1;
        const ParHipID version            = vertex_weights_bit | edge_weights_bit;

        std::ofstream out(filename, std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
        out.write(reinterpret_cast<const char*>(&version), sizeof(ParHipID));
        out.write(reinterpret_cast<const char*>(&num_global_vertices), sizeof(ParHipID));
        out.write(reinterpret_cast<const char*>(&num_global_edges), sizeof(ParHipID));
    }

    // Generate the CSR data
    const SInt            num_local_vertices = vertex_range_.second - vertex_range_.first;
    std::vector<ParHipID> offset(num_local_vertices + 1);
    std::vector<ParHipID> edges(edges_.size());

    SInt cur_offset = (3 + num_global_vertices + 1) * sizeof(ParHipID); // 3 = header size
    SInt cur_edge   = 0;
    SInt cur_vertex = 0;
    for (SInt from = vertex_range_.first; from < vertex_range_.second; ++from) {
        offset[cur_vertex]         = cur_offset;
        const SInt cur_edge_before = cur_edge;

        while (cur_edge < edges_.size() && std::get<0>(edges_[cur_edge]) == from) {
            edges[cur_edge] = std::get<1>(edges_[cur_edge]);
            ++cur_edge;
        }

        const SInt degree = cur_edge - cur_edge_before;
        cur_offset += degree * sizeof(ParHipID);
        ++cur_vertex;
    }
    offset[cur_vertex] = cur_offset;

    // Write graph to binary file
    auto output_loop = [&](auto&& action) {
        for (PEID pe = 0; pe < size; ++pe) {
            if (pe == rank) {
                std::ofstream out(filename, std::ios_base::out | std::ios_base::binary | std::ios_base::app);
                action(out);
            }
            MPI_Barrier(comm_);
        }
    };

    // Write offset array
    output_loop([&offset](std::ofstream& out) {
        out.write(reinterpret_cast<const char*>(offset.data()), offset.size() * sizeof(ParHipID));
    });

    // Write edges array
    output_loop([&edges](std::ofstream& out) {
        out.write(reinterpret_cast<const char*>(edges.data()), edges.size() * sizeof(ParHipID));
    });

    // Write  weights
    static_assert(std::is_same_v<SSInt, ParHipWeight>); // @todo create copy if the data types do not match

    if (HasVertexWeights()) {
        output_loop([this](std::ofstream& out) {
            out.write(
                reinterpret_cast<const char*>(vertex_weights_.data()), vertex_weights_.size() * sizeof(ParHipWeight));
        });
    }
    if (HasEdgeWeights()) {
        output_loop([this](std::ofstream& out) {
            out.write(reinterpret_cast<const char*>(edge_weights_.data()), edge_weights_.size() * sizeof(ParHipWeight));
        });
    }
}

namespace {
using ParHipID     = unsigned long long;
using ParHipWeight = SSInt;

SInt OffsetToEdge(const SInt n, const SInt offset) {
    return (offset / sizeof(ParHipID)) - 3 - (n + 1);
}

SInt ReadFirstEdge(std::ifstream& in, const SInt n, const SInt u) {
    const SInt offset = (3 + u) * sizeof(ParHipID);
    in.seekg(static_cast<std::streamsize>(offset));

    ParHipID entry = 0;
    in.read(reinterpret_cast<char*>(&entry), sizeof(ParHipID));

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

BinaryParhipReader::BinaryParhipReader(const std::string& filename) : in_(filename) {
    std::array<ParHipID, 3> size;
    in_.read(reinterpret_cast<char*>(size.data()), 3 * sizeof(ParHipID));
    version_ = size[0];
    n_       = size[1];
    m_       = size[2];
}

GraphSize BinaryParhipReader::ReadSize() {
    return {n_, m_};
}

Graph BinaryParhipReader::Read(
    const SInt from, SInt to_node, const SInt to_edge, const GraphRepresentation representation) {
    if (to_node > n_) {
        to_node = FindNodeByEdge(to_edge);
    }

    // Read xadj array of the CSR representation
    const SInt num_local_nodes = to_node - from;
    in_.seekg((3 + from) * sizeof(ParHipID));
    auto xadj = ReadVector<ParHipID>(in_, num_local_nodes + 1);

    const SInt first_edge_offset = xadj.front();
    for (auto& entry: xadj) { // Transform file offsets to edge offsets
        entry = (entry - first_edge_offset) / sizeof(ParHipID);
    }

    // Read adjncy array
    const SInt num_local_edges = xadj.back();
    in_.seekg(first_edge_offset);
    auto adjncy = ReadVector<ParHipID>(in_, num_local_edges);

    Graph ans;
    ans.vertex_range   = {from, from + num_local_nodes};
    ans.representation = representation;

    // Build edge list array from xadj and adjncy
    if (representation == GraphRepresentation::EDGE_LIST) {
        ans.edges.reserve(num_local_edges);
        for (SInt u = 0; u < num_local_nodes; ++u) {
            const SInt first_edge         = xadj[u];
            const SInt first_invalid_edge = xadj[u + 1];
            for (SInt e = first_edge; e < first_invalid_edge; ++e) {
                const SInt v = adjncy[e];
                ans.edges.emplace_back(from + u, v);
            }
        }
    } else {
        ans.xadj   = std::move(xadj);
        ans.adjncy = std::move(adjncy);
    }

    // Load weights if the graph is weighted
    if (HasVertexWeights(version_)) {
        ans.vertex_weights.resize(num_local_nodes);

        const SInt offset = (3 + n_ + m_ + 1 + from) * sizeof(ParHipID);
        in_.seekg(offset);
        in_.read(reinterpret_cast<char*>(ans.vertex_weights.data()), num_local_nodes * sizeof(ParHipWeight));
    }
    if (HasEdgeWeights(version_)) {
        ans.edge_weights.resize(num_local_edges);

        SInt offset = first_edge_offset + m_ * sizeof(ParHipID);
        if (HasVertexWeights(version_)) {
            offset += n_ * sizeof(ParHipID);
        }
        in_.seekg(offset);
        in_.read(reinterpret_cast<char*>(ans.edge_weights.data()), num_local_edges * sizeof(ParHipWeight));
    }

    return ans;
}

SInt BinaryParhipReader::FindNodeByEdge(const SInt edge) {
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
} // namespace kagen
