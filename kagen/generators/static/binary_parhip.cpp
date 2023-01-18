#include "kagen/generators/static/binary_parhip.h"

#include <array>
#include <fstream>

#include "kagen/definitions.h"

namespace kagen::staticgraph {
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
    return ~(version & 2);
}

bool HasEdgeWeights(const SInt version) {
    return ~(version & 1);
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

        const SInt offset = (3 + n_ + m_ + from) * sizeof(ParHipID);
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
} // namespace kagen::staticgraph
