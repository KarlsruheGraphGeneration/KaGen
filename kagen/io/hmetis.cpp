#include "kagen/io/hmetis.h"

#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/io/buffered_writer.h"
#include "kagen/io/hmetis_hypergraph.h"

#include <unordered_map>

namespace kagen {
HmetisWriter::HmetisWriter(
    const bool directed, const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank,
    const PEID size)
    : StandardGraphWriter(config, graph, info, rank, size),
      directed_(directed) {}

void HmetisWriter::WriteHeader(const std::string& filename, const SInt n, const SInt m) {
    BufferedTextOutput<> out(tag::append, filename);
    if (graph_.IsHypergraph()) {
        // TODO(clickup)[2026-05-10]: Weights for Hypergraphs not supported
        IgnoresEdgeWeights();
        IgnoresVertexWeights();

        const SInt num_hyperedges =
            graph_.hyperedge_offsets.empty() ? 0 : static_cast<SInt>(graph_.hyperedge_offsets.size() - 1);

        out.WriteInt(num_hyperedges).WriteChar(' ').WriteInt(n).WriteChar('\n').Flush();
        return;
    }

    out.WriteInt(m / 2).WriteChar(' ').WriteInt(n);
    if (info_.has_vertex_weights || info_.has_edge_weights) {
        out.WriteChar(' ');
        if (info_.has_vertex_weights) {
            out.WriteChar('1');
        }
        if (info_.has_edge_weights) {
            out.WriteChar('1');
        } else {
            out.WriteChar('0');
        }
    }
    out.WriteChar('\n').Flush();
}

bool HmetisWriter::WriteBody(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);

    if (graph_.IsHypergraph()) {
        IgnoresEdgeWeights();

        if (graph_.hyperedge_offsets.empty()) {
            return info_.has_vertex_weights;
        }

        for (SInt e = 0; e + 1 < static_cast<SInt>(graph_.hyperedge_offsets.size()); ++e) {
            const SInt begin = graph_.hyperedge_offsets[e];
            const SInt end   = graph_.hyperedge_offsets[e + 1];

            for (SInt i = begin; i < end; ++i) {
                if (i > begin) {
                    out.WriteChar(' ');
                }

                // hMetis uses 1-based vertex IDs
                out.WriteInt(graph_.hyperedge_pins[i] + 1);
            }

            out.WriteChar('\n').Flush();
        }

        return info_.has_vertex_weights;
    }

    for (SInt e = 0; e < graph_.edges.size(); ++e) {
        const auto& [from, to] = graph_.edges[e];

        if (info_.has_edge_weights) {
            out.WriteInt(static_cast<SInt>(graph_.edge_weights[e])).WriteChar(' ');
        }

        // Edges should only occur in one direction
        if (directed_ || from < to) {
            out.WriteInt(from + 1).WriteChar(' ').WriteInt(to + 1).WriteChar('\n').Flush();
        }
    }

    // Only call WriteFooter() if we have vertex weights
    return info_.has_vertex_weights;
}

void HmetisWriter::WriteFooter(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);
    for (const SSInt& weight: graph_.vertex_weights) {
        out.WriteInt(weight).WriteChar('\n').Flush();
    }
}

std::unique_ptr<GraphWriter> HmetisFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    if (IsHypergraph(graph)) {
        return std::make_unique<HmetisHypergraphWriter>(config, graph, info, rank, size);
    }
    return std::make_unique<HmetisWriter>(false, config, graph, info, rank, size);
}

std::unique_ptr<GraphWriter> DirectedHmetisFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<HmetisWriter>(true, config, graph, info, rank, size);
}

HmetisEpWriter::HmetisEpWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size)
    : StandardGraphWriter(config, graph, info, rank, size) {
    if (size > 1) {
        throw IOError("HmetisEpWriter does not support distributed/chunked output");
    }
}

void HmetisEpWriter::WriteHeader(const std::string& filename, const SInt n, const SInt m) {
    IgnoresEdgeWeights();
    IgnoresVertexWeights();

    BufferedTextOutput<> out(tag::append, filename);
    out.WriteInt(n).WriteChar(' ').WriteInt(m / 2).WriteChar('\n').Flush();
}

bool HmetisEpWriter::WriteBody(const std::string& filename) {
    graph_.SortEdgelist();

    BufferedTextOutput<> out(tag::append, filename);

    // Node and edge weights are ignored
    std::unordered_map<SInt, SInt> hnodes;

    SInt next_hnode        = 0;
    auto map_edge_to_hnode = [&](const SInt u, const SInt v) {
        const SInt key    = std::min(u, v) * info_.global_n + std::max(u, v);
        const auto key_it = hnodes.find(key);
        if (key_it == hnodes.end()) {
            hnodes[key] = ++next_hnode;
            return next_hnode;
        } else {
            return key_it->second;
        }
    };

    SInt cur_edge = 0;
    for (SInt from = graph_.vertex_range.first; from < graph_.vertex_range.second; ++from) {
        bool has_more = cur_edge < graph_.edges.size() && std::get<0>(graph_.edges[cur_edge]) == from;
        while (has_more) {
            const SInt to    = std::get<1>(graph_.edges[cur_edge]);
            const SInt hnode = map_edge_to_hnode(from, to);
            out.WriteInt(hnode);

            ++cur_edge;
            has_more = cur_edge < graph_.edges.size() && std::get<0>(graph_.edges[cur_edge]) == from;
            if (has_more) {
                out.WriteChar(' ').Flush();
            }
        }
        out.WriteChar('\n').Flush();
    }

    return false;
}
std::unique_ptr<GraphWriter> HmetisEpFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<HmetisEpWriter>(config, graph, info, rank, size);
}
} // namespace kagen
