#include "kagen/io/freight-netl.h"

#include "kagen/io/buffered_writer.h"

#include <unordered_map>

namespace kagen {
FreightNetlWriter::FreightNetlWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size)
    : StandardGraphWriter(config, graph, info, rank, size) {
    if (size > 1) {
        throw IOError("FreightNetlWriter does not support distributed/chunked output");
    }
}

void FreightNetlWriter::WriteHeader(const std::string& filename, const SInt n, const SInt m) {
    BufferedTextOutput<> out(tag::append, filename);
    out.WriteInt(n).WriteChar(' ').WriteInt(m / 2).WriteString(" 11\n").Flush();
}

bool FreightNetlWriter::WriteBody(const std::string& filename) {
    graph_.SortEdgelist();

    BufferedTextOutput<> out(tag::append, filename);

    std::unordered_map<SInt, SInt> nets;

    SInt next_net        = 0;
    auto map_edge_to_net = [&](const SInt u, const SInt v) {
        const SInt key    = std::min(u, v) * info_.global_n + std::max(u, v);
        const auto key_it = nets.find(key);
        if (key_it == nets.end()) {
            nets[key] = ++next_net;
            return next_net;
        } else {
            return key_it->second;
        }
    };

    SInt cur_edge = 0;
    for (SInt from = graph_.vertex_range.first; from < graph_.vertex_range.second; ++from) {
        if (info_.has_vertex_weights) {
            out.WriteInt(graph_.vertex_weights[from - graph_.vertex_range.first]);
        } else {
            out.WriteChar('1');
        }

        while (cur_edge < graph_.edges.size() && std::get<0>(graph_.edges[cur_edge]) == from) {
            const SInt to  = std::get<1>(graph_.edges[cur_edge]);
            const SInt net = map_edge_to_net(from, to);

            out.WriteChar(' ').WriteInt(net).WriteChar(' ');

            if (info_.has_edge_weights) {
                out.WriteInt(graph_.edge_weights[cur_edge]);
            } else {
                out.WriteChar('1');
            }
            out.Flush();

            ++cur_edge;
        }

        out.WriteChar('\n').Flush();
    }

    return false;
}

std::unique_ptr<GraphWriter> FreightNetlFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<FreightNetlWriter>(config, graph, info, rank, size);
}

FreightNetlEpWriter::FreightNetlEpWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size)
    : StandardGraphWriter(config, graph, info, rank, size) {}

void FreightNetlEpWriter::WriteHeader(const std::string& filename, const SInt n, const SInt m) {
    IgnoresEdgeWeights();
    IgnoresVertexWeights();

    BufferedTextOutput<> out(tag::append, filename);
    out.WriteInt(m / 2).WriteChar(' ').WriteInt(n).WriteString(" 11\n").Flush();
}

bool FreightNetlEpWriter::WriteBody(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);

    for (SInt e = 0; e < graph_.edges.size(); ++e) {
        const auto& [from, to] = graph_.edges[e];
        if (from >= to) {
            continue;
        }

        out.WriteChar('1')
            .WriteChar(' ')
            .WriteInt(from)
            .WriteChar(' ')
            .WriteChar('1')
            .WriteChar(' ')
            .WriteInt(to)
            .WriteChar(' ')
            .WriteChar('1')
            .WriteChar('\n')
            .Flush();
    }

    return false;
}

std::unique_ptr<GraphWriter> FreightNetlEpFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<FreightNetlEpWriter>(config, graph, info, rank, size);
}
} // namespace kagen
