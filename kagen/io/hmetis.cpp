#include "kagen/io/hmetis.h"

#include "kagen/io/buffered_writer.h"

namespace kagen {
HmetisWriter::HmetisWriter(
    const bool directed, const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank,
    const PEID size)
    : StandardGraphWriter(config, graph, info, rank, size),
      directed_(directed) {}

void HmetisWriter::WriteHeader(const std::string& filename, const SInt n, const SInt m) {
    BufferedTextOutput<> out(tag::append, filename);
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
    return std::make_unique<HmetisWriter>(false, config, graph, info, rank, size);
}

std::unique_ptr<GraphWriter> DirectedHmetisFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<HmetisWriter>(true, config, graph, info, rank, size);
}
} // namespace kagen
