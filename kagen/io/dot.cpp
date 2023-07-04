#include "kagen/io/dot.h"

#include "kagen/io/buffered_writer.h"

namespace kagen {
DotWriter::DotWriter(
    const bool directed, const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank,
    const PEID size)
    : StandardGraphWriter(config, graph, info, rank, size),
      directed_(directed) {}

void DotWriter::WriteHeader(const std::string& filename, SInt, SInt) {
    IgnoresVertexWeights();
    IgnoresEdgeWeights();

    BufferedTextOutput<> out(tag::append, filename);
    const char*          type = directed_ ? "digraph" : "graph";
    out.WriteString(type).WriteString(" G{\n").Flush();
}

bool DotWriter::WriteBody(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);

    if (!graph_.coordinates.first.empty()) {
        auto& coordinates = graph_.coordinates.first; // 2D
        for (SInt node = graph_.vertex_range.first; node < graph_.vertex_range.second; ++node) {
            const auto& [x, y] = coordinates[node - graph_.vertex_range.first];
            out.WriteInt(node + 1)
                .WriteString("[pos=\"")
                .WriteFloat(x * 10)
                .WriteChar(',')
                .WriteFloat(y * 10)
                .WriteString("!\"]\n")
                .Flush();
        }
    }

    const char* arrow = directed_ ? "->" : "--";
    for (const auto& [from, to]: graph_.edges) {
        if (directed_ || from < to) { // need edges only once
            out.WriteInt(from + 1).WriteString(arrow).WriteInt(to + 1).WriteChar('\n').Flush();
        }
    }

    return true;
}

void DotWriter::WriteFooter(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);
    out.WriteString("}\n").Flush();
}

std::unique_ptr<GraphWriter> DotFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<DotWriter>(false, config, graph, info, rank, size);
}

std::unique_ptr<GraphWriter> DirectedDotFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<DotWriter>(true, config, graph, info, rank, size);
}
} // namespace kagen
