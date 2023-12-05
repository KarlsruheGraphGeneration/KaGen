#include "kagen/io/freight-netl.h"

#include "kagen/io/buffered_writer.h"

namespace kagen {
FreightNetlEpWriter::FreightNetlEpWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size)
    : StandardGraphWriter(config, graph, info, rank, size) {}

void FreightNetlEpWriter::WriteHeader(const std::string& filename, const SInt n, const SInt m) {
    BufferedTextOutput<> out(tag::append, filename);
    out.WriteInt(m / 2).WriteChar(' ').WriteInt(n).WriteChar(' ').WriteString("11").WriteChar('\n').Flush();
}

bool FreightNetlEpWriter::WriteBody(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);

    // Ignores edge and node weights

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
