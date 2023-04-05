#include "kagen/io/dot.h"

#include "kagen/io/buffered_writer.h"
#include "kagen/io/seq_graph_writer.h"

namespace kagen {
DotWriter::DotWriter(const bool directed, const OutputGraphConfig& config, Graph& graph, MPI_Comm comm)
    : SequentialGraphWriter(config, graph, comm),
      directed_(directed) {}

int DotWriter::Requirements() const {
    return Requirement::NO_VERTEX_WEIGHTS | Requirement::NO_EDGE_WEIGHTS;
}

void DotWriter::AppendHeaderTo(const std::string& filename, SInt, SInt) {
    BufferedTextOutput<> out(tag::append, filename);
    const char*          type = directed_ ? "digraph" : "graph";
    out.WriteString(type).WriteString(" G{\n").Flush();
}

void DotWriter::AppendTo(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);

    if (!coordinates_.first.empty()) {
        auto& coordinates = coordinates_.first; // 2D
        for (SInt node = vertex_range_.first; node < vertex_range_.second; ++node) {
            const auto& [x, y] = coordinates[node - vertex_range_.first];
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

    for (const auto& [from, to]: edges_) {
        if (directed_ || from < to) { // need edges only once
            out.WriteInt(from + 1).WriteString(arrow).WriteInt(to + 1).WriteChar('\n').Flush();
        }
    }
}

void DotWriter::AppendFooterTo(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);
    out.WriteString("}\n").Flush();
}

std::unique_ptr<GraphReader> DotFactory::CreateReader(const InputGraphConfig&) const {
    return nullptr;
}

std::unique_ptr<GraphWriter>
DotFactory::CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const {
    return std::make_unique<DotWriter>(false, config, graph, comm);
}

std::unique_ptr<GraphReader> DirectedDotFactory::CreateReader(const InputGraphConfig&) const {
    return nullptr;
}

std::unique_ptr<GraphWriter>
DirectedDotFactory::CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const {
    return std::make_unique<DotWriter>(true, config, graph, comm);
}
} // namespace kagen
