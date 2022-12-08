#include "kagen/io/dot.h"

#include "kagen/io/buffered_writer.h"

namespace kagen {
DotWriter::DotWriter(Graph& graph, const bool directed_output, MPI_Comm comm)
    : SequentialGraphWriter(graph, comm),
      directed_output_(directed_output) {}

std::string DotWriter::DefaultExtension() const {
    return "dot";
}

int DotWriter::Requirements() const {
    return Requirement::NO_VERTEX_WEIGHTS | Requirement::NO_EDGE_WEIGHTS;
}

void DotWriter::AppendHeaderTo(const std::string& filename, SInt, SInt) {
    BufferedTextOutput<> out(tag::append, filename);
    const char*          type = directed_output_ ? "digraph" : "graph";
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

    const char* arrow = directed_output_ ? "->" : "--";

    for (const auto& [from, to]: edges_) {
        if (directed_output_ || from < to) { // need edges only once
            out.WriteInt(from + 1).WriteString(arrow).WriteInt(to + 1).WriteChar('\n').Flush();
        }
    }
}

void DotWriter::AppendFooterTo(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);
    out.WriteString("}\n").Flush();
}
} // namespace kagen
