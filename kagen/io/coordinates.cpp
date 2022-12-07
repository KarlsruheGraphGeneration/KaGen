#include "kagen/io/coordinates.h"

#include "kagen/io/buffered_writer.h"
#include "kagen/io/graph_writer.h"

namespace kagen {
CoordinatesWriter::CoordinatesWriter(Graph& graph, MPI_Comm comm) : SequentialGraphWriter(graph, comm) {}

std::string CoordinatesWriter::DefaultExtension() const {
    return "xyz";
}

void CoordinatesWriter::AppendHeaderTo(const std::string&, SInt, SInt) {}

void CoordinatesWriter::AppendTo(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);

    for (const auto& [x, y]: coordinates_.first) {
        out.WriteFloat(x).WriteChar(' ').WriteFloat(y).WriteChar(' ').WriteFloat(0.0).WriteChar('\n').Flush();
    }
    for (const auto& [x, y, z]: coordinates_.second) {
        out.WriteFloat(x).WriteChar(' ').WriteFloat(y).WriteChar(' ').WriteFloat(z).WriteChar('\n').Flush();
    }
}

SequentialGraphWriter::Requirement CoordinatesWriter::Requirements() const {
    return SequentialGraphWriter::Requirement::COORDINATES;
}
} // namespace kagen
