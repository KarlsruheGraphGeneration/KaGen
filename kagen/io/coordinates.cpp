#include "kagen/io/coordinates.h"

#include "kagen/io/buffered_writer.h"
#include "kagen/io/graph_format.h"

namespace kagen {
CoordinatesWriter::CoordinatesWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm)
    : SequentialGraphWriter(config, graph, comm) {}

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

int CoordinatesWriter::Requirements() const {
    return SequentialGraphWriter::Requirement::COORDINATES;
}

std::unique_ptr<GraphWriter>
CoordinatesFactory::CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const {
    return std::make_unique<CoordinatesWriter>(config, graph, comm);
}
} // namespace kagen
