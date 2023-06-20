#include "kagen/io/coordinates.h"

#include "kagen/io/buffered_writer.h"
#include "kagen/io/graph_format.h"

namespace kagen {
CoordinatesWriter::CoordinatesWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size)
    : StandardGraphWriter(config, graph, info, rank, size) {}

void CoordinatesWriter::WriteHeader(const std::string&, SInt, SInt) {}

bool CoordinatesWriter::WriteBody(const std::string& filename) {
    RequiresCoordinates();

    BufferedTextOutput<> out(tag::append, filename);

    for (const auto& [x, y]: graph_.coordinates.first) {
        out.WriteFloat(x).WriteChar(' ').WriteFloat(y).WriteChar(' ').WriteFloat(0.0).WriteChar('\n').Flush();
    }
    for (const auto& [x, y, z]: graph_.coordinates.second) {
        out.WriteFloat(x).WriteChar(' ').WriteFloat(y).WriteChar(' ').WriteFloat(z).WriteChar('\n').Flush();
    }

    return false;
}

std::unique_ptr<GraphWriter> CoordinatesFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<CoordinatesWriter>(config, graph, info, rank, size);
}
} // namespace kagen
