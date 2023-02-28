#include "kagen/io/io.h"

#include <algorithm>
#include <cassert>
#include <memory>

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/io/binary_parhip.h"
#include "kagen/io/buffered_writer.h"
#include "kagen/io/coordinates.h"
#include "kagen/io/dot.h"
#include "kagen/io/edgelist.h"
#include "kagen/io/graph_writer.h"
#include "kagen/io/hmetis.h"
#include "kagen/io/metis.h"
#include "kagen/tools/statistics.h"

namespace kagen {
namespace {
std::unique_ptr<GraphWriter> CreateGraphWriter(const OutputFormat format, Graph& graph, MPI_Comm comm) {
    switch (format) {
        case OutputFormat::NONE:
            return std::make_unique<NoopWriter>(graph, comm);

        case OutputFormat::EDGE_LIST:
            return std::make_unique<EdgeListWriter>(graph, comm);

        case OutputFormat::BINARY_EDGE_LIST:
        case OutputFormat::XTRAPULP64:
            return std::make_unique<BinaryEdgeListWriter>(graph, comm, 64, format == OutputFormat::BINARY_EDGE_LIST);

        case OutputFormat::BINARY_EDGE_LIST32:
        case OutputFormat::XTRAPULP32:
            return std::make_unique<BinaryEdgeListWriter>(graph, comm, 32, format == OutputFormat::BINARY_EDGE_LIST32);

        case OutputFormat::METIS:
            return std::make_unique<MetisWriter>(graph, comm);

        case OutputFormat::HMETIS:
            return std::make_unique<HMetisWriter>(graph, comm);

        case OutputFormat::DOT:
            return std::make_unique<DotWriter>(graph, false, comm);

        case OutputFormat::DOT_DIRECTED:
            return std::make_unique<DotWriter>(graph, true, comm);

        case OutputFormat::COORDINATES:
            return std::make_unique<CoordinatesWriter>(graph, comm);

        case OutputFormat::BINARY_PARHIP:
            return std::make_unique<BinaryParHipWriter>(graph, comm);
    }

    __builtin_unreachable();
}
} // namespace

void WriteGraph(const PGeneratorConfig& config, Graph& graph, MPI_Comm comm) {
    assert(graph.representation == GraphRepresentation::EDGE_LIST && "graph must be in edge list representation");
    auto writer = CreateGraphWriter(config.output_format, graph, comm);
    writer->Write(config);
}
} // namespace kagen
