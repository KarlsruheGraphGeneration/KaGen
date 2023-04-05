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
std::unique_ptr<GraphWriter>
CreateGraphWriter(const OutputFormat format, const int width, Graph& graph, MPI_Comm comm) {
    switch (format) {
        case OutputFormat::NONE:
            return std::make_unique<NoopWriter>(graph, comm);

        case OutputFormat::EDGE_LIST:
        case OutputFormat::EDGE_LIST_UNDIRECTED: {
            const bool undirected = format == OutputFormat::EDGE_LIST_UNDIRECTED;
            return std::make_unique<EdgeListWriter>(graph, comm, true, undirected);
        }

        case OutputFormat::BINARY_EDGE_LIST:
        case OutputFormat::BINARY_EDGE_LIST_UNDIRECTED:
        case OutputFormat::XTRAPULP: {
            const bool undirected = format != OutputFormat::BINARY_EDGE_LIST;
            const bool header     = format != OutputFormat::XTRAPULP;
            return std::make_unique<BinaryEdgeListWriter>(graph, comm, width, header, undirected);
        }

        case OutputFormat::METIS:
            return std::make_unique<MetisWriter>(graph, comm);

        case OutputFormat::HMETIS:
        case OutputFormat::HMETIS_DIRECTED: {
            const bool directed = format == OutputFormat::HMETIS_DIRECTED;
            return std::make_unique<HMetisWriter>(graph, comm, directed);
        }

        case OutputFormat::DOT:
        case OutputFormat::DOT_DIRECTED: {
            const bool directed = format == OutputFormat::DOT_DIRECTED;
            return std::make_unique<DotWriter>(graph, comm, directed);
        }

        case OutputFormat::COORDINATES:
            return std::make_unique<CoordinatesWriter>(graph, comm);

        case OutputFormat::PARHIP:
            return std::make_unique<BinaryParHipWriter>(graph, comm);
    }

    __builtin_unreachable();
}
} // namespace

void WriteGraph(const PGeneratorConfig& config, Graph& graph, MPI_Comm comm) {
    assert(graph.representation == GraphRepresentation::EDGE_LIST && "graph must be in edge list representation");

    for (const OutputFormat& format: config.output.formats) {
        auto writer = CreateGraphWriter(format, config.output.width, graph, comm);
        writer->Write(config);
    }
}
} // namespace kagen
