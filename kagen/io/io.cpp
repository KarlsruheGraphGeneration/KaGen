#include "kagen/io/io.h"

#include <algorithm>
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
std::unique_ptr<GraphWriter> CreateGraphWriter(const GraphFormat format, Graph& graph, MPI_Comm comm) {
    switch (format) {
        case GraphFormat::NONE:
            return std::make_unique<NoopWriter>(graph, comm);

        case GraphFormat::EDGE_LIST:
            return std::make_unique<EdgeListWriter>(graph, comm);

        case GraphFormat::BINARY_EDGE_LIST:
            return std::make_unique<BinaryEdgeListWriter>(graph, comm);

        case GraphFormat::METIS:
            return std::make_unique<MetisWriter>(graph, comm);

        case GraphFormat::HMETIS:
            return std::make_unique<HMetisWriter>(graph, comm);

        case GraphFormat::DOT:
            return std::make_unique<DotWriter>(graph, false, comm);

        case GraphFormat::DOT_DIRECTED:
            return std::make_unique<DotWriter>(graph, true, comm);

        case GraphFormat::COORDINATES:
            return std::make_unique<CoordinatesWriter>(graph, comm);
            
        case GraphFormat::BINARY_PARHIP:
            return std::make_unique<BinaryParHipWriter>(graph, comm);
    }

    __builtin_unreachable();
}
} // namespace

void WriteGraph(const PGeneratorConfig& config, Graph& graph, MPI_Comm comm) {
    auto writer = CreateGraphWriter(config.output_format, graph, comm);
    writer->Write(config);
}
} // namespace kagen
