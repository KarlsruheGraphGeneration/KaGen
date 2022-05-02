#include "kagen/io/io.h"

#include <algorithm>
#include <memory>

#include <mpi.h>

#include "kagen/definitions.h"
#include "kagen/io/buffered_writer.h"
#include "kagen/io/dot.h"
#include "kagen/io/edgelist.h"
#include "kagen/io/graph_writer.h"
#include "kagen/io/hmetis.h"
#include "kagen/io/metis.h"
#include "kagen/tools/statistics.h"

namespace kagen {
namespace {
std::unique_ptr<GraphWriter> CreateGraphWriter(
    const OutputFormat format, EdgeList& edges, const VertexRange vertex_range, Coordinates& coordinates,
    MPI_Comm comm) {
    switch (format) {
        case OutputFormat::NONE:
            return std::make_unique<NoopWriter>(edges, vertex_range, coordinates, comm);

        case OutputFormat::EDGE_LIST:
            return std::make_unique<EdgeListWriter>(edges, vertex_range, coordinates, comm);

        case OutputFormat::BINARY_EDGE_LIST:
            return std::make_unique<BinaryEdgeListWriter>(edges, vertex_range, coordinates, comm);

        case OutputFormat::METIS:
            return std::make_unique<MetisWriter>(edges, vertex_range, coordinates, comm);

        case OutputFormat::HMETIS:
            return std::make_unique<HMetisWriter>(edges, vertex_range, coordinates, comm);

        case OutputFormat::DOT:
            return std::make_unique<DotWriter>(edges, vertex_range, coordinates, comm);
    }

    __builtin_unreachable();
}
} // namespace

void WriteGraph(
    const PGeneratorConfig& config, EdgeList& edges, const VertexRange vertex_range, Coordinates& coordinates,
    MPI_Comm comm) {
    auto writer = CreateGraphWriter(config.output_format, edges, vertex_range, coordinates, comm);
    writer->Write(config);
}
} // namespace kagen
