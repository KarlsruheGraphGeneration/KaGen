#include "kagen/io/graph_format.h"

#include "kagen/definitions.h"
#include "kagen/io.h"

#include <mpi.h>

namespace kagen {
GraphInfo::GraphInfo(const Graph& graph, MPI_Comm comm)
    : local_n(graph.NumberOfLocalVertices()),
      local_m(graph.NumberOfLocalEdges()),
      global_n(graph.NumberOfLocalVertices()),
      global_m(graph.NumberOfLocalEdges()),
      has_vertex_weights(!graph.vertex_weights.empty()),
      has_edge_weights(!graph.edge_weights.empty()) {
    MPI_Allreduce(MPI_IN_PLACE, &global_n, 1, KAGEN_MPI_SINT, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &global_m, 1, KAGEN_MPI_SINT, MPI_SUM, comm);
    MPI_Exscan(&local_n, &offset_n, 1, KAGEN_MPI_SINT, MPI_SUM, comm);
    MPI_Exscan(&local_m, &offset_m, 1, KAGEN_MPI_SINT, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &has_vertex_weights, 1, MPI_C_BOOL, MPI_LOR, comm);
    MPI_Allreduce(MPI_IN_PLACE, &has_edge_weights, 1, MPI_C_BOOL, MPI_LOR, comm);
}

GraphWriter::GraphWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size)
    : config_(config),
      info_(info),
      graph_(graph),
      rank_(rank),
      size_(size) {}

void GraphWriter::RequiresCoordinates() const {
    const SInt local_n = graph_.vertex_range.second - graph_.vertex_range.first;
    if (graph_.coordinates.first.size() != local_n && graph_.coordinates.second.size() != local_n) {
        throw IOError("output format requires coordinates, but the graph was generated without coordinates");
    }
}

void GraphWriter::Requires2DCoordinates() const {
    if (graph_.coordinates.first.size() != graph_.vertex_range.second - graph_.vertex_range.first) {
        throw IOError("output format requires 2D coordinates, but the graph was generated without 2D coordinates");
    }
}

void GraphWriter::Requires3DCoordinates() const {
    if (graph_.coordinates.second.size() != graph_.vertex_range.second - graph_.vertex_range.first) {
        throw IOError("output format requires 3D coordinates, but the graph was generated without 3D coordinates");
    }
}

void GraphWriter::IgnoresVertexWeights() const {
    if (rank_ == ROOT && info_.has_vertex_weights) {
        std::cerr
            << "Warning: output format does not support vertex weights, but the graph was generated with vertex weights"
            << std::endl;
    }
}

void GraphWriter::IgnoresEdgeWeights() const {
    if (rank_ == ROOT && info_.has_edge_weights) {
        std::cerr
            << "Warning: output format does not support edge weights, but the graph was generated with edge weights"
            << std::endl;
    }
}

StandardGraphWriter::StandardGraphWriter(
    const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size)
    : GraphWriter(config, graph, info, rank, size) {}

bool StandardGraphWriter::Write(const int pass, const std::string& filename) {
    const bool write_header_footer = [&] {
        if (config_.distributed) {
            return (rank_ == ROOT && config_.header == OutputHeader::ROOT) || config_.header == OutputHeader::ALWAYS;
        } else {
            return rank_ == ROOT;
        }
    }();

    if (pass == 0) {
        if (write_header_footer) {
            const SInt n =
                config_.distributed ? graph_.vertex_range.second - graph_.vertex_range.first : info_.global_n;
            const SInt m = config_.distributed ? graph_.edges.size() : info_.global_m;
            WriteHeader(filename, n, m);
        }
        return WriteBody(filename);
    }
    if (write_header_footer) {
        WriteFooter(filename);
    }
    return false;
}

void StandardGraphWriter::WriteFooter(const std::string&) {}
} // namespace kagen
