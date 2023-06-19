#include "kagen/io/graph_format.h"

namespace kagen {
GraphWriter::GraphWriter(const OutputGraphConfig& config, Graph& graph, const PEID rank, const PEID size)
    : config_(config),
      graph_(graph),
      rank_(rank),
      size_(size) {}

bool GraphWriter::HasVertexWeights() const {
    bool has_vertex_weights =
        !vertex_weights_.empty() && vertex_weights_.size() == vertex_range_.second - vertex_range_.first;
    MPI_Allreduce(MPI_IN_PLACE, &has_vertex_weights, 1, MPI_CXX_BOOL, MPI_LAND, comm_);
    return has_vertex_weights;
}

bool GraphWriter::HasEdgeWeights() const {
    bool has_edge_weights = !edge_weights_.empty() && edge_weights_.size() == edges_.size();
    MPI_Allreduce(MPI_IN_PLACE, &has_edge_weights, 1, MPI_CXX_BOOL, MPI_LAND, comm_);
    return has_edge_weights;
}
} // namespace kagen
