#include "kagen/generators/geometric/rgg/rgg_2d.h"

#include "kagen/generators/geometric/geometric_2d.h"
#include "kagen/tools/geometry.h"

namespace kagen {
RGG2D::RGG2D(const PGeneratorConfig& config, const PEID rank, const PEID size) : SpatialGrid2D(config, rank, size) {}

void RGG2D::GenerateEdges(const SInt chunk_row, const SInt chunk_column) {
    // Iterate grid cells
    for (SInt cell_row = 0; cell_row < cells_per_dim_; ++cell_row) {
        for (SInt cell_column = 0; cell_column < cells_per_dim_; ++cell_column) {
            // Iterate neighboring cells
            for (SSInt i = -1; i <= 1; i++) {
                SSInt neighbor_row = cell_row + i;
                for (SSInt j = -1; j <= 1; j++) {
                    SSInt neighbor_column = cell_column + j;

                    // Compute diffs
                    int horizontal_diff = 0;
                    int vertical_diff   = 0;
                    if (neighbor_column < 0)
                        horizontal_diff = -1;
                    else if (neighbor_column >= (SSInt)cells_per_dim_)
                        horizontal_diff = 1;
                    if (neighbor_row < 0)
                        vertical_diff = -1;
                    else if (neighbor_row >= (SSInt)cells_per_dim_)
                        vertical_diff = 1;

                    // Get correct grid cells
                    SInt neighbor_cell_row = (neighbor_row % (SSInt)cells_per_dim_ + cells_per_dim_) % cells_per_dim_;
                    SInt neighbor_cell_column =
                        (neighbor_column % (SSInt)cells_per_dim_ + cells_per_dim_) % cells_per_dim_;

                    // Skip invalid cells
                    if ((SSInt)chunk_row + vertical_diff < 0 || chunk_row + vertical_diff >= chunks_per_dim_
                        || (SSInt)chunk_column + horizontal_diff < 0
                        || chunk_column + horizontal_diff >= chunks_per_dim_)
                        continue;

                    // Get grid buckets for each cell
                    SInt chunk_id         = Encode(chunk_column, chunk_row);
                    SInt cell_id          = cell_row * cells_per_dim_ + cell_column;
                    SInt neighbor_id      = Encode(chunk_column + horizontal_diff, chunk_row + vertical_diff);
                    SInt neighbor_cell_id = neighbor_cell_row * cells_per_dim_ + neighbor_cell_column;

                    // If neighbor is local chunk skip
                    if (chunk_id > neighbor_id && IsLocalChunk(neighbor_id))
                        continue;
                    // Skip grid cells with lower id
                    if (chunk_id == neighbor_id && cell_id > neighbor_cell_id)
                        continue;

                    GenerateGridEdges(chunk_id, cell_id, neighbor_id, neighbor_cell_id);
                }
            }
        }
    }
}

void RGG2D::GenerateGridEdges(
    const SInt first_chunk_id, const SInt first_cell_id, const SInt second_chunk_id, const SInt second_cell_id) {
    // Check if vertices not generated
    // SInt first_global_cell_id = ComputeGlobalCellId(first_chunk_id, first_cell_id);
    // SInt second_global_cell_id =
    //    ComputeGlobalCellId(second_chunk_id, second_cell_id);

    // Gather vertices (temp)
    std::vector<Vertex> vertices_first;
    std::vector<Vertex> vertices_second;
    GenerateVertices(first_chunk_id, first_cell_id, vertices_first);
    GenerateVertices(second_chunk_id, second_cell_id, vertices_second);
    if (vertices_first.size() == 0)
        return;
    if (vertices_second.size() == 0)
        return;

    // Gather vertices
    // GenerateVertices(first_chunk_id, first_cell_id);
    // GenerateVertices(second_chunk_id, second_cell_id);
    // if (vertices_.find(first_global_cell_id) == end(vertices_)) return;
    // if (vertices_.find(second_global_cell_id) == end(vertices_)) return;
    // const std::vector<Vertex> &vertices_first = vertices_[first_global_cell_id];
    // const std::vector<Vertex> &vertices_second = vertices_[second_global_cell_id];
    // Generate edges
    // Same cell

    if (first_chunk_id == second_chunk_id && first_cell_id == second_cell_id) {
        for (SInt i = 0; i < vertices_first.size(); ++i) {
            const Vertex& v1 = vertices_first[i];
            for (SInt j = i + 1; j < vertices_second.size(); ++j) {
                const Vertex& v2               = vertices_second[j];
                const auto    squared_distance = PGGeometry<LPFloat>::SquaredEuclideanDistance(v1, v2);
                if (squared_distance <= target_r_) {
                    PushEdge(std::get<2>(v1), std::get<2>(v2));
                    PushWeightIfRequested(config_.edge_weights, squared_distance, target_r_);
                    PushEdge(std::get<2>(v2), std::get<2>(v1));
                    PushWeightIfRequested(config_.edge_weights, squared_distance, target_r_);
                }
            }
        }
    } else {
        for (SInt i = 0; i < vertices_first.size(); ++i) {
            const Vertex& v1 = vertices_first[i];
            for (SInt j = 0; j < vertices_second.size(); ++j) {
                const Vertex& v2               = vertices_second[j];
                const auto    squared_distance = PGGeometry<LPFloat>::SquaredEuclideanDistance(v1, v2);
                if (squared_distance <= target_r_) {
                    PushEdge(std::get<2>(v1), std::get<2>(v2));
                    PushWeightIfRequested(config_.edge_weights, squared_distance, target_r_);
                    if (IsLocalChunk(second_chunk_id)) {
                        PushEdge(std::get<2>(v2), std::get<2>(v1));
                        PushWeightIfRequested(config_.edge_weights, squared_distance, target_r_);
                    }
                }
            }
        }
    }
}

void RGG2D::GenerateEdgeList() {
    GenerateGeometry();
}

} // namespace kagen
