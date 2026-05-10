#include "kagen/generators/hyper/h_geometric/hyper_rgg2d.h"

#include "kagen/context.h"
#include "kagen/tools/geometry.h"

namespace kagen {

HyperRGG2D::HyperRGG2D(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : Geometric2D(config, rank, size),
      SpatialGrid2D(config, rank, size) {}

void HyperRGG2D::GenerateBallHyperedge(const Vertex& center, const SInt center_chunk_id, const SInt center_cell_id) {
    const SInt center_id = std::get<2>(center);

    std::vector<SInt> pins;
    pins.push_back(std::get<2>(center));

    SInt cell_x;
    SInt cell_y;
    DecodeCell(center_cell_id, cell_x, cell_y);

    SInt chunk_x;
    SInt chunk_y;
    Decode(center_chunk_id, chunk_x, chunk_y);

    for (SSInt dx = -1; dx <= 1; ++dx) {
        for (SSInt dy = -1; dy <= 1; ++dy) {
            const SSInt neighbor_cell_x = cell_x + dx;
            const SSInt neighbor_cell_y = cell_y + dy;

            int chunk_dx = 0;
            int chunk_dy = 0;

            if (neighbor_cell_x < 0) {
                chunk_dx = -1;
            } else if (neighbor_cell_x >= static_cast<SSInt>(cells_per_dim_)) {
                chunk_dx = 1;
            }

            if (neighbor_cell_y < 0) {
                chunk_dy = -1;
            } else if (neighbor_cell_y >= static_cast<SSInt>(cells_per_dim_)) {
                chunk_dy = 1;
            }

            const SSInt neighbor_chunk_x = static_cast<SSInt>(chunk_x) + chunk_dx;
            const SSInt neighbor_chunk_y = static_cast<SSInt>(chunk_y) + chunk_dy;

            if (neighbor_chunk_x < 0 || neighbor_chunk_x >= static_cast<SSInt>(chunks_per_dim_) || neighbor_chunk_y < 0
                || neighbor_chunk_y >= static_cast<SSInt>(chunks_per_dim_)) {
                continue;
            }

            const SInt wrapped_cell_x =
                (neighbor_cell_x % static_cast<SSInt>(cells_per_dim_) + cells_per_dim_) % cells_per_dim_;
            const SInt wrapped_cell_y =
                (neighbor_cell_y % static_cast<SSInt>(cells_per_dim_) + cells_per_dim_) % cells_per_dim_;

            const SInt neighbor_chunk_id = Encode(neighbor_chunk_x, neighbor_chunk_y);
            const SInt neighbor_cell_id  = EncodeCell(wrapped_cell_x, wrapped_cell_y);

            std::vector<Vertex> candidates;
            GenerateVertices(neighbor_chunk_id, neighbor_cell_id, candidates);

            for (const Vertex& candidate: candidates) {
                const SInt candidate_id = std::get<2>(candidate);

                if (candidate_id == center_id) {
                    continue;
                }

                const auto squared_distance = PGGeometry<LPFloat>::SquaredEuclideanDistance(center, candidate);

                if (squared_distance <= target_r_) {
                    pins.push_back(candidate_id);
                }
            }
        }
    }

    if (pins.size() > 1) {
        std::sort(pins.begin(), pins.end());
        pins.erase(std::unique(pins.begin(), pins.end()), pins.end());

        PushHyperedge(pins);
    }
}

void HyperRGG2D::GenerateEdges(SInt chunk_row, SInt chunk_column) {
    const SInt chunk_id = Encode(chunk_column, chunk_row);

    if (!IsLocalChunk(chunk_id)) {
        return;
    }

    for (SInt cell_row = 0; cell_row < cells_per_dim_; ++cell_row) {
        for (SInt cell_column = 0; cell_column < cells_per_dim_; ++cell_column) {
            const SInt cell_id = EncodeCell(cell_column, cell_row);

            std::vector<Vertex> centers;
            GenerateVertices(chunk_id, cell_id, centers);

            for (const Vertex& center: centers) {
                GenerateBallHyperedge(center, chunk_id, cell_id);
            }
        }
    }
}

void HyperRGG2D::GenerateHypergraph() {
    GenerateGeometry();
}

void HyperRGG2D::FinalizeHypergraph(MPI_Comm) {}

} // namespace kagen