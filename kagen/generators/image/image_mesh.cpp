#include "kagen/generators/image/image_mesh.h"

#include "kagen/generators/generator.h"
#include "kagen/generators/image/kargb.h"
#include "kagen/generators/image/weight_models.h"
#include "kagen/kagen.h"

#include <array>
#include <cmath>
#include <iostream>

namespace kagen {
namespace {
constexpr bool kDebug = false;
}

PGeneratorConfig
ImageMeshFactory::NormalizeParameters(PGeneratorConfig config, PEID, const PEID size, const bool output) const {
    ImageMeshConfig& iconfig = config.image_mesh;

    {
        bool exists = false;
        bool kargb  = false;
        CheckKARGB(iconfig.filename, exists, kargb);
        if (!exists) {
            throw ConfigurationError("input file does not exist");
        }
        if (!kargb) {
            throw ConfigurationError("input file is not in KARGB format");
        }
    }

    if (iconfig.grid_x == 0 && iconfig.max_grid_x == 0) {
        iconfig.max_grid_x = std::sqrt(size);
    }
    if (iconfig.grid_y == 0 && iconfig.max_grid_y == 0) {
        iconfig.max_grid_y = size / std::sqrt(size);
    }

    // Use the whole grid if not specified otherwise
    if (iconfig.grid_x == 0) {
        iconfig.grid_x = iconfig.max_grid_x;
    } else if (iconfig.max_grid_x == 0) {
        iconfig.max_grid_x = iconfig.grid_x;
    }
    if (iconfig.grid_y == 0) {
        iconfig.grid_y = iconfig.max_grid_y;
    } else if (iconfig.max_grid_y == 0) {
        iconfig.max_grid_y = iconfig.grid_y;
    }

    // Compute number of cols / rows per PE:
    // If either parameter is set, deduce the other one
    // Otherwise, we cut rows and assign multiple cells of just one row to each PE;
    // or, if there are more rows than PEs, we assign whole rows to PEs
    if (iconfig.cols_per_pe == 0 && iconfig.rows_per_pe == 0) {
        iconfig.rows_per_pe = std::max<SInt>(1, iconfig.grid_y / size);
        iconfig.cols_per_pe = iconfig.grid_x / std::max<SInt>(1, size / iconfig.grid_y);
    } else if (iconfig.cols_per_pe == 0) {
        iconfig.cols_per_pe = (iconfig.grid_x * iconfig.grid_y) / (size * iconfig.rows_per_pe);
    } else if (iconfig.rows_per_pe == 0) {
        iconfig.rows_per_pe = (iconfig.grid_x * iconfig.grid_y) / (size * iconfig.cols_per_pe);
    }

    if (output) {
        std::cout << "Grid summary:\n";
        std::cout << "  Divide the image by a " << iconfig.max_grid_x << "x" << iconfig.max_grid_y << " grid";
        if (iconfig.grid_x != iconfig.max_grid_x || iconfig.grid_y != iconfig.max_grid_y) {
            std::cout << ", but only use the top-left " << iconfig.grid_x << "x" << iconfig.grid_y << " subgrid";
        }
        std::cout << "\n";

        std::cout << "  Assign a " << iconfig.cols_per_pe << "x" << iconfig.rows_per_pe << " subgrid to each PE\n";
    }

    // Rectangles must tile the whole grid
    if (size * iconfig.cols_per_pe * iconfig.rows_per_pe != iconfig.grid_x * iconfig.grid_y) {
        throw ConfigurationError("grid too small or too large for the given number of PEs and grid cells per PE");
    }

    // Number of PEs per cell/row must fit
    if (iconfig.grid_x % iconfig.cols_per_pe != 0) {
        throw ConfigurationError(
            "number of used columns must be dividable by the number of columns assigned to each PE");
    }
    if (iconfig.grid_y % iconfig.rows_per_pe != 0) {
        throw ConfigurationError("number of used rows must be dividable by the number of rows assigned to each PE");
    }

    return config;
}

std::unique_ptr<Generator> ImageMeshFactory::Create(const PGeneratorConfig& config, PEID rank, PEID size) const {
    return std::make_unique<ImageMesh>(config, rank, size);
}

namespace {
enum GridDirection {
    //             ULDR
    RIGHT      = 0b0001,
    DOWN_RIGHT = 0b0011,
    DOWN       = 0b0010,
    DOWN_LEFT  = 0b0110,
    LEFT       = 0b0100,
    UP_LEFT    = 0b1100,
    UP         = 0b1000,
    UP_RIGHT   = 0b1001,
    MAX        = 0b1100,
};

struct PEInfo {
    PEInfo() = default;

    PEInfo(
        const PEID rank, const SInt num_pixel_rows, const SInt num_pixel_cols, const SInt rows_per_pe,
        const SInt cols_per_pe, const SInt grid_x, const SInt grid_y, const SInt max_grid_x, const SInt max_grid_y)
        : rank(rank),
          rows_per_pe(rows_per_pe),
          cols_per_pe(cols_per_pe),
          grid_rows(grid_y),
          grid_cols(grid_x) {
        const SInt rows_per_cell     = num_pixel_rows / max_grid_y;
        const SInt rows_per_cell_rem = num_pixel_rows % max_grid_y;
        const SInt cols_per_cell     = num_pixel_cols / max_grid_x;
        const SInt cols_per_cell_rem = num_pixel_cols % max_grid_x;
        const SInt pes_per_row       = grid_x / cols_per_pe;

        grid_start_row = (rank / pes_per_row) * rows_per_pe;
        grid_end_row   = grid_start_row + rows_per_pe;
        grid_start_col = (rank % pes_per_row) * cols_per_pe;
        grid_end_col   = grid_start_col + cols_per_pe;

        total_pixel_rows = grid_y * rows_per_cell + std::min<SInt>(grid_y, rows_per_cell_rem);
        total_pixel_cols = grid_x * cols_per_cell + std::min<SInt>(grid_x, cols_per_cell_rem);

        pixel_start_row = grid_start_row * rows_per_cell + std::min<SInt>(grid_start_row, rows_per_cell_rem);
        pixel_end_row   = (grid_start_row + rows_per_pe) * rows_per_cell
                        + std::min<SInt>(grid_start_row + rows_per_pe, rows_per_cell_rem);
        pixel_start_col = grid_start_col * cols_per_cell + std::min<SInt>(grid_start_col, cols_per_cell_rem);
        pixel_end_col   = (grid_start_col + cols_per_pe) * cols_per_cell
                        + std::min<SInt>(grid_start_col + cols_per_pe, cols_per_cell_rem);

        if constexpr (kDebug) {
            std::cout << "rank=" << rank << ", pixel_start_row=" << pixel_start_row
                      << ", pixel_end_row=" << pixel_end_row << ", pixel_start_col=" << pixel_start_col
                      << ", pixel_end_col=" << pixel_end_col << ", grid_start_row=" << grid_start_row
                      << ", grid_end_row=" << grid_end_row << ", grid_start_col=" << grid_start_col
                      << ", grid_end_col=" << grid_end_col << ", rows_per_cell=" << rows_per_cell
                      << ", rows_per_cell_rem=" << rows_per_cell_rem << ", cols_per_cell=" << cols_per_cell
                      << ", cols_per_cell_rem=" << cols_per_cell_rem << ", pes_per_row=" << pes_per_row << std::endl;
        }
    }

    SSInt NumPixelRows() const {
        return pixel_end_row - pixel_start_row;
    }

    SSInt NumPixelCols() const {
        return pixel_end_col - pixel_start_col;
    }

    SInt LocalPixelToGlobalVertex(const SInt row, const SInt col) const {
        return FirstGlobalPixel() + LocalPixelToLocalVertex(row, col);
    }

    SInt GlobalPixelToGlobalVertex(const SInt row, const SInt col) const {
        const SInt local_row = row - pixel_start_row;
        const SInt local_col = col - pixel_start_col;
        return FirstGlobalPixel() + LocalPixelToLocalVertex(local_row, local_col);
    }

    SInt LocalPixelToLocalVertex(const SInt row, const SInt col) const {
        return row * NumPixelCols() + col;
    }

    SInt NumLocalPixels() const {
        return NumPixelRows() * NumPixelCols();
    }

    SInt FirstGlobalPixel() const {
        return pixel_start_row * total_pixel_cols + NumPixelRows() * pixel_start_col;
    }

    bool IsRightmost() const {
        return grid_end_col == grid_cols;
    }

    bool IsBottommost() const {
        return grid_end_row == grid_rows;
    }

    bool IsLeftmost() const {
        return grid_start_col == 0;
    }

    bool IsTopmost() const {
        return grid_start_row == 0;
    }

    PEID GetNeighboringRank(const GridDirection direction) const {
        const PEID up   = rank - grid_cols / cols_per_pe;
        const PEID down = rank + grid_cols / cols_per_pe;

        switch (direction) {
            case GridDirection::RIGHT:
                return rank + 1;
            case GridDirection::DOWN_RIGHT:
                return down + 1;
            case GridDirection::DOWN:
                return down;
            case GridDirection::DOWN_LEFT:
                return down - 1;
            case GridDirection::LEFT:
                return rank - 1;
            case GridDirection::UP_LEFT:
                return up - 1;
            case GridDirection::UP:
                return up;
            case GridDirection::UP_RIGHT:
                return up + 1;
        }

        return -1;
    }

    PEID rank;

    SInt rows_per_pe;
    SInt cols_per_pe;
    SInt grid_rows;
    SInt grid_cols;

    SInt grid_start_row;
    SInt grid_end_row;
    SInt grid_start_col;
    SInt grid_end_col;

    SInt total_pixel_rows;
    SInt total_pixel_cols;

    SInt pixel_start_row;
    SInt pixel_end_row;
    SInt pixel_start_col;
    SInt pixel_end_col;
};
} // namespace

ImageMesh::ImageMesh(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size) {}

void ImageMesh::GenerateEdgeList() {
    SInt max_pixel_rows, max_pixel_cols;
    std::tie(max_pixel_rows, max_pixel_cols) = ReadDimensions(config_.image_mesh.filename);

    auto make_peinfo = [&](const PEID rank) {
        if (rank < 0 || rank >= size_) {
            return PEInfo();
        }
        return PEInfo(
            rank, max_pixel_rows, max_pixel_cols, config_.image_mesh.rows_per_pe, config_.image_mesh.cols_per_pe,
            config_.image_mesh.grid_x, config_.image_mesh.grid_y, config_.image_mesh.max_grid_x,
            config_.image_mesh.max_grid_y);
    };

    PEInfo my = make_peinfo(rank_);
    SetVertexRange(my.FirstGlobalPixel(), my.FirstGlobalPixel() + my.NumLocalPixels());

    std::array<PEInfo, GridDirection::MAX + 1> neighbors;
    neighbors[GridDirection::RIGHT]      = make_peinfo(my.GetNeighboringRank(GridDirection::RIGHT));
    neighbors[GridDirection::DOWN_RIGHT] = make_peinfo(my.GetNeighboringRank(GridDirection::DOWN_RIGHT));
    neighbors[GridDirection::DOWN]       = make_peinfo(my.GetNeighboringRank(GridDirection::DOWN));
    neighbors[GridDirection::DOWN_LEFT]  = make_peinfo(my.GetNeighboringRank(GridDirection::DOWN_LEFT));
    neighbors[GridDirection::UP_LEFT]    = make_peinfo(my.GetNeighboringRank(GridDirection::UP_LEFT));
    neighbors[GridDirection::UP]         = make_peinfo(my.GetNeighboringRank(GridDirection::UP));
    neighbors[GridDirection::UP_RIGHT]   = make_peinfo(my.GetNeighboringRank(GridDirection::UP_RIGHT));
    neighbors[GridDirection::LEFT]       = make_peinfo(my.GetNeighboringRank(GridDirection::LEFT));

    // If requested, generate coordinates
    if (config_.coordinates) {
        for (SInt cur_row = my.pixel_start_row; cur_row < my.pixel_end_row; cur_row++) {
            for (SInt cur_col = my.pixel_start_col; cur_col < my.pixel_end_col; ++cur_col) {
                PushCoordinate(1.0 * cur_col / max_pixel_cols, 1.0 * cur_row / max_pixel_rows);
            }
        }
    }

    const SSInt     border = (config_.image_mesh.neighborhood > 8) ? 2 : 1;
    const ImageRect img    = ReadRect(
        config_.image_mesh.filename, my.pixel_start_row, my.pixel_start_col, my.pixel_end_row, my.pixel_end_col,
        border);
    const bool diagonal_edges     = config_.image_mesh.neighborhood == 8;
    const bool second_layer_edges = config_.image_mesh.neighborhood == 24;

    auto edge = [&](auto weight_model, const SSInt from_row, const SSInt from_col, const SSInt to_row,
                    const SSInt to_col) {
        const std::uint8_t right = to_col >= my.NumPixelCols();
        const std::uint8_t down  = to_row >= my.NumPixelRows();
        const std::uint8_t left  = to_col < 0;
        const std::uint8_t up    = to_row < 0;

        // Do not generated edges to outside the image
        if ((left && my.IsLeftmost()) || (up && my.IsTopmost()) || (down && my.IsBottommost())
            || (right && my.IsRightmost())) {
            return;
        }

        const RGB&   lhs = img.GetPixel(from_row, from_col);
        const RGB&   rhs = img.GetPixel(to_row, to_col);
        const double weight =
            weight_model(lhs, rhs) * config_.image_mesh.weight_multiplier + config_.image_mesh.weight_offset;

        // Do not generate edges if their weight is below / above the threshold
        if (weight < config_.image_mesh.weight_min_threshold || weight > config_.image_mesh.weight_max_threshold) {
            return;
        }

        const SInt from = my.LocalPixelToGlobalVertex(from_row, from_col);
        const SInt to   = [&] {
            const std::uint8_t direction = (up << 3) | (left << 2) | (down << 1) | right;

            if (direction) {
                // External edge
                const SSInt global_to_row = static_cast<SSInt>(my.pixel_start_row) + to_row;
                const SSInt global_to_col = static_cast<SSInt>(my.pixel_start_col) + to_col;
                return neighbors[direction].GlobalPixelToGlobalVertex(global_to_row, global_to_col);
            } else {
                // Internal edge
                return my.LocalPixelToGlobalVertex(to_row, to_col);
            }
        }();

        PushEdge(from, to);
        PushEdgeWeight(weight);
    };

    auto generate_graph = [&](auto weight_model) {
        for (SSInt row = 0; row < static_cast<SSInt>(my.NumPixelRows()); ++row) {
            for (SSInt col = 0; col < static_cast<SSInt>(my.NumPixelCols()); ++col) {
                edge(weight_model, row, col, row, col + 1);
                edge(weight_model, row, col, row + 1, col);
                edge(weight_model, row, col, row, col - 1);
                edge(weight_model, row, col, row - 1, col);

                if (diagonal_edges) {
                    edge(weight_model, row, col, row + 1, col + 1);
                    edge(weight_model, row, col, row + 1, col - 1);
                    edge(weight_model, row, col, row - 1, col - 1);
                    edge(weight_model, row, col, row - 1, col + 1);
                }

                if (second_layer_edges) {
                    edge(weight_model, row, col, row, col + 2);
                    edge(weight_model, row, col, row + 1, col + 2);
                    edge(weight_model, row, col, row + 2, col + 2);
                    edge(weight_model, row, col, row + 2, col + 1);
                    edge(weight_model, row, col, row + 2, col);
                    edge(weight_model, row, col, row + 2, col - 1);
                    edge(weight_model, row, col, row + 2, col - 2);
                    edge(weight_model, row, col, row + 1, col - 2);
                    edge(weight_model, row, col, row, col - 2);
                    edge(weight_model, row, col, row - 1, col - 2);
                    edge(weight_model, row, col, row - 2, col - 2);
                    edge(weight_model, row, col, row - 2, col - 1);
                    edge(weight_model, row, col, row - 2, col);
                    edge(weight_model, row, col, row - 2, col + 1);
                    edge(weight_model, row, col, row - 2, col + 2);
                    edge(weight_model, row, col, row - 1, col + 2);
                }
            }
        }
    };

    switch (config_.image_mesh.weight_model) {
        case ImageMeshWeightModel::L2:
            generate_graph(L2WeightModel{});
            break;

        case ImageMeshWeightModel::INV_L2:
            generate_graph(InvL2WeightModel{});
            break;

        case ImageMeshWeightModel::RATIO:
            generate_graph(RatioWeightModel{});
            break;

        case ImageMeshWeightModel::INV_RATIO:
            generate_graph(InvRatioWeightModel{});
            break;

        case ImageMeshWeightModel::SIMILARITY:
            generate_graph(SimilarityWeightModel{});
            break;
    }
}
} // namespace kagen
