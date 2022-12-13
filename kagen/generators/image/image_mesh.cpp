#include "kagen/generators/image/image_mesh.h"

#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>

#include "kagen/definitions.h"
#include "kagen/generators/generator.h"

namespace kagen {
namespace {
constexpr bool kDebug = false;
}

PGeneratorConfig ImageMeshFactory::NormalizeParameters(PGeneratorConfig config, PEID size, bool output) const {
    ImageMeshConfig& iconfig = config.image_mesh;

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
        std::cout << "  Divide the image by a " << iconfig.max_grid_x << "x" << iconfig.max_grid_y << " grid\n";
        if (iconfig.grid_x != iconfig.max_grid_x || iconfig.grid_y != iconfig.max_grid_y) {
            std::cout << "  -> but only use the top-left " << iconfig.grid_x << "x" << iconfig.grid_y << " subgrid\n";
        }
        std::cout << "  Assign a " << iconfig.cols_per_pe << "x" << iconfig.rows_per_pe << " subgrid to each PE\n";
    }

    // Rectangles must tile the whole grid
    if (size * iconfig.cols_per_pe * iconfig.rows_per_pe != iconfig.grid_x * iconfig.grid_y) {
        throw ConfigurationError("PE rectangles do not cover the whole grid");
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
struct RGB {
    RGB() = default;
    RGB(const std::uint8_t r, const std::uint8_t g, const std::uint8_t b) : r(r), g(g), b(b) {}
    std::uint8_t r;
    std::uint8_t g;
    std::uint8_t b;
};

using Pixels                                 = std::vector<RGB>;
constexpr std::size_t kKargbIdentifierLength = 5;
constexpr std::size_t kKargbHeaderLength     = kKargbIdentifierLength + 2 * sizeof(std::uint64_t);

std::pair<SInt, SInt> ReadDimensions(const std::string& filename) {
    std::uint64_t                                rows;
    std::uint64_t                                cols;
    std::array<char, kKargbIdentifierLength + 1> identifier;

    std::ifstream in(filename, std::ios_base::binary);
    in.read(identifier.data(), kKargbIdentifierLength * sizeof(char));
    in.read(reinterpret_cast<char*>(&rows), sizeof(std::uint64_t));
    in.read(reinterpret_cast<char*>(&cols), sizeof(std::uint64_t));
    identifier[kKargbIdentifierLength] = 0;

    if (std::strcmp(identifier.data(), "KARGB")) {
        std::cerr << "Error: invalid input file; use tools/img2kargb to convert input image\n";
        std::exit(1);
    }

    return {rows, cols};
}

class ImageRect {
public:
    ImageRect(Pixels pixels, const SInt num_cols, const SInt overlap)
        : pixels_(std::move(pixels)),
          num_cols_(num_cols),
          overlap_(overlap) {}

    const RGB& GetPixel(const SSInt row, const SSInt col) const {
        return pixels_[(row + overlap_) * (num_cols_ + 2 * overlap_) + (col + overlap_)];
    }

private:
    Pixels pixels_;
    SInt   num_cols_;
    SInt   overlap_;
};

ImageRect ReadRect(
    const std::string& filename, const SInt from_row, const SInt from_col, const SInt to_row, const SInt to_col,
    const SInt overlap) {
    const SSInt actual_from_row = static_cast<SSInt>(from_row) - static_cast<SSInt>(overlap);
    const SSInt actual_from_col = static_cast<SSInt>(from_col) - static_cast<SSInt>(overlap);
    const SSInt actual_to_row   = static_cast<SSInt>(to_row + overlap);
    const SSInt actual_to_col   = static_cast<SSInt>(to_col + overlap);
    const SSInt actual_num_rows = actual_to_row - actual_from_row;
    const SSInt actual_num_cols = actual_to_col - actual_from_col;

    std::vector<RGB> pixels;
    pixels.reserve(actual_num_rows * actual_num_cols);

    std::uint64_t num_rows_in_file;
    std::uint64_t num_cols_in_file;
    std::ifstream in(filename, std::ios_base::binary);
    in.seekg(kKargbIdentifierLength * sizeof(char));
    in.read(reinterpret_cast<char*>(&num_rows_in_file), sizeof(std::uint64_t));
    in.read(reinterpret_cast<char*>(&num_cols_in_file), sizeof(std::uint64_t));

    auto push_row = [&in, &pixels, num_cols_in_file](const SInt row, const SInt from_col, const SInt to_col) {
        const SInt row_start_pos = row * num_cols_in_file;
        const SInt col_start_pos = row_start_pos + from_col;
        in.seekg(kKargbHeaderLength + col_start_pos * 3 * sizeof(std::uint8_t));
        for (SInt cur_col = from_col; cur_col < to_col; ++cur_col) {
            std::uint8_t r, g, b;
            in.read(reinterpret_cast<char*>(&r), sizeof(std::uint8_t));
            in.read(reinterpret_cast<char*>(&g), sizeof(std::uint8_t));
            in.read(reinterpret_cast<char*>(&b), sizeof(std::uint8_t));
            pixels.emplace_back(r, g, b);
        }
    };

    auto push_blank_row = [&pixels](const SInt num_cols) {
        for (SInt cur_col = 0; cur_col < num_cols; ++cur_col) {
            pixels.emplace_back(0, 0, 0);
        }
    };

    SSInt cur_row = actual_from_row;
    for (; cur_row < 0; ++cur_row) {
        push_blank_row(actual_num_cols);
    }
    for (; cur_row < std::min<SSInt>(num_rows_in_file, actual_to_row); ++cur_row) {
        for (SSInt cur_col = actual_from_col; cur_col < 0; ++cur_col) {
            pixels.emplace_back(0, 0, 0);
        }
        push_row(cur_row, std::max<SSInt>(0, actual_from_col), std::min<SSInt>(num_cols_in_file, actual_to_col));
        for (SSInt cur_col = num_cols_in_file; cur_col < actual_to_col; ++cur_col) {
            pixels.emplace_back(0, 0, 0);
        }
    }
    for (; cur_row < actual_to_row; ++cur_row) {
        push_blank_row(actual_num_cols);
    }

    return {std::move(pixels), static_cast<SInt>(to_col - from_col), overlap};
}

std::uint8_t Delta(const std::uint8_t lhs, const std::uint8_t rhs) {
    return std::max(lhs, rhs) - std::min(lhs, rhs);
}

struct L2WeightModel {
    double operator()(const RGB& lhs, const RGB& rhs) const {
        const std::uint8_t dr = Delta(lhs.r, rhs.r);
        const std::uint8_t dg = Delta(lhs.g, rhs.g);
        const std::uint8_t db = Delta(lhs.b, rhs.b);
        return std::sqrt(dr * dr + dg * dg + db * db);
    }
};

struct InvL2WeightModel {
    double operator()(const RGB& lhs, const RGB& rhs) const {
        return max_value_ - l2_(lhs, rhs);
    }

private:
    L2WeightModel l2_{};
    double        max_value_ = 255 * std::sqrt(3) + 1;
};

double MaxMinRatio(const std::uint8_t lhs, const std::uint8_t rhs) {
    return 1.0 * std::max(lhs, rhs) / std::min(lhs, rhs);
}

struct InvRatioWeightModel {
    double operator()(const RGB& lhs, const RGB& rhs) const {
        return 1.0 / MaxMinRatio(lhs.r, rhs.r) * 1.0 / MaxMinRatio(lhs.g, rhs.g) * 1.0 / MaxMinRatio(lhs.b, rhs.b);
    }
};

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
        grid_end_col   = grid_start_col + 1;

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
                      << ", pixel_end_col=" << pixel_end_col << "grid_start_row=" << grid_start_row
                      << ", grid_end_row=" << grid_end_row << ", grid_start_col=" << grid_start_col
                      << ", grid_end_col=" << grid_end_col << "rows_per_cell=" << rows_per_cell
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

void ImageMesh::GenerateImpl() {
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

    std::array<PEInfo, GridDirection::MAX> neighbors;
    if (!my.IsRightmost()) {
        neighbors[GridDirection::RIGHT] = make_peinfo(my.GetNeighboringRank(GridDirection::RIGHT));
    }
    if (!my.IsBottommost()) {
        if (!my.IsRightmost()) {
            neighbors[GridDirection::DOWN_RIGHT] = make_peinfo(my.GetNeighboringRank(GridDirection::DOWN_RIGHT));
        }
        neighbors[GridDirection::DOWN] = make_peinfo(my.GetNeighboringRank(GridDirection::DOWN));
        if (!my.IsLeftmost()) {
            neighbors[GridDirection::DOWN_LEFT] = make_peinfo(my.GetNeighboringRank(GridDirection::DOWN_LEFT));
        }
    }
    if (!my.IsTopmost()) {
        if (!my.IsLeftmost()) {
            neighbors[GridDirection::UP_LEFT] = make_peinfo(my.GetNeighboringRank(GridDirection::UP_LEFT));
        }
        neighbors[GridDirection::UP] = make_peinfo(my.GetNeighboringRank(GridDirection::UP));
        if (!my.IsRightmost()) {
            neighbors[GridDirection::UP_RIGHT] = make_peinfo(my.GetNeighboringRank(GridDirection::UP_RIGHT));
        }
    }
    if (!my.IsLeftmost()) {
        neighbors[GridDirection::LEFT] = make_peinfo(my.GetNeighboringRank(GridDirection::LEFT));
    }

    // If requested, generate coordinates
    if (config_.coordinates) {
        for (SInt cur_row = my.pixel_start_row; cur_row < my.pixel_end_row; cur_row++) {
            for (SInt cur_col = my.pixel_start_col; cur_col < my.pixel_end_col; ++cur_col) {
                PushCoordinate(1.0 * cur_col / max_pixel_cols, 1.0 * cur_row / max_pixel_rows);
            }
        }
    }

    const ImageRect img = ReadRect(
        config_.image_mesh.filename, my.pixel_start_row, my.pixel_start_col, my.pixel_end_row, my.pixel_end_col, 1);
    const bool diagonal_edges = config_.image_mesh.neighborhood == 8;

    auto edge = [&](auto weight_model, const SSInt from_row, const SSInt from_col, const SSInt to_row,
                    const SSInt to_col) {
        const std::uint8_t right = to_col >= my.NumPixelCols();
        const std::uint8_t down  = to_row >= my.NumPixelRows();
        const std::uint8_t left  = to_col < 0;
        const std::uint8_t up    = to_row < 0;

        // Do not generated edges to outside the image
        if ((left && my.IsLeftmost()) || (up && my.IsTopmost()) || (down && my.IsBottommost())
            || (right && my.IsRightmost())) {
            // std::cout << from_row << "," << from_col << " --> " << to_row << "," << to_col << " is cut edge"
            //           << std::endl;
            return;
        }

        const RGB&   lhs    = img.GetPixel(from_row, from_col);
        const RGB&   rhs    = img.GetPixel(to_row, to_col);
        const double weight = weight_model(lhs, rhs);

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

        // std::cout << from << " --> " << to << std::endl;

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

        case ImageMeshWeightModel::INV_RATIO:
            generate_graph(InvRatioWeightModel{});
            break;
    }
}
} // namespace kagen

