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

std::vector<RGB> ReadRect(
    const std::string& filename, const SSInt first_row, const SSInt first_col, const SSInt first_invalid_row,
    const SSInt first_invalid_col) {
    const SSInt num_rows = first_invalid_row - first_row;
    const SSInt num_cols = first_invalid_col - first_col;

    std::vector<RGB> pixels;
    pixels.reserve(num_rows * num_cols);

    std::uint64_t num_rows_in_file;
    std::uint64_t num_cols_in_file;
    std::ifstream in(filename, std::ios_base::binary);
    in.seekg(kKargbIdentifierLength * sizeof(char));
    in.read(reinterpret_cast<char*>(&num_rows_in_file), sizeof(std::uint64_t));
    in.read(reinterpret_cast<char*>(&num_cols_in_file), sizeof(std::uint64_t));

    auto push_row = [&](const SInt row, const SInt from_col, const SInt to_col) {
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

    auto push_blank_row = [&](const SInt num_cols) {
        for (SInt cur_col = 0; cur_col < num_cols; ++cur_col) {
            pixels.emplace_back(0, 0, 0);
        }
    };

    SSInt cur_row = first_row;
    for (; cur_row < 0; ++cur_row) {
        push_blank_row(num_cols);
    }
    for (; cur_row < std::min<SSInt>(num_rows, first_invalid_row); ++cur_row) {
        for (SSInt cur_col = first_col; cur_col < 0; ++cur_col) {
            pixels.emplace_back(0, 0, 0);
        }
        push_row(cur_row, first_col, std::min<SSInt>(num_cols, first_invalid_col));
        for (SSInt cur_col = num_cols; cur_col < first_invalid_col; ++cur_col) {
            pixels.emplace_back(0, 0, 0);
        }
    }
    for (; cur_row < num_rows; ++cur_row) {
        push_blank_row(num_cols);
    }

    return pixels;
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

enum GridDirection { RIGHT = 0, DOWN_RIGHT = 1, DOWN = 2, DOWN_LEFT = 3, LEFT = 4, UP_LEFT = 5, UP = 6, UP_RIGHT = 7 };

struct PEInfo {
    PEInfo() = default;

    PEInfo(
        const PEID rank, const SInt num_pixel_rows, const SInt num_pixel_cols, const SInt rows_per_pe,
        const SInt cols_per_pe, const SInt grid_x, const SInt grid_y, const SInt max_grid_x, const SInt max_grid_y)
        : rank(rank),
          max_pixel_rows(num_pixel_rows),
          max_pixel_cols(num_pixel_cols),
          rows_per_pe(rows_per_pe),
          cols_per_pe(cols_per_pe),
          grid_rows(grid_y),
          grid_cols(grid_x),
          max_grid_rows(max_grid_y),
          max_grid_cols(max_grid_x) {
        const SInt rows_per_cell     = num_pixel_rows / max_grid_y;
        const SInt rows_per_cell_rem = num_pixel_rows % max_grid_y;
        const SInt cols_per_cell     = num_pixel_cols / max_grid_x;
        const SInt cols_per_cell_rem = num_pixel_cols % max_grid_x;
        const SInt pes_per_row       = grid_x / cols_per_pe;

        grid_start_row = (rank / pes_per_row) * rows_per_pe;
        grid_end_row   = grid_start_row + rows_per_pe;
        grid_start_col = (rank % pes_per_row) * cols_per_pe;
        grid_end_row   = grid_start_col + 1;

        total_pixel_rows = grid_y * rows_per_cell + std::min<SInt>(grid_y, rows_per_cell_rem);
        total_pixel_cols = grid_x * cols_per_cell + std::min<SInt>(grid_x, cols_per_cell_rem);

        pixel_start_row = grid_start_row * rows_per_cell + std::min<SInt>(grid_start_row, rows_per_cell_rem);
        pixel_end_row   = (grid_start_row + rows_per_pe) * rows_per_cell
                        + std::min<SInt>(grid_start_row + rows_per_pe, rows_per_cell_rem);
        pixel_start_col = grid_start_col * cols_per_cell + std::min<SInt>(grid_start_col, cols_per_cell_rem);
        pixel_end_col   = (grid_start_col * cols_per_pe) * cols_per_cell
                        + std::min<SInt>(grid_start_col + cols_per_pe, cols_per_cell_rem);
    }

    SInt NumPixelRows() const {
        return pixel_end_row - pixel_start_row;
    }

    SInt NumPixelCols() const {
        return pixel_end_col - pixel_start_col;
    }

    SInt VirtualPixelStartRow() const {
        return std::max<SInt>(1, pixel_start_row) - 1;
    }

    SInt VirtualPixelStartCol() const {
        return std::max<SInt>(1, pixel_start_col) - 1;
    }

    SInt VirtualPixelEndRow() const {
        return std::min<SInt>(total_pixel_rows, pixel_end_row + 1);
    }

    SInt VirtualPixelEndCol() const {
        return std::min<SInt>(total_pixel_cols, pixel_end_col + 1);
    }

    SInt NumVirtualPixelRows() const {
        return VirtualPixelEndRow() - VirtualPixelStartRow();
    }

    SInt NumVirtualPixelCols() const {
        return VirtualPixelEndCol() - VirtualPixelStartCol();
    }

    SInt GlobalPixelToLocalVertex(const SInt row, const SInt col) const {
        const SInt local_row = row - pixel_start_row;
        const SInt local_col = col - pixel_start_col;
        return LocalPixelToLocalVertex(local_row, local_col);
    }

    SInt LocalPixelToLocalVertex(const SInt row, const SInt col) const {
        return row * NumPixelCols() + col;
    }

    SInt LocalPixelToGlobalVertex(const SInt row, const SInt col) const {
        return FirstGlobalPixel() + LocalPixelToLocalVertex(row, col);
    }

    SInt GlobalPixelToGlobalVertex(const SInt row, const SInt col) const {
        return FirstGlobalPixel() + GlobalPixelToLocalVertex(row, col);
    }

    SInt NumLocalPixels() const {
        return NumPixelRows() * NumPixelCols();
    }

    SInt FirstGlobalPixel() const {
        return pixel_start_row * total_pixel_cols + NumPixelRows() * pixel_start_col;
    }

    SInt FirstInvalidGlobalPixel() const {
        return FirstGlobalPixel() + NumLocalPixels();
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

    SInt max_pixel_rows;
    SInt max_pixel_cols;
    SInt rows_per_pe;
    SInt cols_per_pe;
    SInt grid_rows;
    SInt grid_cols;
    SInt max_grid_rows;
    SInt max_grid_cols;

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

    // Number neighbors as follows:
    // 5 6 7
    // 4 * 0
    // 3 2 1
    std::array<PEInfo, 8> neighbors;
    neighbors[GridDirection::RIGHT]      = make_peinfo(my.GetNeighboringRank(GridDirection::RIGHT));
    neighbors[GridDirection::DOWN_RIGHT] = make_peinfo(my.GetNeighboringRank(GridDirection::DOWN_RIGHT));
    neighbors[GridDirection::DOWN]       = make_peinfo(my.GetNeighboringRank(GridDirection::DOWN));
    neighbors[GridDirection::DOWN_LEFT]  = make_peinfo(my.GetNeighboringRank(GridDirection::DOWN_LEFT));
    neighbors[GridDirection::LEFT]       = make_peinfo(my.GetNeighboringRank(GridDirection::LEFT));
    neighbors[GridDirection::UP_LEFT]    = make_peinfo(my.GetNeighboringRank(GridDirection::UP_LEFT));
    neighbors[GridDirection::UP]         = make_peinfo(my.GetNeighboringRank(GridDirection::UP));
    neighbors[GridDirection::UP_RIGHT]   = make_peinfo(my.GetNeighboringRank(GridDirection::UP_RIGHT));

    // If requested, generate coordinates
    if (config_.coordinates) {
        for (SInt cur_row = my.pixel_start_row; cur_row < my.pixel_end_row; cur_row++) {
            for (SInt cur_col = my.pixel_start_col; cur_col < my.pixel_end_col; ++cur_col) {
                PushCoordinate(1.0 * cur_col / max_pixel_cols, 1.0 * cur_row / max_pixel_rows);
            }
        }
    }

    const std::vector<RGB> pixels = ReadRect(
        config_.image_mesh.filename, static_cast<SSInt>(my.pixel_start_row) - 1,
        static_cast<SSInt>(my.pixel_start_col) - 1, my.NumPixelRows() + 2, my.NumPixelCols() + 2);
    const bool diagonal_edges = config_.image_mesh.neighborhood == 8;

    auto generate_graph = [&](auto weight_model) {
        auto internal_edge = [&](const SInt from_local_row, const SInt from_local_col, const SInt to_local_row,
                                 const SInt to_local_col) {
            const SInt   from_rgb_pos = (from_local_row + 1) * (my.NumPixelRows() + 2) + (from_local_col + 1);
            const SInt   to_rgb_pos   = (to_local_row + 1) * (my.NumPixelRows() + 2) + (to_local_col + 1);
            const RGB&   from_rgb     = pixels[from_rgb_pos];
            const RGB&   to_rgb       = pixels[to_rgb_pos];
            const double weight       = weight_model(from_rgb, to_rgb);

            if (weight >= config_.image_mesh.weight_min_threshold
                && weight <= config_.image_mesh.weight_max_threshold) {
                PushEdgeWeight(static_cast<SSInt>(weight));

                const SInt from = my.LocalPixelToGlobalVertex(from_local_row, from_local_col);
                const SInt to   = my.LocalPixelToGlobalVertex(to_local_row, to_local_col);
                PushEdge(from, to);
            }
        };
        auto external_edge = [&](const SInt from_local_row, const SInt from_local_col, const GridDirection direction) {
            const SInt from_rgb_pos = (from_local_row + 1) * (my.NumPixelRows() + 2) + (from_local_col + 1);
            const SInt from         = my.LocalPixelToGlobalVertex(from_local_row, from_local_col);

            SInt to_rgb_pos = 0;
            SInt to         = 0;

            switch (direction) {
                case GridDirection::RIGHT:
                    if (my.IsRightmost()) {
                        return;
                    }
                    to_rgb_pos = (from_local_row + 1) * (my.NumPixelRows() + 2) + (from_local_col + 2);
                    to         = neighbors[GridDirection::RIGHT].LocalPixelToGlobalVertex(from_local_row, 0);
                    break;

                case GridDirection::DOWN_RIGHT:
                    if (my.IsRightmost() || my.IsBottommost()) {
                        return;
                    }
                    to_rgb_pos = (from_local_row + 2) * (my.NumPixelRows() + 2) + (from_local_col + 2);
                    to         = neighbors[GridDirection::DOWN_RIGHT].LocalPixelToGlobalVertex(0, 0);
                    break;

                case GridDirection::DOWN:
                    if (my.IsBottommost()) {
                        return;
                    }
                    to_rgb_pos = (from_local_row + 2) * (my.NumPixelRows() + 2) + (from_local_col + 1);
                    to         = neighbors[GridDirection::DOWN].LocalPixelToGlobalVertex(0, from_local_col);
                    break;

                case GridDirection::DOWN_LEFT:
                    if (my.IsLeftmost() || my.IsBottommost()) {
                        return;
                    }
                    to_rgb_pos = (from_local_row + 2) * (my.NumPixelRows() + 2) + from_local_col;
                    to         = neighbors[GridDirection::DOWN_LEFT].LocalPixelToGlobalVertex(0, from_local_col - 1);
                    break;
            }

            const RGB&   from_rgb = pixels[from_rgb_pos];
            const RGB&   to_rgb   = pixels[to_rgb_pos];
            const double weight   = weight_model(from_rgb, to_rgb);

            if (weight >= config_.image_mesh.weight_min_threshold
                && weight <= config_.image_mesh.weight_max_threshold) {
                PushEdgeWeight(static_cast<SSInt>(weight));
                PushEdge(from, to);
            }
        };

        { // Top row
            const SInt row = 0;
            { // Top left corner
                const SInt col = 0;
                internal_edge(row, col, row, col + 1);
                internal_edge(row, col, row + 1, col);
                external_edge(row, col, GridDirection::LEFT);
                external_edge(row, col, GridDirection::UP);
                if (diagonal_edges) {
                    internal_edge(row, col, row + 1, col + 1);
                    external_edge(row, col, GridDirection::DOWN_LEFT);
                    external_edge(row, col, GridDirection::UP_LEFT);
                    external_edge(row, col, GridDirection::UP_RIGHT);
                }
            }
            for (SInt col = 1; col + 1 < my.NumPixelCols(); ++col) { // Top row interior
                internal_edge(row, col, row, col + 1);
                internal_edge(row, col, row + 1, col);
                internal_edge(row, col, row, col - 1);
                external_edge(row, col, GridDirection::UP);
                if (diagonal_edges) {
                    internal_edge(row, col, row + 1, col + 1);
                    internal_edge(row, col, row + 1, col - 1);
                    external_edge(row, col, GridDirection::UP_LEFT);
                    external_edge(row, col, GridDirection::UP_RIGHT);
                }
            }
            { // Top right corner
                const SInt col = my.NumPixelCols() - 1;
                external_edge(row, col, GridDirection::RIGHT);
                internal_edge(row, col, row + 1, col);
                internal_edge(row, col, row, col - 1);
                external_edge(row, col, GridDirection::UP);
                if (diagonal_edges) {
                    external_edge(row, col, GridDirection::DOWN_RIGHT);
                    internal_edge(row, col, row - 1, col - 1);
                    external_edge(row, col, GridDirection::UP_LEFT);
                    external_edge(row, col, GridDirection::UP_RIGHT);
                }
            }
        }
        // Middle part
        for (SInt row = 1; row + 1 < my.NumPixelRows(); ++row) {
            // Leftmost column
            // Middle
            for (SInt col = 1; col + 1 < my.NumPixelCols(); ++col) {
            }
            // Rightmost column
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
    /*
        auto generate_edge = [&](const auto& weight_model, const SInt row1, const SInt col1, const SInt row2,
                                 const SInt col2) {
            const RGB&   rgb1 = pixels[(row1 - my_virtual_start_row) * my_num_virtual_rows + (col1 -
       my_virtual_start_col)]; const RGB&   rgb2 = pixels[(row2 - my_virtual_start_row) * my_num_virtual_rows + (col2 -
       my_virtual_start_col)]; const double weight = weight_model(rgb1, rgb2);

            const SInt vertex1 = row1 * num_cols + col1;
            const SInt vertex2 = row2 * num_cols + col2;

            if constexpr (kDebug) {
                std::cout << vertex1 << "[" << static_cast<int>(rgb1.r) << "," << static_cast<int>(rgb1.g) << ","
                          << static_cast<int>(rgb1.b) << "] --> " << vertex2 << "[" << static_cast<int>(rgb2.r) << ","
                          << static_cast<int>(rgb2.g) << "," << static_cast<int>(rgb2.b) << "] = " << weight <<
       std::endl;
            }

            if (weight >= config_.image_mesh.weight_min_threshold && weight <= config_.image_mesh.weight_max_threshold)
       { PushEdgeWeight(static_cast<SSInt>(weight_model(rgb1, rgb2))); PushEdge(vertex1, vertex2);
            }
        };

        auto generate_edges = [&](auto weight_model) {
            // Special treatment for first row
            if (my_start_row == 0 && my_end_row != 0) { // && -> catch special case where PE 0 gets no vertices
                const SInt cur_row = my_start_row;
                if (my_start_col == 0 && my_end_col != 0) {
                    const SInt cur_col = my_start_col;
                    generate_edge(weight_model, cur_row, cur_col, cur_row, cur_col + 1);
                    generate_edge(weight_model, cur_row, cur_col, cur_row + 1, cur_col);
                    if (config_.image_mesh.neighborhood == 8) {
                        generate_edge(weight_model, cur_row, cur_col, cur_row + 1, cur_col + 1);
                    }
                }
                for (SInt cur_col = my_virtual_start_col + 1; cur_col + 1 < my_virtual_end_col; ++cur_col) {
                    generate_edge(weight_model, cur_row, cur_col, cur_row, cur_col + 1);
                    generate_edge(weight_model, cur_row, cur_col, cur_row, cur_col - 1);
                    generate_edge(weight_model, cur_row, cur_col, cur_row + 1, cur_col);
                    if (config_.image_mesh.neighborhood == 8) {
                        generate_edge(weight_model, cur_row, cur_col, cur_row + 1, cur_col + 1);
                        generate_edge(weight_model, cur_row, cur_col, cur_row + 1, cur_col - 1);
                    }
                }
                if (my_end_col == num_cols && my_end_col > 0) {
                    const SInt cur_col = my_end_col - 1;
                    generate_edge(weight_model, cur_row, cur_col, cur_row + 1, cur_col);
                    generate_edge(weight_model, cur_row, cur_col, cur_row, cur_col - 1);
                    if (config_.image_mesh.neighborhood == 8) {
                        generate_edge(weight_model, cur_row, cur_col, cur_row + 1, cur_col - 1);
                    }
                }
            }

            // Special treatment for first column
            if (my_start_col == 0 && my_end_col > 0) {
                const SInt cur_col = my_start_col;
                for (SInt cur_row = my_virtual_start_row + 1; cur_row + 1 < my_virtual_end_row; ++cur_row) {
                    generate_edge(weight_model, cur_row, cur_col, cur_row, cur_col + 1);
                    generate_edge(weight_model, cur_row, cur_col, cur_row + 1, cur_col);
                    generate_edge(weight_model, cur_row, cur_col, cur_row - 1, cur_col);
                    if (config_.image_mesh.neighborhood == 8) {
                        generate_edge(weight_model, cur_row, cur_col, cur_row + 1, cur_col + 1);
                        generate_edge(weight_model, cur_row, cur_col, cur_row - 1, cur_col + 1);
                    }
                }
            }

            // Special treatment for last column
            if (my_end_col == num_cols && my_start_col < num_cols) {
                const SInt cur_col = my_end_col - 1;
                for (SInt cur_row = my_virtual_start_row + 1; cur_row + 1 < my_virtual_end_row; ++cur_row) {
                    generate_edge(weight_model, cur_row, cur_col, cur_row + 1, cur_col);
                    generate_edge(weight_model, cur_row, cur_col, cur_row, cur_col - 1);
                    generate_edge(weight_model, cur_row, cur_col, cur_row - 1, cur_col);
                    if (config_.image_mesh.neighborhood == 8) {
                        generate_edge(weight_model, cur_row, cur_col, cur_row + 1, cur_col - 1);
                        generate_edge(weight_model, cur_row, cur_col, cur_row - 1, cur_col - 1);
                    }
                }
            }

            // Special treatment for last row
            if (my_end_row == num_rows && my_start_row < num_rows) {
                const SInt cur_row = my_end_row - 1;
                if (my_start_col == 0 && my_end_col != 0) {
                    const SInt cur_col = my_start_col;
                    generate_edge(weight_model, cur_row, cur_col, cur_row, cur_col + 1);
                    generate_edge(weight_model, cur_row, cur_col, cur_row - 1, cur_col);
                    if (config_.image_mesh.neighborhood == 8) {
                        generate_edge(weight_model, cur_row, cur_col, cur_row - 1, cur_col + 1);
                    }
                }
                for (SInt cur_col = my_virtual_start_col + 1; cur_col + 1 < my_virtual_end_col; ++cur_col) {
                    generate_edge(weight_model, cur_row, cur_col, cur_row, cur_col + 1);
                    generate_edge(weight_model, cur_row, cur_col, cur_row, cur_col - 1);
                    generate_edge(weight_model, cur_row, cur_col, cur_row - 1, cur_col);
                    if (config_.image_mesh.neighborhood == 8) {
                        generate_edge(weight_model, cur_row, cur_col, cur_row - 1, cur_col + 1);
                        generate_edge(weight_model, cur_row, cur_col, cur_row - 1, cur_col - 1);
                    }
                }
                if (my_end_col == num_cols && my_end_col > 0) {
                    const SInt cur_col = my_end_col - 1;
                    generate_edge(weight_model, cur_row, cur_col, cur_row - 1, cur_col);
                    generate_edge(weight_model, cur_row, cur_col, cur_row, cur_col - 1);
                    if (config_.image_mesh.neighborhood == 8) {
                        generate_edge(weight_model, cur_row, cur_col, cur_row - 1, cur_col - 1);
                    }
                }
            }

            for (SInt cur_row = my_virtual_start_row + 1; cur_row + 1 < my_virtual_end_row; ++cur_row) {
                for (SInt cur_col = my_virtual_start_col + 1; cur_col + 1 < my_virtual_end_col; ++cur_col) {
                    generate_edge(weight_model, cur_row, cur_col, cur_row, cur_col + 1);
                    generate_edge(weight_model, cur_row, cur_col, cur_row + 1, cur_col);
                    generate_edge(weight_model, cur_row, cur_col, cur_row, cur_col - 1);
                    generate_edge(weight_model, cur_row, cur_col, cur_row - 1, cur_col);

                    if (config_.image_mesh.neighborhood == 8) {
                        generate_edge(weight_model, cur_row, cur_col, cur_row + 1, cur_col + 1);
                        generate_edge(weight_model, cur_row, cur_col, cur_row + 1, cur_col - 1);
                        generate_edge(weight_model, cur_row, cur_col, cur_row - 1, cur_col - 1);
                        generate_edge(weight_model, cur_row, cur_col, cur_row - 1, cur_col + 1);
                    }
                }
            }
        };

        switch (config_.image_mesh.weight_model) {
            case ImageMeshWeightModel::L2:
                generate_edges(L2WeightModel{});
                break;

            case ImageMeshWeightModel::INV_L2:
                generate_edges(InvL2WeightModel{});
                break;

            case ImageMeshWeightModel::INV_RATIO:
                generate_edges(InvRatioWeightModel{});
                break;
        }*/
}
} // namespace kagen

