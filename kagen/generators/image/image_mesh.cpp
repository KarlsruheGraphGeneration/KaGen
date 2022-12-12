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
constexpr bool kDebug = true;
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

constexpr std::size_t kKargbHeaderLength = 5;

std::pair<SInt, SInt> ReadDimensions(const std::string& filename) {
    std::uint64_t                            rows;
    std::uint64_t                            cols;
    std::array<char, kKargbHeaderLength + 1> identifier;

    std::ifstream in(filename, std::ios_base::binary);
    in.read(identifier.data(), kKargbHeaderLength * sizeof(char));
    in.read(reinterpret_cast<char*>(&rows), sizeof(std::uint64_t));
    in.read(reinterpret_cast<char*>(&cols), sizeof(std::uint64_t));
    identifier[kKargbHeaderLength] = 0;

    if (std::strcmp(identifier.data(), "KARGB")) {
        std::cerr << "Error: invalid input file; use tools/img2kargb to convert input image\n";
        std::exit(1);
    }

    return {rows, cols};
}

std::vector<RGB>
ReadRect(const std::string& filename, const SInt row, const SInt col, const SInt num_rows, const SInt num_cols) {
    std::vector<RGB> pixels;
    pixels.reserve(num_rows * num_cols);

    std::uint64_t num_rows_in_file;
    std::uint64_t num_cols_in_file;

    std::ifstream in(filename, std::ios_base::binary);
    in.seekg(kKargbHeaderLength * sizeof(char));
    in.read(reinterpret_cast<char*>(&num_rows_in_file), sizeof(std::uint64_t));
    in.read(reinterpret_cast<char*>(&num_cols_in_file), sizeof(std::uint64_t));

    for (SInt cur_row = row; cur_row < row + num_rows; ++cur_row) {
        const SInt row_start_pos = cur_row * num_cols_in_file;
        const SInt col_start_pos = row_start_pos + col;
        in.seekg(col_start_pos * 3 * sizeof(std::uint8_t));

        for (SInt cur_col = col; cur_col < col + num_cols; ++cur_col) {
            std::uint8_t r, g, b;
            in.read(reinterpret_cast<char*>(&r), sizeof(std::uint8_t));
            in.read(reinterpret_cast<char*>(&g), sizeof(std::uint8_t));
            in.read(reinterpret_cast<char*>(&b), sizeof(std::uint8_t));
            pixels.emplace_back(r, g, b);
        }
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
} // namespace

ImageMesh::ImageMesh(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size) {}

void ImageMesh::GenerateImpl() {
    SInt num_rows, num_cols;
    std::tie(num_rows, num_cols) = ReadDimensions(config_.image_mesh.filename);

    const SInt rows_per_cell     = num_rows / config_.image_mesh.max_grid_y;
    const SInt rows_per_cell_rem = num_rows % config_.image_mesh.max_grid_y;
    const SInt cols_per_cell     = num_cols / config_.image_mesh.max_grid_x;
    const SInt cols_per_cell_rem = num_cols % config_.image_mesh.max_grid_x;
    const SInt my_start_grid_row =
        rank_ / (config_.image_mesh.max_grid_x / config_.image_mesh.cols_per_pe) * config_.image_mesh.rows_per_pe;
    const SInt my_start_grid_col =
        rank_ % (config_.image_mesh.max_grid_x / config_.image_mesh.cols_per_pe) * config_.image_mesh.cols_per_pe;

    // Compute the first and last row / col that is read by this PE
    const SInt my_start_row = my_start_grid_row * rows_per_cell + std::min<SInt>(my_start_grid_row, rows_per_cell_rem);
    const SInt my_start_col = my_start_grid_col * cols_per_cell + std::min<SInt>(my_start_grid_col, cols_per_cell_rem);
    const SInt my_end_row   = (my_start_grid_row + config_.image_mesh.rows_per_pe) * rows_per_cell
                            + std::min<SInt>(my_start_grid_row + config_.image_mesh.rows_per_pe, rows_per_cell_rem);
    const SInt my_end_col = (my_start_grid_col + config_.image_mesh.cols_per_pe) * cols_per_cell
                            + std::min<SInt>(my_start_grid_col + config_.image_mesh.cols_per_pe, cols_per_cell_rem);

    // Compute our vertex range
    const SInt my_first_pixel         = my_start_row * num_cols + my_start_col;
    const SInt my_first_invalid_pixel = my_end_row * num_cols + my_end_col;
    SetVertexRange(my_first_pixel, my_first_invalid_pixel);

    // If we want coordinates, use pixel positions as coordinates
    if (config_.coordinates) {
        for (SInt cur_row = my_start_row; cur_row < my_end_row; ++cur_row) {
            for (SInt cur_col = my_start_col; cur_col < my_end_col; ++cur_col) {
                PushCoordinate(1.0 * cur_col / num_cols, 1.0 * cur_row / num_rows);
            }
        }
    }

    // If we are not at the border, overlap our rect by one row/col with pixels from neighboring PEs
    const SInt my_virtual_start_row = std::max<SInt>(1, my_start_row) - 1;
    const SInt my_virtual_end_row   = std::min<SInt>(num_rows - 1, my_end_row) + 1;
    const SInt my_virtual_start_col = std::max<SInt>(1, my_start_col) - 1;
    const SInt my_virtual_end_col   = std::min<SInt>(num_cols - 1, my_end_col) + 1;
    const SInt my_num_virtual_rows  = my_virtual_end_row - my_virtual_start_row;
    const SInt my_num_virtual_cols  = my_virtual_end_col - my_virtual_start_col;

    if constexpr (kDebug) {
        std::cout << "PE " << rank_ << "/" << size_ << ": " << my_virtual_start_col << "x" << my_virtual_start_row
                  << " <-> " << my_virtual_end_col << "x" << my_virtual_end_row << std::endl;
    }

    const std::vector<RGB> pixels = ReadRect(
        config_.image_mesh.filename, my_virtual_start_row, my_virtual_start_col, my_num_virtual_rows,
        my_num_virtual_cols);

    auto generate_edge = [&](const auto& weight_model, const SInt row1, const SInt col1, const SInt row2,
                             const SInt col2) {
        const RGB&   rgb1 = pixels[(row1 - my_virtual_start_row) * my_num_virtual_rows + (col1 - my_virtual_start_col)];
        const RGB&   rgb2 = pixels[(row2 - my_virtual_start_row) * my_num_virtual_rows + (col2 - my_virtual_start_col)];
        const double weight = weight_model(rgb1, rgb2);

        if (weight >= config_.image_mesh.weight_min_threshold && weight <= config_.image_mesh.weight_max_threshold) {
            PushEdgeWeight(static_cast<SSInt>(weight_model(rgb1, rgb2)));

            const SInt vertex1 = row1 * num_cols + col1;
            const SInt vertex2 = row2 * num_cols + col2;
            PushEdge(vertex1, vertex2);
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
                if (config_.image_mesh.neighborhood) {
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
    }
}
} // namespace kagen

