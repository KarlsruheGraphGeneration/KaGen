#include "kagen/generators/image/image_mesh.h"

#include <cstring>

#include "kagen/definitions.h"
#include "kagen/generators/generator.h"

namespace kagen {
PGeneratorConfig ImageMeshFactory::NormalizeParameters(PGeneratorConfig config, PEID size, bool output) const {
    ImageMeshConfig& iconfig = config.image_mesh;

    // Use the whole grid if not specified otherwise
    if (iconfig.grid_x == 0) {
        iconfig.grid_x = iconfig.max_grid_x;
    }
    if (iconfig.grid_y == 0) {
        iconfig.grid_y = iconfig.max_grid_y;
    }

    // Compute number of cols / rows per PE: 
    // If either parameter is set, deduce the other one 
    // Otherwise, we cut rows and assign multiple cells of just one row to each PE; 
    // or, if there are more rows than PEs, we assign whole rows to PEs
    bool output_pe_rect = true;
    if (iconfig.cols_per_pe == 0 && iconfig.rows_per_pe == 0) {
        iconfig.rows_per_pe = std::max<SInt>(1, iconfig.grid_y / size);
        iconfig.cols_per_pe = iconfig.grid_x / std::max<SInt>(1, size / iconfig.grid_y);
    } else if (iconfig.cols_per_pe == 0) {
        iconfig.cols_per_pe = (iconfig.grid_x * iconfig.grid_y) / (size * iconfig.rows_per_pe);
    } else if (iconfig.rows_per_pe == 0) {
        iconfig.rows_per_pe = (iconfig.grid_x * iconfig.grid_y) / (size * iconfig.cols_per_pe);
    } else {
        output_pe_rect = false;
    }
    if (output_pe_rect && output) {
        std::cout << "Assigning " << iconfig.cols_per_pe << "x" << iconfig.rows_per_pe << " tiles of the grid to PEs\n";
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

ImageMesh::ImageMesh(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size) {}

void ImageMesh::GenerateImpl() {}
} // namespace kagen

