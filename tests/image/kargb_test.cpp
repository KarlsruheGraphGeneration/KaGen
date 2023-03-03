#include <gtest/gtest.h>

#include <cmath>

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/image/image_mesh.h"
#include "kagen/generators/image/kargb.h"

using namespace kagen;

const char* CHECKERBOARD = "tests/data/image/rgcheckerboard.kargb";

TEST(KARGB, reads_checkerboard_dimension) {
    const auto [nrows, ncols] = ReadDimensions(CHECKERBOARD);
    EXPECT_EQ(nrows, 16);
    EXPECT_EQ(ncols, 16);
}

TEST(KARGB, reads_checkerboard) {
    const SInt nrows = 16;
    const SInt ncols = 16;

    // G R G R G R G R
    // R G R G R G R G
    // ...
    const auto pixels = ReadRect(CHECKERBOARD, 0, 0, nrows, ncols, 0);
    for (SInt row = 0; row < nrows; ++row) {
        for (SInt col = 0; col < ncols; ++col) {
            const RGB  pixel            = pixels.GetPixel(row, col);
            const SInt checkerboard_row = row / 2;
            const SInt checkerboard_col = col / 2;
            if ((checkerboard_row + checkerboard_col) % 2 == 0) {
                EXPECT_EQ(pixel.r, 0);
                EXPECT_EQ(pixel.g, 255);
                EXPECT_EQ(pixel.b, 0);
            } else {
                EXPECT_EQ(pixel.r, 255);
                EXPECT_EQ(pixel.g, 0);
                EXPECT_EQ(pixel.b, 0);
            }
        }
    }
}

TEST(KARGB, reads_inner_checkerboard) {
    const SInt from_row = 2;
    const SInt from_col = 2;
    const SInt to_row   = 14;
    const SInt to_col   = 14;

    // G R G R G R
    // R G R G R G
    // ...
    const auto pixels = ReadRect(CHECKERBOARD, from_row, from_col, to_row, to_col, 0);

    for (SInt row = 0; row < to_row - from_row; ++row) {
        for (SInt col = 0; col < to_col - from_col; ++col) {
            const RGB  pixel            = pixels.GetPixel(row, col);
            const SInt checkerboard_row = row / 2;
            const SInt checkerboard_col = col / 2;
            if ((checkerboard_row + checkerboard_col) % 2 == 0) {
                EXPECT_EQ(pixel.r, 0);
                EXPECT_EQ(pixel.g, 255);
                EXPECT_EQ(pixel.b, 0);
            } else {
                EXPECT_EQ(pixel.r, 255);
                EXPECT_EQ(pixel.g, 0);
                EXPECT_EQ(pixel.b, 0);
            }
        }
    }
}

TEST(KARGB, reads_checkerboard_with_overlap) {
    const SSInt from_row = 2;
    const SSInt from_col = 2;
    const SSInt to_row   = 14;
    const SSInt to_col   = 14;

    // G R G R G R G R
    // R G R G R G R G
    // ...
    const auto pixels = ReadRect(CHECKERBOARD, from_row, from_col, to_row, to_col, 2);
    for (SSInt row = -2; row < (to_row - from_row) + 2; ++row) {
        for (SSInt col = -2; col < (to_col - from_col) + 2; ++col) {
            const RGB  pixel            = pixels.GetPixel(row, col);
            const SInt checkerboard_row = (row + 2) / 2;
            const SInt checkerboard_col = (col + 2) / 2;
            if ((checkerboard_row + checkerboard_col) % 2 == 0) {
                EXPECT_EQ(pixel.r, 0);
                EXPECT_EQ(pixel.g, 255);
                EXPECT_EQ(pixel.b, 0);
            } else {
                EXPECT_EQ(pixel.r, 255);
                EXPECT_EQ(pixel.g, 0);
                EXPECT_EQ(pixel.b, 0);
            }
        }
    }
}

TEST(KARGB, reads_checkerboard_with_virtual_overlap) {
    const SSInt from_row = 0;
    const SSInt from_col = 0;
    const SSInt to_row   = 2;
    const SSInt to_col   = 2;

    // G R G R G R G R
    // R G R G R G R G
    // ...
    const auto pixels = ReadRect(CHECKERBOARD, from_row, from_col, to_row, to_col, 2);
    for (SSInt row = 0; row < (to_row - from_row) + 2; ++row) {
        for (SSInt col = 0; col < (to_col - from_col) + 2; ++col) {
            const RGB  pixel            = pixels.GetPixel(row, col);
            const SInt checkerboard_row = (row + 2) / 2;
            const SInt checkerboard_col = (col + 2) / 2;
            if ((checkerboard_row + checkerboard_col) % 2 == 0) {
                EXPECT_EQ(pixel.r, 0);
                EXPECT_EQ(pixel.g, 255);
                EXPECT_EQ(pixel.b, 0);
            } else {
                EXPECT_EQ(pixel.r, 255);
                EXPECT_EQ(pixel.g, 0);
                EXPECT_EQ(pixel.b, 0);
            }
        }
    }
}

TEST(KARGB, generates_graph_on_one_PE) {
    PGeneratorConfig config;
    config.image_mesh.filename             = CHECKERBOARD;
    config.image_mesh.grid_x               = 1;
    config.image_mesh.grid_y               = 1;
    config.image_mesh.neighborhood         = 4;
    config.image_mesh.weight_model         = ImageMeshWeightModel::L2;
    config.image_mesh.weight_multiplier    = 1;
    config.image_mesh.weight_min_threshold = 1;

    ImageMeshFactory factory;
    config         = factory.NormalizeParameters(config, 0, 1, false);
    auto generator = factory.Create(config, 0, 1);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    const auto result = generator->Take();

    // Vertices: 16x16 image -> 256
    ASSERT_EQ(result.vertex_range.second - result.vertex_range.first, 16 * 16);

    // Edges: 16 rows a 7 R <-> G edges + 16 columns a 7 R <-> G edges
    ASSERT_EQ(result.edges.size(), 2 * 16 * 7 * 2);
    ASSERT_EQ(result.edge_weights.size(), 2 * 16 * 7 * 2);

    // Edge weights should all be sqrt(2) * 255
    for (const SSInt weight: result.edge_weights) {
        EXPECT_EQ(weight, static_cast<SSInt>(std::sqrt(2) * 255));
    }
}
