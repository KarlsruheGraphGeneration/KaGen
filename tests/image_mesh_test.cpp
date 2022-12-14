#include <gtest/gtest.h>

#include "kagen/generators/image/kargb.h"

using namespace kagen;

TEST(KARGBTest, reads_checkerboard_dimension) {
    const auto [nrows, ncols] = ReadDimensions("data/rgcheckerboard.kargb");
    EXPECT_EQ(nrows, 16);
    EXPECT_EQ(ncols, 16);
}

TEST(KARGBTest, reads_full_checkerboard) {
    const SInt nrows = 16;
    const SInt ncols = 16;

    // G R G R G R G R
    // R G R G R G R G
    // ...
    const auto pixels = ReadRect("data/rgcheckerboard.kargb", 0, 0, nrows, ncols, 0);
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

TEST(KARGBTest, reads_inner_checkerboard) {
    const SInt from_row = 2;
    const SInt from_col = 2;
    const SInt to_row   = 14;
    const SInt to_col   = 14;

    // G R G R G R
    // R G R G R G
    // ...
    const auto pixels = ReadRect("data/rgcheckerboard.kargb", from_row, from_col, to_row, to_col, 0);

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

TEST(KARGBTest, reads_full_checkerboard_through_overlap) {
    const SSInt from_row = 2;
    const SSInt from_col = 2;
    const SSInt to_row   = 14;
    const SSInt to_col   = 14;

    // G R G R G R G R
    // R G R G R G R G
    // ...
    const auto pixels = ReadRect("data/rgcheckerboard.kargb", from_row, from_col, to_row, to_col, 2);
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

TEST(KARGBTest, reads_checkerboard_with_virtual_overlap) {
    const SSInt from_row = 0;
    const SSInt from_col = 0;
    const SSInt to_row   = 2;
    const SSInt to_col   = 2;

    // G R G R G R G R
    // R G R G R G R G
    // ...
    const auto pixels = ReadRect("data/rgcheckerboard.kargb", from_row, from_col, to_row, to_col, 2);
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
