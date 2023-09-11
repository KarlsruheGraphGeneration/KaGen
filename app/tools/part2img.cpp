#include "app/CLI11.h"

#include <Magick++.h>

#include <array>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    // @todo: support output formats other than TGA
    // @todo: support other methods to visualize the (pseudo-)partition file
    // @todo: support grid parameters

    Magick::InitializeMagick(*argv);

    std::string input_image_filename;
    std::string input_partition_filename;
    std::string output_image_filename;

    CLI::App app(
        "part2img: visualize the partition result of an image mesh graph by marking pixels in the original image");
    app.add_option("input image", input_image_filename)->check(CLI::ExistingFile)->required();
    app.add_option("input partition", input_partition_filename)->check(CLI::ExistingFile)->required();
    app.add_option("-o,--output", output_image_filename);
    CLI11_PARSE(app, argc, argv);

    if (output_image_filename.empty()) {
        output_image_filename = input_image_filename + ".tga";
    }

    // Load input image into memory
    std::cout << "Reading " << input_image_filename << " ... " << std::flush;

    Magick::Image img(input_image_filename);
    img.colorSpaceType(Magick::sRGBColorspace); // ?
    img.type(Magick::TrueColorType);            // ?

    const std::uint64_t num_rows = img.rows();
    const std::uint64_t num_cols = img.columns();

    std::cout << "OK (" << num_cols << "x" << num_rows << ", " << img.magick() << ", " << img.channels()
              << " channels)\n";

    std::cout << "Reading " << input_partition_filename << " ... " << std::flush;
    std::vector<std::uint64_t> partition;
    partition.reserve(img.rows() * img.columns());
    {
        std::ifstream partition_in(input_partition_filename);
        std::uint64_t block;
        while (partition_in >> block) {
            partition.push_back(block);
        }
        std::cout << "OK (read " << partition.size() << " IDs)" << std::endl;
        if (partition.size() != img.rows() * img.columns()) {
            std::cerr << "Error: expected " << img.rows() * img.columns() << " IDs\n";
            std::exit(1);
        }
    }

    std::array<std::uint8_t, 18> tga;
    std::fill(tga.begin(), tga.end(), 0);
    tga[2]  = 2; // Uncompressed true-color image
    tga[12] = num_cols & 255;
    tga[13] = num_cols >> 8;
    tga[14] = num_rows & 255;
    tga[15] = num_rows >> 8;
    tga[16] = 24; // 3 * 8 bits per pixel
    tga[17] = 32; // Top-to-bottom, left-to-right

    std::cout << "Writing output to " << output_image_filename << " ... " << std::flush;
    std::ofstream out(output_image_filename, std::ios_base::binary | std::ios_base::trunc);
    out.write(reinterpret_cast<const char*>(tga.data()), sizeof(std::uint8_t) * tga.size());

    const auto*               pixels = img.getConstPixels(0, 0, num_cols, num_rows);
    std::vector<std::uint8_t> data(num_cols * 3);

    for (std::size_t row = 0; row < num_rows; ++row) {
        for (std::size_t col = 0; col < num_cols; ++col) {
            data[3 * col]     = pixels[2];
            data[3 * col + 1] = pixels[1];
            data[3 * col + 2] = pixels[0];
            pixels += img.channels();

            const std::uint64_t my_block           = partition[row * num_cols + col];
            bool                is_boundary_vertex = false;
            if (col + 1 < num_cols) {
                is_boundary_vertex |= my_block != partition[row * num_cols + col + 1];
            }
            if (row > 0) {
                is_boundary_vertex |= my_block != partition[(row - 1) * num_cols + col];
            }
            if (col > 0) {
                is_boundary_vertex |= my_block != partition[row * num_cols + col - 1];
            }
            if (row + 1 < num_rows) {
                is_boundary_vertex |= my_block != partition[(row + 1) * num_cols + col];
            }

            if (is_boundary_vertex) {
                data[3 * col]     = 0;
                data[3 * col + 1] = 0;
                data[3 * col + 2] = 255;
            }
        }

        out.write(reinterpret_cast<const char*>(data.data()), sizeof(std::uint8_t) * data.size());
    }

    std::cout << "OK" << std::endl;
}
