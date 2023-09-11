#include "app/CLI11.h"

#include <Magick++.h>

#include <array>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    Magick::InitializeMagick(*argv);

    std::string input_filename;
    std::string output_filename;
    std::string debug_output_filename;

    CLI::App app("img2kargb: convert images to plain RGB files compatible with the Image Mesh Generator");
    app.add_option("input filename", input_filename)->check(CLI::ExistingFile)->required();
    app.add_option("-o,--output", output_filename);
    app.add_option("-d,--debug-output", debug_output_filename);
    CLI11_PARSE(app, argc, argv);

    if (output_filename.empty()) {
        output_filename = input_filename + ".kargb";
    }

    const bool write_debug_tga = !debug_output_filename.empty();

    // Load input image into memory
    std::cout << "Reading " << input_filename << " ... " << std::flush;

    Magick::Image img(input_filename);
    img.colorSpaceType(Magick::sRGBColorspace); // ?
    img.type(Magick::TrueColorType);            // ?

    const std::uint64_t num_rows = img.rows();
    const std::uint64_t num_cols = img.columns();

    std::cout << "OK (" << num_cols << "x" << num_rows << ", " << img.magick() << ", " << img.channels()
              << " channels)\n";

    if (write_debug_tga) {
        std::array<std::uint8_t, 18> tga;
        std::fill(tga.begin(), tga.end(), 0);
        tga[2]  = 2; // Uncompressed true-color image
        tga[12] = num_cols & 255;
        tga[13] = num_cols >> 8;
        tga[14] = num_rows & 255;
        tga[15] = num_rows >> 8;
        tga[16] = 24; // 3 * 8 bits per pixel
        tga[17] = 32; // Top-to-bottom, left-to-right

        std::ofstream debug_out(debug_output_filename, std::ios_base::binary | std::ios_base::trunc);
        debug_out.write(reinterpret_cast<const char*>(tga.data()), sizeof(std::uint8_t) * tga.size());
    }

    std::cout << "Writing RGB channels to " << output_filename << " ... " << std::flush;

    std::ofstream out(output_filename, std::ios_base::binary | std::ios_base::trunc);
    out.write("KARGB", 5 * sizeof(char));
    out.write(reinterpret_cast<const char*>(&num_rows), sizeof(std::uint64_t));
    out.write(reinterpret_cast<const char*>(&num_cols), sizeof(std::uint64_t));

    std::vector<std::uint8_t> data(num_cols * 3);

    const auto* pixels = img.getConstPixels(0, 0, num_cols, num_rows);

    for (std::size_t row = 0; row < num_rows; ++row) {
        for (std::size_t col = 0; col < num_cols; ++col) {
            data[3 * col]     = pixels[0];
            data[3 * col + 1] = pixels[1];
            data[3 * col + 2] = pixels[2];
            pixels += img.channels();
        }

        out.write(reinterpret_cast<const char*>(data.data()), sizeof(std::uint8_t) * data.size());
        if (write_debug_tga) {
            // Swap R and B channels for TGA debug output
            for (std::size_t col = 0; col < num_cols; ++col) {
                std::swap(data[3 * col], data[3 * col + 2]);
            }
            std::ofstream debug_out(debug_output_filename, std::ios_base::binary | std::ios_base::app);
            debug_out.write(reinterpret_cast<const char*>(data.data()), sizeof(std::uint8_t) * data.size());
        }
    }

    out.close();
    std::cout << "OK" << std::endl;
}
