#include "app/CLI11.h"

#include <array>
#include <fstream>
#include <iostream>
#include <string>

template <typename T, std::size_t length>
T read(std::ifstream& in) {
    if constexpr (std::is_same_v<T, std::string>) {
        std::string str(length, ' ');
        in.read(str.data(), length);
        return str;
    } else {
        std::array<std::uint8_t, length> data;
        in.read(reinterpret_cast<char*>(data.data()), length);
        T var = 0;
        for (std::size_t i = 0; i < length; ++i) {
            var <<= 8;
            var |= data[i];
        }
        return var;
    }
}

struct RGB {
    std::uint64_t             width;
    std::uint64_t             height;
    std::vector<std::uint8_t> R;
    std::vector<std::uint8_t> G;
    std::vector<std::uint8_t> B;
};

RGB parse_simple_psb(const std::string& filename, const bool quiet) {
    std::ifstream in(filename);
    in.seekg(0, std::ios_base::end);
    auto size = in.tellg();
    in.seekg(0);

    if (!quiet) {
        std::cout << "File size: " << size << std::endl;
    }

    auto signature    = read<std::string, 4>(in);
    auto version      = read<std::uint16_t, 2>(in);
    auto reserved     = read<std::uint64_t, 6>(in);
    auto num_channels = read<std::uint16_t, 2>(in);
    auto height       = read<std::uint32_t, 4>(in);
    auto width        = read<std::uint32_t, 4>(in);
    auto depth        = read<std::uint16_t, 2>(in);
    auto color_mode   = read<std::uint16_t, 2>(in);
    if (!quiet) {
        std::cout << "Signature: " << signature << std::endl;
        std::cout << "Version: " << version << std::endl;
        std::cout << "Reserved: " << reserved << std::endl;
        std::cout << "# channels: " << num_channels << std::endl;
        std::cout << "Height: " << height << std::endl;
        std::cout << "Width: " << width << std::endl;
        std::cout << "Depth: " << depth << std::endl;
        std::cout << "Color mode: " << color_mode << std::endl;
    }

    if (version != 2) {
        std::cerr << "Error: unexpected file format " << version
                  << " (only supporting 2 = PSB; add support for 1 = PSD when needed)\n";
        std::exit(1);
    }
    if (depth != 8) {
        std::cerr << "Error: unexpected channel depth " << depth << " (only supporting 8)\n";
        std::exit(1);
    }
    if (num_channels != 3) {
        std::cerr << "Error: unexpected number of channels " << num_channels << " (only supporting 3)\n";
        std::exit(1);
    }
    if (color_mode != 3) {
        std::cerr << "Error: unexpected color mode " << color_mode << " (only supporting 3 = RGB)\n";
        std::exit(1);
    }

    auto color_section_length = read<std::uint32_t, 4>(in);
    if (!quiet)
        std::cout << "Color Mode Data Length: " << color_section_length << " (ignoring)" << std::endl;
    in.seekg(color_section_length, std::ios_base::cur);

    auto image_resources_length = read<std::uint32_t, 4>(in);
    if (!quiet)
        std::cout << "Image Resources Length: " << image_resources_length << " (ignoring)" << std::endl;
    in.seekg(image_resources_length, std::ios_base::cur);

    auto layer_mask_length = read<std::uint64_t, 8>(in); // PSD = 4, PSB = 8
    if (!quiet)
        std::cout << "Layer and Mask Information Length: " << layer_mask_length << " (ignoring)" << std::endl;
    in.seekg(layer_mask_length, std::ios_base::cur);

    auto compression = read<std::uint16_t, 2>(in);
    if (!quiet) {
        std::cout << "Image Data Compression: " << compression << std::endl;
        std::cout << " -- current position in file: " << in.tellg() << std::endl;
        std::cout << " -- remaining bytes: " << size - in.tellg() << std::endl;
        std::cout << " -- div. by # channels: " << ((size - in.tellg()) % num_channels == 0 ? "yes" : "no")
                  << std::endl;
        std::cout << " -- entries per channel: " << (size - in.tellg()) / num_channels << std::endl;
        std::cout << " --  ... expected: " << static_cast<std::uint64_t>(width) * static_cast<std::uint64_t>(height)
                  << std::endl;
    }

    const std::uint64_t       num_pixels = static_cast<std::uint64_t>(width) * static_cast<std::uint64_t>(height);
    std::vector<std::uint8_t> red(num_pixels);
    std::vector<std::uint8_t> green(num_pixels);
    std::vector<std::uint8_t> blue(num_pixels);
    in.read(reinterpret_cast<char*>(red.data()), num_pixels);
    in.read(reinterpret_cast<char*>(green.data()), num_pixels);
    in.read(reinterpret_cast<char*>(blue.data()), num_pixels);
    if (!quiet) {
        std::cout << "Finished reading, position: " << in.tellg() << std::endl;
    }

    return {width, height, std::move(red), std::move(green), std::move(blue)};
}

enum class WeightModel {
    L2,
    IL2,
};

int main(int argc, char* argv[]) {
    std::string input_filename;
    std::string output_filename;
    std::string debug_output_filename;

    CLI::App app(
        "upsb2kargb: convert uncompressed PSB files to plain RGB files compatible with the Image Mesh Generator");
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

    RGB               rgb      = parse_simple_psb(input_filename, false);
    const std::size_t num_cols = rgb.width;
    const std::size_t num_rows = rgb.height;
    const auto&       R        = rgb.R;
    const auto&       G        = rgb.G;
    const auto&       B        = rgb.B;

    std::cout << "OK (" << num_cols << "x" << num_rows << ")\n";

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

    for (std::size_t row = 0; row < num_rows; ++row) {
        for (std::size_t col = 0; col < num_cols; ++col) {
            data[3 * col]     = R[row * num_cols + col];
            data[3 * col + 1] = G[row * num_cols + col];
            data[3 * col + 2] = B[row * num_cols + col];
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
