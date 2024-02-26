#include "kagen/generators/image/kargb.h"

#include "kagen/kagen.h"

#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace kagen {
constexpr std::size_t kKargbIdentifierLength = 5;
constexpr std::size_t kKargbHeaderLength     = kKargbIdentifierLength + 2 * sizeof(std::uint64_t);

void CheckKARGB(const std::string& filename, bool& exists, bool& kargb) {
    exists = false;
    kargb  = false;

    std::ifstream in(filename, std::ios_base::binary);
    exists = !!in;
    if (!exists) {
        return;
    }

    std::array<char, kKargbIdentifierLength + 1> identifier;
    in.read(identifier.data(), kKargbIdentifierLength * sizeof(char));
    identifier[kKargbIdentifierLength] = 0;
    kargb                              = !std::strcmp(identifier.data(), "KARGB");
}

std::pair<SInt, SInt> ReadDimensions(const std::string& filename) {
    std::uint64_t                                rows;
    std::uint64_t                                cols;
    std::array<char, kKargbIdentifierLength + 1> identifier;

    std::ifstream in(filename, std::ios_base::binary);
    if (!in) {
        std::cerr << "Error: input file cannot  be read\n";
        std::exit(1);
    }

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

ImageRect::ImageRect(Pixels pixels, const SInt num_cols, const SInt overlap)
    : pixels_(std::move(pixels)),
      num_cols_(num_cols),
      overlap_(overlap) {}

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
} // namespace kagen
