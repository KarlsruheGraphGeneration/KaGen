#pragma once

#include "kagen/kagen.h"

#include <cstdint>
#include <string>
#include <vector>

namespace kagen {
struct RGB {
    RGB() = default;
    RGB(const std::uint8_t r, const std::uint8_t g, const std::uint8_t b) : r(r), g(g), b(b) {}

    std::uint8_t r;
    std::uint8_t g;
    std::uint8_t b;
};

using Pixels = std::vector<RGB>;

class ImageRect {
public:
    ImageRect(Pixels pixels, SInt num_cols, SInt overlap);

    inline const RGB& GetPixel(const SSInt row, const SSInt col) const {
        return pixels_[(row + overlap_) * (num_cols_ + 2 * overlap_) + (col + overlap_)];
    }

private:
    Pixels pixels_;
    SInt   num_cols_;
    SInt   overlap_;
};

void CheckKARGB(const std::string& filename, bool& exists, bool& kargb);

std::pair<SInt, SInt> ReadDimensions(const std::string& filename);

ImageRect ReadRect(const std::string& filename, SInt from_row, SInt from_col, SInt to_row, SInt to_col, SInt overlap);
} // namespace kagen
