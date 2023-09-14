#pragma once

#include "kagen/generators/image/kargb.h"

#include <cmath>
#include <cstdint>

namespace kagen {
inline std::uint8_t Delta(const std::uint8_t lhs, const std::uint8_t rhs) {
    return std::max(lhs, rhs) - std::min(lhs, rhs);
}

inline double MinMaxRatio(const std::uint8_t lhs, const std::uint8_t rhs) {
    return lhs == rhs ? 1.0 : 1.0 * std::min(lhs, rhs) / std::max(lhs, rhs);
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

struct InvRatioWeightModel {
    double operator()(const RGB& lhs, const RGB& rhs) const {
        return 1.0 / (MinMaxRatio(lhs.r, rhs.r) * MinMaxRatio(lhs.g, rhs.g) * MinMaxRatio(lhs.b, rhs.b));
    }
};

struct RatioWeightModel {
    double operator()(const RGB& lhs, const RGB& rhs) const {
        return MinMaxRatio(lhs.r + 1, rhs.r + 1) * MinMaxRatio(lhs.g + 1, rhs.g + 1)
               * MinMaxRatio(lhs.b + 1, rhs.b + 1);
    }
};

struct SimilarityWeightModel {
    double operator()(const RGB& lhs, const RGB& rhs) const {
        constexpr static double sigma = 0.1;
        return std::exp(-l2_(lhs, rhs) / (2.0 * sigma * sigma));
    }

private:
    L2WeightModel l2_{};
};
} // namespace kagen
