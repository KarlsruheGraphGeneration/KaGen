// kagen/generators/hyper/h_hyperbolic/circular_interval.h
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>

namespace kagen::circular_interval {

template <typename Double>
struct SplitInterval {
    std::array<std::pair<Double, Double>, 2> parts{};
    int                                      count = 0;
};

template <typename Double>
Double NormalizePhi(Double phi) {
    const Double two_pi = Double{2.0 * M_PI};

    while (phi < Double{0.0}) {
        phi += two_pi;
    }

    while (phi >= two_pi) {
        phi -= two_pi;
    }

    return phi;
}

template <typename Double>
bool AngleInInterval(Double phi, Double min_phi, Double max_phi) {
    phi     = NormalizePhi(phi);
    min_phi = NormalizePhi(min_phi);
    max_phi = NormalizePhi(max_phi);

    if (min_phi <= max_phi) {
        return min_phi <= phi && phi <= max_phi;
    }

    return phi >= min_phi || phi <= max_phi;
}

template <typename Double>
Double AngularDistance(Double a, Double b) {
    a = NormalizePhi(a);
    b = NormalizePhi(b);

    const Double d = std::abs(a - b);
    return std::min(d, Double{2.0 * M_PI} - d);
}

template <typename Double>
Double CircularPhiWidth(Double min_phi, Double max_phi) {
    min_phi = NormalizePhi(min_phi);
    max_phi = NormalizePhi(max_phi);

    if (min_phi <= max_phi) {
        return max_phi - min_phi;
    }

    return (2.0 * M_PI - min_phi) + max_phi;
}

template <typename Double>
Double IntervalOverlap(Double a_begin, Double a_end, Double b_begin, Double b_end) {
    const Double begin = std::max(a_begin, b_begin);
    const Double end   = std::min(a_end, b_end);

    return std::max<Double>(0.0, end - begin);
}

template <typename Double>
SplitInterval<Double> Split(Double begin, Double end, Double eps) {
    const Double two_pi = Double{2.0 * M_PI};
    const Double width  = end - begin;

    SplitInterval<Double> result{};

    if (std::abs(std::abs(width) - two_pi) <= eps || std::abs(width) > two_pi) {
        result.parts[0] = {Double{0.0}, two_pi};
        result.count    = 1;
        return result;
    }

    begin = NormalizePhi(begin);
    end   = NormalizePhi(end);

    if (begin <= end) {
        result.parts[0] = {begin, end};
        result.count    = 1;
        return result;
    }

    result.parts[0] = {begin, two_pi};
    result.parts[1] = {Double{0.0}, end};
    result.count    = 2;
    return result;
}

template <typename Double>
Double CircularOverlap(Double a_begin, Double a_end, Double b_begin, Double b_end, Double eps) {
    Double overlap = Double{0.0};

    const auto a_parts = Split(a_begin, a_end, eps);
    const auto b_parts = Split(b_begin, b_end, eps);

    for (int i = 0; i < a_parts.count; ++i) {
        for (int j = 0; j < b_parts.count; ++j) {
            overlap += IntervalOverlap(
                a_parts.parts[i].first, a_parts.parts[i].second, b_parts.parts[j].first, b_parts.parts[j].second);
        }
    }

    return overlap;
}

} // namespace kagen::circular_interval