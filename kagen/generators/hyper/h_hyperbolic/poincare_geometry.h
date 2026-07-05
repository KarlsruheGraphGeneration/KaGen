// kagen/generators/hyper/h_hyperbolic/poincare_geometry.h
#pragma once

#include "kagen/generators/hyper/h_hyperbolic/circular_interval.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace kagen::poincare_geometry {
template <typename Double>
struct Ball {
    Double cx;
    Double cy;
    Double radius;
    Double radius_sq;
};

template <typename Double>
struct AABB {
    Double min_x;
    Double max_x;
    Double min_y;
    Double max_y;
};
template <typename Double>
Ball<Double> MakeBall(const Double center_r, const Double center_phi, const Double hyperball_radius) {
    const Double s  = std::tanh(center_r / Double{2.0});
    const Double x0 = s * std::sin(center_phi);
    const Double y0 = s * std::cos(center_phi);

    const Double t     = std::tanh(hyperball_radius / Double{2.0});
    const Double denom = Double{1.0} - (s * s * t * t);

    const Double scale            = (Double{1.0} - (t * t)) / denom;
    const Double euclidean_radius = ((Double{1.0} - (s * s)) * t) / denom;

    const Double cx = scale * x0;
    const Double cy = scale * y0;

    return {
        .cx        = cx,
        .cy        = cy,
        .radius    = euclidean_radius,
        .radius_sq = euclidean_radius * euclidean_radius,
    };
}

template <typename Double>
bool AABBOutsideBall(const AABB<Double>& box, const Ball<Double>& ball) {
    const Double closest_x = std::clamp(ball.cx, box.min_x, box.max_x);
    const Double closest_y = std::clamp(ball.cy, box.min_y, box.max_y);

    const Double dx = ball.cx - closest_x;
    const Double dy = ball.cy - closest_y;

    return (dx * dx) + (dy * dy) > ball.radius_sq;
}

template <typename Double>
AABB<Double>
ComputeCellAABB(const Double min_r, const Double max_r, const Double min_phi, const Double max_phi, const Double eps) {
    const bool   full_circle = std::abs((max_phi - min_phi) - Double{2.0 * M_PI}) <= eps;
    const Double rho_min     = std::tanh(min_r / Double{2.0});
    const Double rho_max     = std::tanh(max_r / Double{2.0});

    std::array<Double, 6> phis{};
    int                   count = 0;

    phis[count++] = min_phi;
    phis[count++] = max_phi;

    const Double critical_angles[] = {
        Double{0.0},
        Double{M_PI / 2.0},
        Double{M_PI},
        Double{3.0 * M_PI / 2.0},
    };

    for (const Double angle: critical_angles) {
        if (full_circle || circular_interval::AngleInInterval(angle, min_phi, max_phi)) {
            phis[count++] = angle;
        }
    }

    Double min_x = std::numeric_limits<Double>::infinity();
    Double max_x = -std::numeric_limits<Double>::infinity();
    Double min_y = std::numeric_limits<Double>::infinity();
    Double max_y = -std::numeric_limits<Double>::infinity();

    for (int i = 0; i < count; ++i) {
        for (const Double rho: {rho_min, rho_max}) {
            const Double x = rho * std::sin(phis[i]);
            const Double y = rho * std::cos(phis[i]);

            min_x = std::min(min_x, x);
            max_x = std::max(max_x, x);
            min_y = std::min(min_y, y);
            max_y = std::max(max_y, y);
        }
    }

    return {.min_x = min_x, .max_x = max_x, .min_y = min_y, .max_y = max_y};
}

} // namespace kagen::poincare_geometry