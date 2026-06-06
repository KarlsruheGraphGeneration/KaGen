#pragma once

#include "kagen/kagen.h"

namespace kagen {

// n, x_off, y_off, generated, offset
using Chunk = std::tuple<SInt, LPFloat, LPFloat, bool, SInt>;
// n, x_off, y_off, generated, offset
using Cell = std::tuple<SInt, LPFloat, LPFloat, bool, SInt>;
// x, y, id
using Vertex = std::tuple<LPFloat, LPFloat, SInt>;

struct Point {
    double x;
    double y;
};

inline double Cross(const Point& a, const Point& b) {
    return a.x * b.y - a.y * b.x;
}

inline double PolygonArea(const std::vector<Point>& poly) {
    if (poly.size() < 3) {
        return 0.0;
    }

    double area = 0.0;

    for (size_t i = 0; i < poly.size(); ++i) {
        const Point& a = poly[i];
        const Point& b = poly[(i + 1) % poly.size()];
        area += Cross(a, b);
    }

    return 0.5 * std::abs(area);
}

} // namespace kagen