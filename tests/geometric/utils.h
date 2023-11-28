#pragma once

#include "kagen/context.h"
#include "kagen/tools/geometry.h"

#include <cmath>
#include <numeric>
#include <utility>

namespace kagen::testing {
inline Edgelist CreateExpectedRGG2DEdges(PGeneratorConfig config, const Graph& graph) {
    std::vector<std::pair<SInt, SInt>> edges;
    for (SInt i = 0; i < config.n; i++) {
        auto [x1, y1] = graph.coordinates.first[i];
        for (SInt j = 0; j < config.n; j++) {
            auto [x2, y2] = graph.coordinates.first[j];
            // Comparing all coordinates
            if (i != j && std::hypot(x1 - x2, y1 - y2) <= config.r) {
                edges.emplace_back(i, j);
            }
        }
    }
    return edges;
}

inline Edgelist CreateExpectedRGG3DEdges(PGeneratorConfig config, const Graph& graph) {
    std::vector<std::pair<SInt, SInt>> edges;
    for (SInt i = 0; i < config.n; i++) {
        auto [x1, y1, z1] = graph.coordinates.second[i];
        for (SInt j = 0; j < config.n; j++) {
            auto [x2, y2, z2] = graph.coordinates.second[j];
            // Comparing all coordinates
            if (i != j && std::hypot(x1 - x2, y1 - y2, z1 - z2) < config.r) {
                edges.emplace_back(i, j);
            }
        }
    }
    return edges;
}

template <typename Double>
Edgelist CreateExpectedHyperbolicEdges(PGeneratorConfig config, const Graph& graph) {
    std::vector<std::pair<SInt, SInt>> edges;

    std::vector<std::pair<Double, Double>> polar_coordinates(config.n);
    for (SInt i = 0; i < config.n; ++i) {
        const auto& [x, y] = graph.coordinates.first[i];
        PGGeometry<Double>::CartesianToPolar({x, y}, polar_coordinates[i].first, polar_coordinates[i].second);
        polar_coordinates[i].second = PGGeometry<Double>::EuclideanRadiusToHyperbolic(polar_coordinates[i].second);
    }

    const Double alpha = (config.plexp - 1.0) / 2.0;
    const Double r     = PGGeometry<Double>::GetTargetRadius(config.n, config.n * config.avg_degree / 2, alpha);

    for (SInt i = 0; i < config.n; i++) {
        auto [phi1, r1] = polar_coordinates[i];

        for (SInt j = 0; j < config.n; j++) {
            auto [phi2, r2] = polar_coordinates[j];

            // Comparing all coordinates
            if (i != j && PGGeometry<Double>::HyperbolicDistance(r1, r2, phi1, phi2) <= r) {
                edges.emplace_back(i, j);
            }
        }
    }
    return edges;
}
} // namespace kagen::testing
