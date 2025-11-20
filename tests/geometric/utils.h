#pragma once

#include "kagen/context.h"
#include "kagen/tools/geometry.h"

#include <cmath>
#include <numeric>
#include <utility>

namespace kagen::testing {
inline Edgelist CreateExpectedRGG2DEdges(PGeneratorConfig config, const Graph& graph) {
    std::vector<std::pair<SInt, SInt>> edges;
    for (SInt i = 0; i < graph.NumberOfLocalVertices(); i++) {
        auto [x1, y1] = graph.coordinates.first[i];
        for (SInt j = 0; j < graph.NumberOfLocalVertices(); j++) {
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
    for (SInt i = 0; i < graph.NumberOfLocalVertices(); i++) {
        auto [x1, y1, z1] = graph.coordinates.second[i];
        for (SInt j = 0; j < graph.NumberOfLocalVertices(); j++) {
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
Double ComputeHyperbolicDistance(const Graph& graph, const SInt i, const SInt j) {
    std::pair<Double, Double> coord_i;
    std::pair<Double, Double> coord_j;

    const auto& [ix, iy] = graph.coordinates.first[i];
    PGGeometry<Double>::CartesianToPolar({ix, iy}, coord_i.first, coord_i.second);
    coord_i.second = PGGeometry<Double>::EuclideanRadiusToHyperbolic(coord_i.second);

    const auto& [jx, jy] = graph.coordinates.first[j];
    PGGeometry<Double>::CartesianToPolar({jx, jy}, coord_j.first, coord_j.second);
    coord_j.second = PGGeometry<Double>::EuclideanRadiusToHyperbolic(coord_j.second);

    auto [phi1, r1] = coord_i;
    auto [phi2, r2] = coord_j;
    return PGGeometry<Double>::HyperbolicDistance(r1, r2, phi1, phi2);
}

template <typename Double>
Edgelist CreateExpectedHyperbolicEdges(PGeneratorConfig config, const Graph& graph) {
    std::vector<std::pair<SInt, SInt>> edges;

    std::vector<std::tuple<Double, Double, char>> polar_coordinates(graph.NumberOfLocalVertices());
    for (SInt i = 0; i < graph.NumberOfLocalVertices(); ++i) {
        const auto& [x, y] = graph.coordinates.first[i];
        if (x >= 0 && y >= 0) {
            PGGeometry<Double>::CartesianToPolar({x, y}, std::get<0>(polar_coordinates[i]), std::get<1>(polar_coordinates[i]));
            std::get<1>(polar_coordinates[i]) = PGGeometry<Double>::EuclideanRadiusToHyperbolic(std::get<1>(polar_coordinates[i]));
            std::get<2>(polar_coordinates[i]) = 1;
        }
    }

    const Double alpha = (config.plexp - 1.0) / 2.0;
    const Double r     = PGGeometry<Double>::GetTargetRadius(config.n, config.n * config.avg_degree / 2, alpha);

    for (SInt i = 0; i < graph.NumberOfLocalVertices(); i++) {
        auto [phi1, r1, ok1] = polar_coordinates[i];
        if (!ok1) {
            continue;
        }

        for (SInt j = 0; j < graph.NumberOfLocalVertices(); j++) {
            auto [phi2, r2, ok2] = polar_coordinates[j];
            if (!ok2) {
                continue;
            }

            // Comparing all coordinates
            if (i != j && PGGeometry<Double>::HyperbolicDistance(r1, r2, phi1, phi2) <= r) {
                edges.emplace_back(i, j);
            }
        }
    }
    return edges;
}
} // namespace kagen::testing
