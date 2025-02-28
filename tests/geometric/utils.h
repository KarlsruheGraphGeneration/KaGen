#pragma once

#include "kagen/context.h"
#include "kagen/tools/geometry.h"

#include <cmath>
#include <utility>

namespace kagen::testing {

inline Edgelist GenerateLocalRGG2DEdges(const Graph& graph, const double radius) {
    std::vector<std::pair<SInt, SInt>> edges;

    for (SInt u = graph.vertex_range.first; u < graph.vertex_range.second; ++u) {
        const auto [xu, yu] = graph.coordinates.first[u];

        for (SInt v = u + 1; v < graph.vertex_range.second; ++v) {
            const auto [xv, yv] = graph.coordinates.first[v];

            if (std::hypot(xu - xv, yu - yv) <= radius) {
                edges.emplace_back(u, v);
                edges.emplace_back(v, u);
            }
        }
    }

    std::sort(edges.begin(), edges.end());

    return edges;
}

inline Edgelist GenerateLocalRGG3DEdges(const Graph& graph, const double radius) {
    std::vector<std::pair<SInt, SInt>> edges;

    for (SInt u = graph.vertex_range.first; u < graph.vertex_range.second; ++u) {
        const auto [xu, yu, zu] = graph.coordinates.second[u];

        for (SInt v = u + 1; v < graph.vertex_range.second; ++v) {
            const auto [xv, yv, zv] = graph.coordinates.second[v];

            if (std::hypot(xu - xv, yu - yv, zu - zv) <= radius) {
                edges.emplace_back(u, v);
                edges.emplace_back(v, u);
            }
        }
    }

    std::sort(edges.begin(), edges.end());

    return edges;
}

template <typename Double>
Edgelist
GenerateLocalRHGEdges(PGeneratorConfig config, const Graph& graph, const double plexp, const double avg_degree) {
    std::vector<std::pair<SInt, SInt>> edges;

    std::vector<std::pair<Double, Double>> coordinates(graph.NumberOfLocalVertices());
    for (SInt u = graph.vertex_range.first; u < graph.vertex_range.second; ++u) {
        const auto& [x, y] = graph.coordinates.first[u];
        PGGeometry<Double>::CartesianToPolar({x, y}, coordinates[u].first, coordinates[u].second);
        coordinates[u].second = PGGeometry<Double>::EuclideanRadiusToHyperbolic(coordinates[u].second);
    }

    const Double radius = PGGeometry<Double>::GetTargetRadius(
        graph.NumberOfLocalVertices(), graph.NumberOfLocalVertices() * avg_degree / 2, (plexp - 1.0) / 2.0);

    for (SInt u = graph.vertex_range.first; u < graph.vertex_range.second; u++) {
        auto [phiu, ru] = coordinates[u];

        for (SInt v = u + 1; v < graph.vertex_range.second; v++) {
            auto [phiv, rv] = coordinates[v];

            if (PGGeometry<Double>::HyperbolicDistance(ru, rv, phiu, phiv) <= radius) {
                edges.emplace_back(u, v);
                edges.emplace_back(v, u);
            }
        }
    }

    std::sort(edges.begin(), edges.end());

    return edges;
}

} // namespace kagen::testing
