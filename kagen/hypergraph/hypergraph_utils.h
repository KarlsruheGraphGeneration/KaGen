#pragma once

#include "kagen/kagen.h"

namespace kagen {
// phi, r, x, y, gamma, id
using Vertex = std::tuple<double, double, double, double, double, SInt>;
inline bool IsHypergraph(const Graph& graph) {
    return !graph.hyperedge_offsets.empty() || !graph.hyperedge_pins.empty();
}

inline SInt NumberOfLocalHyperedges(const Graph& graph) {
    return graph.hyperedge_offsets.empty() ? SInt{0} : static_cast<SInt>(graph.hyperedge_offsets.size() - 1);
}

inline SInt NumberOfLocalPins(const Graph& graph) {
    return static_cast<SInt>(graph.hyperedge_pins.size());
}
/**
* Checks ```random_radius`` flag of @link PGeneratorConfig
*/
double getSampledOrConstantRadiusSquared(const PGeneratorConfig& config, const SInt identifier);

double getSampledOrConstantRadiusSquared(const PGeneratorConfig& config, const Vertex& center);
} // namespace kagen