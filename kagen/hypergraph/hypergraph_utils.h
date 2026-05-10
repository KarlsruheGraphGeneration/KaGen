#pragma once

#include "kagen/kagen.h"

namespace kagen {

inline bool IsHypergraph(const Graph& graph) {
    return !graph.hyperedge_offsets.empty() || !graph.hyperedge_pins.empty();
}

inline SInt NumberOfLocalHyperedges(const Graph& graph) {
    return graph.hyperedge_offsets.empty() ? SInt{0} : static_cast<SInt>(graph.hyperedge_offsets.size() - 1);
}

inline SInt NumberOfLocalPins(const Graph& graph) {
    return static_cast<SInt>(graph.hyperedge_pins.size());
}
} // namespace kagen