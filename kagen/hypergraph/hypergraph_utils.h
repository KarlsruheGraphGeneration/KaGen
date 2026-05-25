#pragma once

#include "kagen/context.h"
#include "kagen/generators/geometric/geometric_util.h"
#include "kagen/kagen.h"

#include <tuple>
#include <vector>

namespace kagen {
inline bool IsHypergraph(const Graph& graph) {
    return !graph.hyperedge_offsets.empty() || !graph.hyperedge_pins.empty();
}

inline SInt NumberOfLocalHyperedges(Graph& graph) {
    if (!graph.hyperedge_offsets.empty()) {
        return static_cast<SInt>(graph.hyperedge_offsets.size() - 1);
    }

    if (!graph.hyperedge_range_offsets.empty()) {
        return static_cast<SInt>(graph.hyperedge_range_offsets.size() - 1);
    }

    return 0;
}

inline SInt NumberOfLocalPins(const Graph& graph) {
    return static_cast<SInt>(graph.hyperedge_pins.size());
}
/**
 * Checks ```random_radius`` flag of @link PGeneratorConfig
 */
double getSampledOrConstantRadius(const PGeneratorConfig& config, SInt identifier, double lower_bound = 0.0);

bool RandomRadiusChecks(PGeneratorConfig& config);

PinRange getRandomPinRange(
    SInt target_cell_size, SInt range_size, SInt target_cell_offset, SInt seed, const PGeneratorConfig& config);

enum class CellBallRelation : std::uint8_t { INSIDE, PARTIAL, OUTSIDE };
} // namespace kagen