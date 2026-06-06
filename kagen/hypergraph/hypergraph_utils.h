#pragma once

#include "kagen/context.h"
#include "kagen/kagen.h"
#include "kagen/tools/mersenne.h"
#include "kagen/tools/rng_wrapper.h"

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
double getSampledOrConstantRadius(
    const PGeneratorConfig& config, SInt identifier, double actual_lower_bound, double actual_upper_bound);

bool RandomRadiusChecks(PGeneratorConfig& config);

PinRange
getRandomPinRange(SInt target_cell_size, SInt range_size, SInt target_cell_offset, SInt seed, Mersenne& mersenne);

double QuantileOrConstantHyperedgeRadius(const PGeneratorConfig& config);

enum class CellBallRelation : std::uint8_t { INSIDE, PARTIAL, OUTSIDE };
} // namespace kagen