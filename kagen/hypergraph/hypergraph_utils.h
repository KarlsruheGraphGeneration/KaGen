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
    const PGeneratorConfig& config, SInt identifier, double actual_lower_bound, double actual_upper_bound,
    RNGWrapper<>& rng);

bool RandomRadiusChecks(PGeneratorConfig& config);

PinRange
getRandomPinRange(SInt target_cell_size, SInt range_size, SInt target_cell_offset, SInt seed, Mersenne& mersenne);

double QuantileOrConstantHyperedgeRadius(const PGeneratorConfig& config);

template <typename RNG>
void FloydSampleGeometricAppend(
    const SInt universe_offset, const SInt universe_size, const SInt sample_size, RNG& rng, SInt& seed,
    std::vector<SInt>& out, std::unordered_set<SInt>& selected) {
    if (sample_size > universe_size) {
        throw ConfigurationError("Cannot sample more pins than available vertices");
    }
    selected.clear();

    if (selected.bucket_count() < sample_size * 2) {
        selected.reserve(static_cast<std::size_t>(sample_size) * 2);
    }

    if (sample_size == 0) {
        return;
    }

    if (sample_size <= 32) {
        std::array<SInt, 32> local;
        SInt                 count = 0;

        while (count < sample_size) {
            const long double x = std::min<long double>(
                static_cast<long double>(rng.GenerateCanonicalDoubleStream()), std::nextafter(1.0L, 0.0L));
            ++seed;

            const SInt candidate = universe_offset + static_cast<SInt>(x * static_cast<long double>(universe_size));

            bool duplicate = false;
            for (SInt i = 0; i < count; ++i) {
                if (local[i] == candidate) {
                    duplicate = true;
                    break;
                }
            }

            if (!duplicate) {
                local[count++] = candidate;
            }
        }

        out.insert(out.end(), local.begin(), local.begin() + count);
        return;
    }

    const SInt start = universe_size - sample_size;

    for (SInt j = start; j < universe_size; ++j) {
        const long double x = std::min<long double>(
            static_cast<long double>(rng.GenerateCanonicalDoubleStream()), std::nextafter(1.0L, 0.0L));
        ++seed;

        const SInt t = static_cast<SInt>(x * static_cast<long double>(j + 1));

        const SInt candidate = universe_offset + t;
        const SInt fallback  = universe_offset + j;

        selected.insert(selected.contains(candidate) ? fallback : candidate);
    }

    out.insert(out.end(), selected.begin(), selected.end());
}

enum class CellBallRelation : std::uint8_t { INSIDE, PARTIAL, OUTSIDE };
} // namespace kagen