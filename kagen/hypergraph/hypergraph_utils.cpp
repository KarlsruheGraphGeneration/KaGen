
#include "kagen/hypergraph/hypergraph_utils.h"

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"
#include "kagen/tools/rng_wrapper.h"

#include <algorithm>
#include <cstdlib>
#include <format>

namespace kagen {

double SampleHyperedgeRadius(
    const SInt identifier, const PGeneratorConfig& config, const double lower_bound, const double upper_bound) {
    if (lower_bound <= 0.0 || upper_bound <= 0.0 || lower_bound > upper_bound) {
        throw ConfigurationError(
            std::format("Invalid hyperedge radius bounds: lower bound:{} upper bound:{}", lower_bound, upper_bound));
    }

    const SInt seed = sampling::Spooky::hash(config.seed + (7919 * config.n) + identifier);

    RNGWrapper<> rng(config);

    const double uniformRandom     = rng.GenerateUniform<double>(seed);
    const double alpha             = config.hyperedge_radius_exponent;
    const double transformed_lower = std::pow(lower_bound, -alpha);
    const double transformed_upper = std::pow(upper_bound, -alpha);

    double sampled =
        std::pow(transformed_lower - (uniformRandom * (transformed_lower - transformed_upper)), -1.0 / alpha);

    sampled = std::min(sampled, upper_bound);
    sampled = std::max(sampled, lower_bound);

    if (!std::isfinite(sampled) || sampled <= 0.0) {
        throw ConfigurationError("Invalid sampled hyperedge radius");
    }

    return sampled;
}

double getSampledOrConstantRadius(
    const PGeneratorConfig& config, const SInt identifier, const double actual_lower_bound,
    const double actual_upper_bound) {
    if (!config.random_radius) {
        return config.r;
    }

    return SampleHyperedgeRadius(identifier, config, actual_lower_bound, actual_upper_bound);
}

/**
 * Performs checks to make sure that config is valid in terms of radius data.
 */
bool RandomRadiusChecks(PGeneratorConfig& config) {
    if (config.random_radius) {
        if (config.min_hyperedge_radius <= 0.0 || config.max_hyperedge_radius <= 0.0
            || config.hyperedge_radius_exponent <= 0.0) {
            throw ConfigurationError("non-constant hyperball radii require min radius, max radius and exponent > 0");
        }
        if (config.min_hyperedge_radius > config.max_hyperedge_radius) {
            throw ConfigurationError("min hyperedge radius must not exceed max hyperedge radius");
        }
    } else if (config.r <= 0) {
        throw ConfigurationError("hyperbolic hypergraph generator requires hyperedge radius > 0");
    }
    return true;
}

PinRange
getRandomPinRange(SInt target_cell_size, SInt range_size, SInt target_cell_offset, SInt seed, Mersenne& mersenne) {
    const SInt max_start = target_cell_size - range_size;

    const LPFloat u              = mersenne.Random();
    const SInt    interval_start = static_cast<SInt>(std::floor(u * static_cast<LPFloat>(max_start)));

    return {.begin = target_cell_offset + interval_start, .end = target_cell_offset + interval_start + range_size};
}

double QuantileOrConstantHyperedgeRadius(const PGeneratorConfig& config) {
    if (!config.random_radius) {
        return config.r;
    }

    if (config.n <= 0) {
        throw ConfigurationError("Quantile radius requires n > 0");
    }

    const double alpha = config.hyperedge_radius_exponent;

    if (alpha <= 0.0 || !std::isfinite(alpha)) {
        throw ConfigurationError("Quantile radius requires exponent > 0");
    }

    const double lower_bound = std::sqrt(2.0 / (M_PI * static_cast<double>(config.n)));

    const double upper_bound = 1.0;

    const double clamped_q = std::clamp(config.quantile, 0.0, 1.0);

    const double transformed_lower = std::pow(lower_bound, -alpha);
    const double transformed_upper = std::pow(upper_bound, -alpha);

    return std::pow(transformed_lower - clamped_q * (transformed_lower - transformed_upper), -1.0 / alpha);
}

} // namespace kagen