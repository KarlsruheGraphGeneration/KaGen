
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
    SInt identifier, const PGeneratorConfig& config, double lower_bound, double upper_bound, RNGWrapper<>& rng) {
    if (lower_bound <= 0.0 || upper_bound <= 0.0 || lower_bound > upper_bound) {
        throw ConfigurationError(
            std::format("Invalid hyperedge radius bounds: lower bound:{} upper bound:{}", lower_bound, upper_bound));
    }

    const SInt seed = sampling::Spooky::hash(config.seed + (7919 * config.n) + identifier);

    const double uniformRandom = rng.GenerateUniformDouble(seed);
    const double alpha         = config.hyperedge_radius_exponent;

    const double log_lower = -alpha * std::log(lower_bound);
    const double log_upper = -alpha * std::log(upper_bound);

    const double max_log = std::max(log_lower, log_upper);

    const double a = std::exp(log_lower - max_log);
    const double b = std::exp(log_upper - max_log);

    const double mixed = ((1.0 - uniformRandom) * a) + (uniformRandom * b);

    double sampled = std::exp(-(std::log(mixed) + max_log) / alpha);

    sampled = std::clamp(sampled, lower_bound, upper_bound);

    if (!std::isfinite(sampled) || sampled <= 0.0) {
        throw ConfigurationError(
            std::format(
                "Invalid sampled hyperedge radius: sampled={} "
                "u={} alpha={} lower={} upper={} "
                "lower^-alpha={} upper^-alpha={}",
                sampled, uniformRandom, alpha, lower_bound, upper_bound, a, b));
    }

    return sampled;
}

double getSampledOrConstantRadius(
    const PGeneratorConfig& config, const SInt identifier, const double actual_lower_bound,
    const double actual_upper_bound, RNGWrapper<>& rng) {
    if (!config.random_radius) {
        return config.r;
    }

    return SampleHyperedgeRadius(identifier, config, actual_lower_bound, actual_upper_bound, rng);
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

    const double lower = std::sqrt(2.0 / (M_PI * static_cast<double>(config.n)));
    const double upper = 1.0;
    const double q     = std::clamp(config.quantile, 0.0, 1.0);

    if (q <= 0.0) {
        return lower;
    }
    if (q >= 1.0) {
        return upper;
    }

    const double log_lower = -alpha * std::log(lower);
    const double log_upper = -alpha * std::log(upper); // 0 for upper=1

    const double max_log = std::max(log_lower, log_upper);

    const double a = std::exp(log_lower - max_log);
    const double b = std::exp(log_upper - max_log);

    const double mixed = ((1.0 - q) * a) + (q * b);

    return std::exp(-(std::log(mixed) + max_log) / alpha);
}

} // namespace kagen