
#include "kagen/hypergraph/hypergraph_utils.h"

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/sampling/hash.hpp"
#include "kagen/tools/rng_wrapper.h"

#include <vector>

namespace kagen {

double SampleHyperedgeRadius(const SInt identifier, const PGeneratorConfig& config) {
    const SInt seed = sampling::Spooky::hash(
        config.seed + (7919 * config.n) + identifier); // TODO(clickup)[2026-05-12]: Fix seed setup

    RNGWrapper<> rng = RNGWrapper(config);

    const double uniformRandom = rng.GenerateUniform<double>(seed);

    const double r_min = config.min_hyperedge_radius;
    const double r_max = config.max_hyperedge_radius;
    const double alpha = config.hyperedge_radius_exponent;

    const double transformed_lower_bound = std::pow(r_min, -alpha);
    const double transformed_upper_bound = std::pow(r_max, -alpha);

    return std::pow(
        transformed_lower_bound - (uniformRandom * (transformed_lower_bound - transformed_upper_bound)), -1.0 / alpha);
}

double getSampledOrConstantRadius(const PGeneratorConfig& config, const SInt identifier) {
    if (config.random_radius) {
        double sampled_radius = SampleHyperedgeRadius(identifier, config);
        return sampled_radius;
    }
    return config.r * config.r;
}

double getSampledOrConstantRadius(const PGeneratorConfig& config, const Vertex& center) {
    return getSampledOrConstantRadius(config, std::get<5>(center));
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

} // namespace kagen