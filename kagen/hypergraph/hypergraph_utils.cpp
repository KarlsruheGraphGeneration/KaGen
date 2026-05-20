
#include "kagen/hypergraph/hypergraph_utils.h"

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/sampling/hash.hpp"
#include "kagen/tools/rng_wrapper.h"

#include <vector>

namespace kagen {

double SampleHyperedgeRadius(const SInt identifier, const PGeneratorConfig& config) {
    const SInt seed = sampling::Spooky::hash(config.seed + (7919 * config.n) + identifier);

    RNGWrapper<> rng = RNGWrapper(config);

    const double     uniformRandom     = rng.GenerateUniform<double>(seed);
    constexpr double expected_vertices = 32.0;
    double           lower_bound       = 1e-5;
    if (config.generator == GeneratorType::H_RHG) {
        lower_bound = std::sqrt(expected_vertices / static_cast<double>(config.n));
    } else if (config.generator == GeneratorType::HRGG_2D) {
        lower_bound = std::sqrt(expected_vertices / (M_PI * static_cast<double>(config.n)));
    }

    const double upper_bound = 1.0;
    const double alpha       = config.hyperedge_radius_exponent;

    const double transformed_lower_bound = std::pow(lower_bound, -alpha);
    const double transformed_upper_bound = std::pow(upper_bound, -alpha);

    const double sampled = std::pow(
        transformed_lower_bound - uniformRandom * (transformed_lower_bound - transformed_upper_bound), -1.0 / alpha);

    if (!std::isfinite(sampled) || sampled <= 0.0 || sampled > 1.0) {
        throw ConfigurationError("Invalid sampled hyperedge radius");
    }

    return sampled;
}

double getSampledOrConstantRadius(const PGeneratorConfig& config, const SInt identifier) {
    if (config.random_radius) {
        double sampled_radius = SampleHyperedgeRadius(identifier, config);
        return sampled_radius;
    }
    return config.r;
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