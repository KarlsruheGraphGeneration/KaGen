
#include "kagen/hypergraph/hypergraph_utils.h"

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/generators/hyper/h_hyperbolic/hyper_hyperbolic.h"
#include "kagen/sampling/hash.hpp"
#include "kagen/tools/geometry.h"

#include <algorithm>

namespace kagen {

double SampleHyperedgeRadius(const SInt identifier, const PGeneratorConfig& config) {
    const SInt seed =
        sampling::Spooky::hash(config.seed + 7919 * config.n + identifier); // TODO(clickup)[2026-05-12]: Fix seed setup

    RNGWrapper<> rng = RNGWrapper(config);

    const double u = rng.GenerateUniform<double>(seed);

    const double r_min = config.min_hyperedge_radius;
    const double r_max = config.max_hyperedge_radius;
    const double alpha = config.hyperedge_radius_exponent;

    const double a = std::pow(r_min, -alpha);
    const double b = std::pow(r_max, -alpha);

    return std::pow(a - (u * (a - b)), -1.0 / alpha);
}

double getSampledOrConstantRadiusSquared(const PGeneratorConfig& config, const SInt identifier) {
    if (config.random_radius) {
        double sampled_radius = SampleHyperedgeRadius(identifier, config);
        return sampled_radius * sampled_radius;
    }
    return config.r * config.r;
}

double getSampledOrConstantRadiusSquared(const PGeneratorConfig& config, const Vertex& center) {
    return getSampledOrConstantRadiusSquared(config, std::get<5>(center));
}

} // namespace kagen