#include "kagen/generators/hyper/h_geometric/h_rgg.h"

#include "kagen/context.h"
#include "kagen/generators/hyper/h_geometric/hyper_rgg2d.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"

#include <memory>

namespace kagen {

PGeneratorConfig
HyperRGG2DFactory::NormalizeParameters(PGeneratorConfig config, PEID rank, const PEID size, const bool output) const {
    using namespace std::string_literals;

    EnsureSquarePowerOfTwoChunkSize(config, size, output);

    // TODO(clickup)[2026-05-08]: Only supports parameter combination n, r as of now

    if (config.random_radius && config.r == 0) {
        const double expected_vertices = 32.0;

        config.r = std::sqrt(expected_vertices / (M_PI * static_cast<double>(config.n)));
    }

    if (config.k > 1.0 / (config.r * config.r)) {
        throw ConfigurationError(
            "number of chunks (="s + std::to_string(config.k)
            + ") must be smaller than 1/radius^2 (=" + std::to_string(1.0 / (config.r * config.r)) + ")");
    }

    config.is_hypergraph         = true;
    config.validate_simple_graph = false;

    if (config.edge_weights.generator_type != EdgeWeightGeneratorType::DEFAULT) {
        throw ConfigurationError("edge weights are not implemented for hypergraphs yet");
    }

    //
    // Convert expected hyperedge-size bounds to Euclidean radius bounds.
    //
    // For uniformly distributed vertices in the unit square:
    //
    //     E[|e|] = n * pi * r^2
    //
    // so
    //
    //     r = sqrt(E[|e|] / (n * pi)).
    //
    if (config.size_dist_lower_bound > 0 && config.size_dist_upper_bound > 0
        && config.size_dist_lower_bound > config.size_dist_upper_bound) {
        throw ConfigurationError(
            "lower hyperedge size bound "
            "must not exceed upper hyperedge size bound");
    }

    config.min_hyperedge_radius = EuclideanRadiusForExpectedHyperedgeSize(config.size_dist_lower_bound, config.n);

    if (config.size_dist_upper_bound > 0) {
        config.max_hyperedge_radius = EuclideanRadiusForExpectedHyperedgeSize(config.size_dist_upper_bound, config.n);
    }

    if (config.min_hyperedge_radius != -1.0 && config.max_hyperedge_radius != -1.0
        && config.min_hyperedge_radius > config.max_hyperedge_radius) {
        throw ConfigurationError(
            "lower hyperedge size/radius bound "
            "must not exceed upper bound");
    }

    //
    // If a pin budget is given, determine the radius-distribution
    // exponent that yields the requested expected number of pins.
    //
    if (config.size_dist_pin_budget > 0) {
        if (config.n <= 0 || config.m <= 0) {
            throw ConfigurationError("expected total pins requires n > 0 and m > 0");
        }

        if (!config.random_radius) {
            throw ConfigurationError("expected total pins requires random_radius=true");
        }

        double lower_bound = EuclideanRadiusForExpectedHyperedgeSize(2, config.n);

        double upper_bound = 1.0;

        if (config.min_hyperedge_radius != -1.0) {
            lower_bound = config.min_hyperedge_radius;
        }

        if (config.max_hyperedge_radius != -1.0) {
            upper_bound = config.max_hyperedge_radius;
        }

        const double target_expected_r2 = static_cast<double>(config.size_dist_pin_budget)
                                          / (static_cast<double>(config.m) * static_cast<double>(config.n) * M_PI);

        config.hyperedge_radius_exponent =
            SolveRadiusExponentForExpectedPins(target_expected_r2, lower_bound, upper_bound);

        if (rank == 0) {
            std::cout << " Chosen radius exponent = " << config.hyperedge_radius_exponent << '\n';
        }
    }

    return config;
}

std::unique_ptr<Generator>
HyperRGG2DFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<HyperRGG2D>(config, rank, size);
}

void HRGG::PushWeightIfRequested(
    [[maybe_unused]] const EdgeWeightConfig& config, [[maybe_unused]] const double& squared_distance,
    [[maybe_unused]] const double& squared_radius) {
    // TODO(clickup)[2026-05-10]: Not yet implemented
}

} // namespace kagen