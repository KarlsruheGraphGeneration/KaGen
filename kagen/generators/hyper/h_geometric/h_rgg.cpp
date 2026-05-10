#include "kagen/generators/hyper/h_geometric/h_rgg.h"

#include "kagen/context.h"
#include "kagen/generators/hyper/h_geometric/hyper_rgg2d.h"
#include "kagen/kagen.h"

#include <memory>

namespace kagen {

PGeneratorConfig
HyperRGG2DFactory::NormalizeParameters(PGeneratorConfig config, PEID, const PEID size, const bool output) const {
    using namespace std::string_literals;

    EnsureSquarePowerOfTwoChunkSize(config, size, output);

    // TODO(clickup)[2026-05-08]: Only supports parameter combination n, r as of now

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

    return config;
}

std::unique_ptr<Generator>
HyperRGG2DFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<HyperRGG2D>(config, rank, size);
}

void HRGG::PushWeightIfRequested(
    [[maybe_unused]] const EdgeWeightConfig& config, [[maybe_unused]] const double& squared_distance,
    [[maybe_unused]] const double& squared_radius) { // TODO(clickup)[2026-05-10]: Not yet implemented
}

} // namespace kagen