#include "kagen/generators/geometric/delaunay.h"

#include "kagen/generators/geometric/delaunay/delaunay_2d.h"
#include "kagen/generators/geometric/delaunay/delaunay_3d.h"

namespace kagen {
namespace {
PGeneratorConfig NormalizeParametersCommon(PGeneratorConfig config, const double factor, const bool output) {
    if (config.n == 0) {
        if (config.m == 0) {
            throw ConfigurationError("at least one parameter out of {n, m} must be nonzero");
        }

        config.n = config.m / factor + 2;
        if (output) {
            std::cout << "Setting number of nodes to " << config.n << std::endl;
        }
    }

    if (config.streaming) {
        if (config.k < 1) {
            throw ConfigurationError("Number of chunks must be at least 1");
        }
        if (config.k > config.n) {
            throw ConfigurationError("Number of chunks must not exceed number of nodes");
        }
    }

    return config;
}
} // namespace
PGeneratorConfig
Delaunay2DFactory::NormalizeParameters(PGeneratorConfig config, PEID, const PEID size, const bool output) const {
    EnsureSquarePowerOfTwoChunkSize(config, size, output);
    // EnsurePowerOfTwoCommunicatorSize(config, size);

    // I think this is only needed in this case
    if (size > 1) EnsureOneChunkPerPE(config, size);

    return NormalizeParametersCommon(config, 3, output);
}

std::unique_ptr<Generator>
Delaunay2DFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<Delaunay2D>(config, rank, size);
}

PGeneratorConfig
Delaunay3DFactory::NormalizeParameters(PGeneratorConfig config, PEID, const PEID size, const bool output) const {
    EnsureCubicPowerOfTwoChunkSize(config, size, output);
    // EnsurePowerOfTwoCommunicatorSize(config, size);

    return NormalizeParametersCommon(config, 15.53 / 2.0, output);
}

std::unique_ptr<Generator>
Delaunay3DFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<Delaunay3D>(config, rank, size);
}
} // namespace kagen
