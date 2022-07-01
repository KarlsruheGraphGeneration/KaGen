#include "kagen/generators/geometric/delaunay.h"

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

    return config;
}
} // namespace

int Delaunay2DFactory::Requirements() const {
    return GeneratorRequirement::SQUARE_CHUNKS | GeneratorRequirement::POWER_OF_TWO_COMMUNICATOR_SIZE
           | GeneratorRequirement::ONE_CHUNK_PER_PE;
}

PGeneratorConfig Delaunay2DFactory::NormalizeParameters(PGeneratorConfig config, const bool output) const {
    return NormalizeParametersCommon(config, 3, output);
}

std::unique_ptr<Generator>
Delaunay2DFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<Delaunay2D>(config, rank, size);
}

int Delaunay3DFactory::Requirements() const {
    return GeneratorRequirement::CUBIC_CHUNKS | GeneratorRequirement::POWER_OF_TWO_COMMUNICATOR_SIZE;
}

PGeneratorConfig Delaunay3DFactory::NormalizeParameters(PGeneratorConfig config, const bool output) const {
    return NormalizeParametersCommon(config, 15.53 / 2.0, output);
}

std::unique_ptr<Generator>
Delaunay3DFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<Delaunay3D>(config, rank, size);
}
} // namespace kagen
