#include "kagen/generators/generator.h"
#include "kagen/context.h"
#include "kagen/definitions.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>

namespace kagen {
Generator::~Generator() = default;

void Generator::Generate() {
    edges_.clear();
    coordinates_.first.clear();
    coordinates_.second.clear();
    GenerateImpl();
}

void Generator::Finalize(MPI_Comm) {}

const EdgeList& Generator::GetEdges() const {
    return edges_;
}

Graph Generator::TakeResult() {
    return {
        std::move(edges_), vertex_range_, std::move(coordinates_), std::move(vertex_weights_),
        std::move(edge_weights_)};
}

void Generator::SetVertexRange(const VertexRange vertex_range) {
    vertex_range_ = vertex_range;
}

void Generator::FilterDuplicateEdges() {
    std::sort(edges_.begin(), edges_.end());
    auto it = std::unique(edges_.begin(), edges_.end());
    edges_.erase(it, edges_.end());
}

GeneratorFactory::~GeneratorFactory() = default;

PGeneratorConfig GeneratorFactory::NormalizeParameters(PGeneratorConfig config, const PEID size, bool) const {
    if (config.k == 0) {
        config.k = static_cast<SInt>(size);
    }
    return config;
}

namespace {
bool IsPowerOfTwo(const SInt value) {
    return (value & (value - 1)) == 0;
}

bool IsSquare(const SInt value) {
    const SInt root = std::round(std::sqrt(value));
    return root * root == value;
}

bool IsCubic(const SInt value) {
    const SInt root = std::round(std::cbrt(value));
    return root * root * root == value;
}

SInt FindSquareMultipleOf(const SInt value) {
    if (IsSquare(value)) {
        return value;
    }
    if (IsPowerOfTwo(value)) {
        return 2 * value; // every 2nd power of two is square
    }

    // find smallest square number that is a multiple of value
    const SInt root = std::sqrt(value);
    for (SInt cur = root; cur < value; ++cur) {
        const SInt squared = cur * cur;
        if (squared % value == 0) {
            return squared;
        }
    }
    return value * value;
}

SInt FindCubeMultipleOf(const SInt value) {
    if (IsCubic(value)) {
        return value;
    }
    if (IsPowerOfTwo(value)) {
        return IsCubic(value * 2) ? value * 2 : value * 4;
    }

    // find smallest cubic number that is a multiple of value
    const SInt root = std::cbrt(value);
    for (SInt cur = root; cur < value; ++cur) {
        const SInt cubed = cur * cur * cur;
        if (cubed % value == 0) {
            return cubed;
        }
    }
    return value * value * value;
}
} // namespace

void GeneratorFactory::EnsurePowerOfTwoCommunicatorSize(PGeneratorConfig&, const PEID size) const {
    if (!IsPowerOfTwo(size)) {
        throw ConfigurationError("number of PEs must be a power of two");
    }
}

void GeneratorFactory::EnsureSquareChunkSize(PGeneratorConfig& config, const PEID size) const {
    if (config.k == 0) {
        config.k = FindSquareMultipleOf(size);
    } else if (!IsSquare(config.k)) {
        throw ConfigurationError("number of chunks must be square");
    }
}

void GeneratorFactory::EnsureCubicChunkSize(PGeneratorConfig& config, const PEID size) const {
    if (config.k == 0) {
        config.k = FindCubeMultipleOf(size);
    } else if (!IsCubic(config.k)) {
        throw ConfigurationError("number of chunks must be cubic");
    }
}

void GeneratorFactory::EnsureOneChunkPerPE(PGeneratorConfig& config, const PEID size) const {
    if (config.k != static_cast<SInt>(size)) {
        throw ConfigurationError("number of chunks must match the number of PEs");
    }
}
} // namespace kagen
