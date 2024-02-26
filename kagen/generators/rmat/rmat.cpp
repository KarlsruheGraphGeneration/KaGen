#include "kagen/generators/rmat/rmat.h"

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/generators/rmat/generators/select.hpp"
#include "kagen/generators/rmat/rmat_impl.hpp"

#include <mpi.h>

namespace kagen {
PGeneratorConfig
RMATFactory::NormalizeParameters(PGeneratorConfig config, PEID, const PEID size, const bool output) const {
    if (config.rmat_a < 0 || config.rmat_b < 0 || config.rmat_b < 0) {
        throw ConfigurationError("probabilities may not be negative");
    }
    if (config.rmat_a + config.rmat_b + config.rmat_c > 1) {
        throw ConfigurationError("sum of probabilities may not be larger than 1");
    }

    const SInt log_n = std::log2(config.n);
    if (log_n > 31) {
        throw ConfigurationError("number of vertices is too large (cannot be larger than 31 bits)");
    }

    if (output && config.n != 1ull << log_n) {
        std::cout << "Warning: generator requires the number of vertices to be a power of two" << std::endl;
        std::cout << "  Changing the number of vertices to " << (1ull << log_n) << std::endl;
    }

    if (config.k == 0) {
        config.k = static_cast<SInt>(size);
    }

    return config;
}

std::unique_ptr<Generator> RMATFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<RMAT>(config, rank, size);
}

RMAT::RMAT(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : Graph500Generator(config),
      config_(config),
      rank_(rank) {
    const SInt edges_per_pe    = config_.m / size;
    const SInt remaining_edges = config_.m % size;
    num_edges_                 = edges_per_pe + ((SInt)rank < remaining_edges);
}

void RMAT::GenerateEdgeList() {
    using RNG  = rmat::generators::select_t;
    using RMAT = rmat::rmat<false>;

    const SInt seed  = rank_ + config_.seed;
    const SInt log_n = std::log2(config_.n);
    const SInt depth = std::min<SInt>(9, log_n);

    RNG  gen(seed);
    RNG  gen_scramble(seed + 1000);
    RMAT r(gen_scramble, log_n, config_.rmat_a, config_.rmat_b, config_.rmat_c);
    r.init(depth);

    // Generate local edges
    r.get_edges([&](const auto u, const auto v) { PushLocalEdge(u, v); }, num_edges_, gen);
}
} // namespace kagen
