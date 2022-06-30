#include "kagen/generators/rmat/rmat.h"
#include "kagen/context.h"

#include "kagen/generators/rmat/generators/select.hpp"
#include "kagen/generators/rmat/rmat_impl.hpp"

namespace kagen {
std::unique_ptr<Generator> RMATFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<RMAT>(config, rank, size);
}

RMAT::RMAT(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size) {}

void RMAT::GenerateImpl() {
    using RNG  = rmat::generators::select_t;
    using RMAT = rmat::rmat<false>;

    const SInt seed  = rank_ + 1;
    const SInt log_n = std::log2(config_.n);
    const SInt depth = std::min<SInt>(9, log_n);

    RNG  gen(seed);
    RNG  gen_scramble(seed + 1000);
    RMAT r(gen_scramble, log_n, config_.rmat_a, config_.rmat_b, config_.rmat_c);
    r.init(depth);
    r.get_edges([&](const auto u, const auto v) { PushEdge(u, v); }, config_.m, gen);
}
} // namespace kagen
