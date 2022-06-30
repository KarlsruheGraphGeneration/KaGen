#include "kagen/generators/rmat/rmat.h"
#include "kagen/context.h"

#include "kagen/generators/rmat/rmat_impl.hpp"

namespace kagen {
std::unique_ptr<Generator> RMATFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<RMAT>(config, rank, size);
}

RMAT::RMAT(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size) {}

void RMAT::GenerateImpl() {}
} // namespace kagen
