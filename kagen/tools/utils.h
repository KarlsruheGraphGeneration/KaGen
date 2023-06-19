#pragma once

#include "kagen/definitions.h"

namespace kagen {
inline std::pair<SInt, SInt> ComputeRange(const SInt n, const PEID size, const PEID rank) {
    const SInt chunk = n / size;
    const SInt rem   = n % size;
    const SInt from  = rank * chunk + std::min<SInt>(rank, rem);
    const SInt to    = std::min<SInt>(from + ((static_cast<SInt>(rank) < rem) ? chunk + 1 : chunk), n);
    return {from, to};
}
} // namespace kagen
