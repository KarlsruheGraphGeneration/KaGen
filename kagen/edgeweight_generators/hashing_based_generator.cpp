#include "kagen/edgeweight_generators/hashing_based_generator.h"

#include "kagen/context.h"
#include "kagen/kagen.h"

#ifdef KAGEN_XXHASH_FOUND
    #include "xxhash.h"
#else // KAGEN_XXHASH_FOUND
    #include <stdexcept>
#endif // KAGEN_XXHASH_FOUND

namespace kagen {
HashingBasedEdgeWeightGenerator::HashingBasedEdgeWeightGenerator(EdgeWeightConfig config) : config_(config) {}

SSInt HashingBasedEdgeWeightGenerator::GenerateEdgeWeight(SInt u, SInt v) {
#ifdef KAGEN_XXHASH_FOUND
    const SInt  hash1 = XXH64(&u, 1, 0);
    const SInt  hash2 = XXH64(&v, 1, 0);
    const SSInt combined =
        (hash1 ^ hash2) % (config_.weight_range_end - config_.weight_range_begin) + config_.weight_range_begin;
    return combined;
#else  // KAGEN_XXHASH_FOUND
    ((void)u);
    ((void)v);
    throw std::runtime_error("xxHash is required for hashing based edge weights");
#endif // KAGEN_XXHASH_FOUND
}
} // namespace kagen
