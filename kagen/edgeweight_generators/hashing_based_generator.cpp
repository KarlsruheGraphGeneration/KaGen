#include "kagen/edgeweight_generators/hashing_based_generator.h"

#include "kagen/context.h"
#include "kagen/kagen.h"

#include <array>

#ifdef KAGEN_XXHASH_FOUND
    #include "xxhash.h"
#else // KAGEN_XXHASH_FOUND
    #include <stdexcept>
#endif // KAGEN_XXHASH_FOUND

namespace kagen {
HashingBasedEdgeWeightGenerator::HashingBasedEdgeWeightGenerator(EdgeWeightConfig config, VertexRange vertex_range)
    : PerEdgeWeightGenerator(vertex_range),
      config_(config) {}

SSInt HashingBasedEdgeWeightGenerator::GenerateEdgeWeight(SInt u, SInt v) {
#ifdef KAGEN_XXHASH_FOUND
    if (u > v) {
        std::swap(u, v);
    }
    std::array<SInt, 2> buf{u, v};
    const SInt          range = config_.weight_range_end - config_.weight_range_begin;
    const SInt          hash  = XXH64(buf.data(), sizeof(buf), 0);
    #if defined(__SIZEOF_INT128__)
    const SInt hash_in_range = (static_cast<__uint128_t>(hash) * range >> 64);
    #else
    const SInt hash_in_range = hash % range;
    #endif
    return hash_in_range + config_.weight_range_begin;
#else  // KAGEN_XXHASH_FOUND
    ((void)u);
    ((void)v);
    throw std::runtime_error("xxHash is required for hashing based edge weights");
#endif // KAGEN_XXHASH_FOUND
    return -1;
}
} // namespace kagen
