#include "kagen/edgeweight_generators/hashing_based_generator.h"

#include "kagen/context.h"
#include "kagen/kagen.h"

#include "xxhash.h"

namespace kagen {
HashingBasedEdgeWeightGenerator::HashingBasedEdgeWeightGenerator(EdgeWeightConfig config) : config_(config) {}

SSInt HashingBasedEdgeWeightGenerator::GenerateEdgeWeight(SInt u, SInt v) {
    SSInt hash1, hash2;

    hash1 = XXH64(&u, 1, 0);
    hash2 = XXH64(&v, 1, 0);
    SSInt combined_value =
        (hash1 ^ hash2) % (config_.weight_range_end - config_.weight_range_begin) + config_.weight_range_begin;
    return combined_value;
}
} // namespace kagen
