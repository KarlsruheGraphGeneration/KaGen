#include "kagen/generators/kronecker/kronecker.h"

#include "kagen/generators/generator.h"
#include "kagen/generators/kronecker/splittable_mrg.h"
#include "kagen/generators/kronecker/user_settings.h"
#include "kagen/generators/kronecker/utils.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <cmath>

#ifndef __STDC_FORMAT_MACROS
    #define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

/* Initiator settings: for faster random number generation, the initiator
 * probabilities are defined as fractions (a = INITIATOR_A_NUMERATOR /
 * INITIATOR_DENOMINATOR, b = c = INITIATOR_BC_NUMERATOR /
 * INITIATOR_DENOMINATOR, d = 1 - a - b - c. */
#define INITIATOR_A_NUMERATOR  2500
#define INITIATOR_BC_NUMERATOR 2500
#define INITIATOR_DENOMINATOR  10000

/* If this macro is defined to a non-zero value, use SPK_NOISE_LEVEL /
 * INITIATOR_DENOMINATOR as the noise parameter to use in introducing noise
 * into the graph parameters.  The approach used is from "A Hitchhiker's Guide
 * to Choosing Parameters of Stochastic Kronecker Graphs" by C. Seshadhri, Ali
 * Pinar, and Tamara G. Kolda (http://arxiv.org/abs/1102.5046v1), except that
 * the adjustment here is chosen based on the current level being processed
 * rather than being chosen randomly. */
#define SPK_NOISE_LEVEL 0
/* #define SPK_NOISE_LEVEL 1000 -- in INITIATOR_DENOMINATOR units */

namespace kagen {
std::unique_ptr<Generator>
KroneckerFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<Kronecker>(config, rank, size);
}

Kronecker::Kronecker(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : Graph500Generator(config),
      config_(config),
      size_(size),
      rank_(rank) {
    log_n_                     = std::log2(config_.n);
    const SInt edges_per_pe    = config_.m / size;
    const SInt remaining_edges = config_.m % size;
    num_edges_                 = edges_per_pe + ((SInt)rank < remaining_edges);
}

void Kronecker::GenerateEdgeList() {
    uint_fast32_t seed[5];
    make_mrg_seed(sampling::Spooky::hash((config_.seed + 1) * size_), sampling::Spooky::hash(rank_), seed);

    mrg_state state;

    mrg_seed(&state, seed);

    [[maybe_unused]] uint64_t scramble1_, scramble2_; /* Values for scrambling */
    {
        mrg_state new_state = state;
        mrg_skip(&new_state, 50, 7, 0);
        scramble1_ = mrg_get_uint_orig(&new_state);
        scramble1_ *= UINT64_C(0xFFFFFFFF);
        scramble1_ += mrg_get_uint_orig(&new_state);
        scramble2_ = mrg_get_uint_orig(&new_state);
        scramble2_ *= UINT64_C(0xFFFFFFFF);
        scramble2_ += mrg_get_uint_orig(&new_state);
    }

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (SInt i = 0; i < num_edges_; ++i) {
        mrg_state new_state = state;
        mrg_skip(&new_state, 0, (uint64_t)i, 0);
        GenerateEdge(config_.n, 0, &new_state);
    }
}

int Kronecker::Bernoulli(mrg_state* st, int level, int nlevels) {
#if SPK_NOISE_LEVEL == 0
    /* Avoid warnings */
    (void)level;
    (void)nlevels;
#endif
    /* Generate a pseudorandom number in the range [0, INITIATOR_DENOMINATOR)
     * without modulo bias. */
    static const uint32_t limit = (UINT32_C(0x7FFFFFFF) % INITIATOR_DENOMINATOR);
    uint32_t              val   = mrg_get_uint_orig(st);
    if (/* Unlikely */ val < limit) {
        do {
            val = mrg_get_uint_orig(st);
        } while (val < limit);
    }
#if SPK_NOISE_LEVEL == 0
    int spk_noise_factor = 0;
#else
    int spk_noise_factor = 2 * SPK_NOISE_LEVEL * level / nlevels - SPK_NOISE_LEVEL;
#endif
    unsigned int adjusted_bc_numerator = (unsigned int)(INITIATOR_BC_NUMERATOR + spk_noise_factor);
    val %= INITIATOR_DENOMINATOR;
    if (val < adjusted_bc_numerator)
        return 1;
    val = (uint32_t)(val - adjusted_bc_numerator);
    if (val < adjusted_bc_numerator)
        return 2;
    val = (uint32_t)(val - adjusted_bc_numerator);
#if SPK_NOISE_LEVEL == 0
    if (val < INITIATOR_A_NUMERATOR)
        return 0;
#else
    if (val < INITIATOR_A_NUMERATOR * (INITIATOR_DENOMINATOR - 2 * INITIATOR_BC_NUMERATOR)
                  / (INITIATOR_DENOMINATOR - 2 * adjusted_bc_numerator))
        return 0;
#endif
#if SPK_NOISE_LEVEL == 0
    /* Avoid warnings */
    (void)level;
    (void)nlevels;
#endif
    return 3;
}

/* Reverse bits in a number; this should be optimized for performance
 * (including using bit- or byte-reverse intrinsics if your platform has them).
 * */
inline uint64_t Kronecker::bitreverse(uint64_t x) {
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 3)
    #define USE_GCC_BYTESWAP /* __builtin_bswap* are in 4.3 but not 4.2 */
#endif

#ifdef FAST_64BIT_ARITHMETIC

    /* 64-bit code */
    #ifdef USE_GCC_BYTESWAP
    x = __builtin_bswap64(x);
    #else
    x = (x >> 32) | (x << 32);
    x = ((x >> 16) & UINT64_C(0x0000FFFF0000FFFF)) | ((x & UINT64_C(0x0000FFFF0000FFFF)) << 16);
    x = ((x >> 8) & UINT64_C(0x00FF00FF00FF00FF)) | ((x & UINT64_C(0x00FF00FF00FF00FF)) << 8);
    #endif
    x = ((x >> 4) & UINT64_C(0x0F0F0F0F0F0F0F0F)) | ((x & UINT64_C(0x0F0F0F0F0F0F0F0F)) << 4);
    x = ((x >> 2) & UINT64_C(0x3333333333333333)) | ((x & UINT64_C(0x3333333333333333)) << 2);
    x = ((x >> 1) & UINT64_C(0x5555555555555555)) | ((x & UINT64_C(0x5555555555555555)) << 1);
    return x;

#else

    /* 32-bit code */
    uint32_t h = (uint32_t)(x >> 32);
    uint32_t l = (uint32_t)(x & UINT32_MAX);
    #ifdef USE_GCC_BYTESWAP
    h          = __builtin_bswap32(h);
    l          = __builtin_bswap32(l);
    #else
    h = (h >> 16) | (h << 16);
    l = (l >> 16) | (l << 16);
    h = ((h >> 8) & UINT32_C(0x00FF00FF)) | ((h & UINT32_C(0x00FF00FF)) << 8);
    l = ((l >> 8) & UINT32_C(0x00FF00FF)) | ((l & UINT32_C(0x00FF00FF)) << 8);
    #endif
    h          = ((h >> 4) & UINT32_C(0x0F0F0F0F)) | ((h & UINT32_C(0x0F0F0F0F)) << 4);
    l          = ((l >> 4) & UINT32_C(0x0F0F0F0F)) | ((l & UINT32_C(0x0F0F0F0F)) << 4);
    h          = ((h >> 2) & UINT32_C(0x33333333)) | ((h & UINT32_C(0x33333333)) << 2);
    l          = ((l >> 2) & UINT32_C(0x33333333)) | ((l & UINT32_C(0x33333333)) << 2);
    h          = ((h >> 1) & UINT32_C(0x55555555)) | ((h & UINT32_C(0x55555555)) << 1);
    l          = ((l >> 1) & UINT32_C(0x55555555)) | ((l & UINT32_C(0x55555555)) << 1);
    return ((uint64_t)l << 32) | h; /* Swap halves */

#endif
}

/* Apply a permutation to scramble vertex numbers; a randomly generated
 * permutation is not used because applying it at scale is too expensive. */
inline int64_t Kronecker::Scramble(int64_t v0) {
    uint64_t v = (uint64_t)v0;
    v += scramble1_ + scramble2_;
    v *= (scramble1_ | UINT64_C(0x4519840211493211));
    v = (bitreverse(v) >> (64 - log_n_));
    assert((v >> log_n_) == 0);
    v *= (scramble2_ | UINT64_C(0x3050852102C843A5));
    v = (bitreverse(v) >> (64 - log_n_));
    assert((v >> log_n_) == 0);
    return (int64_t)v;
}

/* Make a single graph edge using a pre-set MRG state. */
void Kronecker::GenerateEdge(int64_t n, int level, mrg_state* st) {
    int64_t base_src = 0, base_tgt = 0;
    while (n > 1) {
        int square     = Bernoulli(st, level, log_n_);
        int src_offset = square / 2;
        int tgt_offset = square % 2;
        assert(base_src <= base_tgt);
        if (base_src == base_tgt) {
            /* Clip-and-flip for undirected graph */
            if (src_offset > tgt_offset) {
                int temp   = src_offset;
                src_offset = tgt_offset;
                tgt_offset = temp;
            }
        }
        n /= 2;
        ++level;
        base_src += n * src_offset;
        base_tgt += n * tgt_offset;
    }

    const auto u = Scramble(base_src);
    const auto v = Scramble(base_tgt);
    PushLocalEdge(u, v);
}
} // namespace kagen
