/* Copyright (C) 2009-2010 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#pragma once

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "hash.hpp"
#include "kagen/definitions.h"
#include "kagen/generator_config.h"
#include "kagen/generators/kronecker/splittable_mrg.h"
#include "kagen/generators/kronecker/user_settings.h"
#include "kagen/generators/kronecker/utils.h"
#include "kagen/io/generator_io.h"

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

#ifdef GENERATOR_USE_PACKED_EDGE_TYPE

typedef struct packed_edge {
    uint32_t v0_low;
    uint32_t v1_low;
    uint32_t high; /* v1 in high half, v0 in low half */
} packed_edge;

static inline int64_t get_v0_from_edge(const packed_edge* p) {
    return (p->v0_low | ((int64_t)((int16_t)(p->high & 0xFFFF)) << 32));
}

static inline int64_t get_v1_from_edge(const packed_edge* p) {
    return (p->v1_low | ((int64_t)((int16_t)(p->high >> 16)) << 32));
}

static inline void write_edge(packed_edge* p, int64_t v0, int64_t v1) {
    p->v0_low = (uint32_t)v0;
    p->v1_low = (uint32_t)v1;
    p->high   = (uint32_t)(((v0 >> 32) & 0xFFFF) | (((v1 >> 32) & 0xFFFF) << 16));
}

#else

typedef struct packed_edge {
    int64_t v0;
    int64_t v1;
} packed_edge;

static inline int64_t get_v0_from_edge(const packed_edge* p) {
    return p->v0;
}

static inline int64_t get_v1_from_edge(const packed_edge* p) {
    return p->v1;
}

static inline void write_edge(packed_edge* p, int64_t v0, int64_t v1) {
    p->v0 = v0;
    p->v1 = v1;
}

#endif

class Kronecker {
public:
    Kronecker(PGeneratorConfig& config, const PEID rank, const PEID /* size */)
        : config_(config),
          rank_(rank),
          io_(config) {
        MPI_Comm_size(MPI_COMM_WORLD, &size_);

        // Init variables
        from_                = 0;
        to_                  = config_.n - 1;
        log_n_               = log2(config_.n);
        edges_per_pe_        = config_.m / config_.k;
        SInt remaining_edges = config_.m % config_.k;
        num_edges_           = edges_per_pe_ + ((SInt)rank < remaining_edges);
    }

    void Generate() {
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

    GeneratorIO& IO() {
        return io_;
    }

    std::pair<SInt, SInt> GetVertexRange() {
        return std::make_pair(from_, to_ + 1);
    }

private:
    // Config
    PGeneratorConfig& config_;
    PEID              size_, rank_;

    // I/O
    GeneratorIO io_;

    // Constants and variables
    int     log_n_;
    SInt    from_, to_;
    SInt    num_edges_;
    int64_t scramble1_, scramble2_;
    SInt    edges_per_pe_;

    int Bernoulli(mrg_state* st, int level, int nlevels) {
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
    static inline uint64_t bitreverse(uint64_t x) {
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
    inline int64_t Scramble(int64_t v0) {
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
    void GenerateEdge(int64_t n, int level, mrg_state* st) {
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
        io_.PushEdge(Scramble(base_src), Scramble(base_tgt));
    }
};
} // namespace kagen
