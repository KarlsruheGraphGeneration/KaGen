/*******************************************************************************
 * rmat/graph500.hpp
 *
 * Graph500 bits used for R-MAT generation
 *
 * Copyright (C) 2018-2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef RMAT_GRAPH500_HEADER
    #define RMAT_GRAPH500_HEADER

    #include "kagen/tlx/attribute_always_inline.hpp"

    // generate unweighted graphs
    #ifdef SSSP
        #undef SSSP
    #endif // SSSP

    #define restrict __restrict__

// graph500 reference code from https://github.com/graph500/graph500
namespace {
    #include "kagen/generators/rmat/graph500/graph_generator.c"
    #include "kagen/generators/rmat/graph500/splittable_mrg.c"
} // namespace

namespace rmat {
namespace graph500 {
struct scramble_state {
    uint64_t val0;
    uint64_t val1;
    int      lgN;
    scramble_state() = delete;
    scramble_state(uint64_t v0, uint64_t v1, int lgN_) : val0(v0), val1(v1), lgN(lgN_) {}
};

/* adapted from generator/graph_generator.c generate_kronecker_range */
inline scramble_state init_scramble_helper(const uint_fast32_t seed[5], int lgN) {
    mrg_state state;
    mrg_seed(&state, seed);

    uint64_t val0, val1; /* Values for scrambling */
    {
        mrg_state new_state = state;
        mrg_skip(&new_state, 50, 7, 0);

        val0 = mrg_get_uint_orig(&new_state);
        val0 *= UINT64_C(0xFFFFFFFF);
        val0 += mrg_get_uint_orig(&new_state);
        val1 = mrg_get_uint_orig(&new_state);
        val1 *= UINT64_C(0xFFFFFFFF);
        val1 += mrg_get_uint_orig(&new_state);
    }

    return scramble_state(val0, val1, lgN);
}

// init scramble state, generating seed from RNG
template <typename RNG>
scramble_state init_scramble_state(RNG& rng, int lgN) {
    uint_fast32_t seed[5]; /* All values in [0, 2^31 - 1), not all zero */
    for (int i = 0; i < 5; i++) {
        seed[i] = static_cast<uint_fast32_t>(rng.next() * ((static_cast<uint_fast32_t>(1) << 31) - 1));
    }

    return init_scramble_helper(seed, lgN);
}

/* from graph500 generator/graph_generator.c */
/* Reverse bits in a number; this should be optimized for performance
 * (including using bit- or byte-reverse intrinsics if your platform has them).
 * */
constexpr uint64_t bitreverse(uint64_t x) {
    /* 64-bit code */
    x = __builtin_bswap64(x);
    x = ((x >> 4) & UINT64_C(0x0F0F0F0F0F0F0F0F)) | ((x & UINT64_C(0x0F0F0F0F0F0F0F0F)) << 4);
    x = ((x >> 2) & UINT64_C(0x3333333333333333)) | ((x & UINT64_C(0x3333333333333333)) << 2);
    x = ((x >> 1) & UINT64_C(0x5555555555555555)) | ((x & UINT64_C(0x5555555555555555)) << 1);
    return x;
}

constexpr void bitreverse_two(uint64_t& x, uint64_t& y) {
    /* 64-bit code */
    x = __builtin_bswap64(x);
    y = __builtin_bswap64(y);
    x = ((x >> 4) & UINT64_C(0x0F0F0F0F0F0F0F0F)) | ((x & UINT64_C(0x0F0F0F0F0F0F0F0F)) << 4);
    y = ((y >> 4) & UINT64_C(0x0F0F0F0F0F0F0F0F)) | ((y & UINT64_C(0x0F0F0F0F0F0F0F0F)) << 4);
    x = ((x >> 2) & UINT64_C(0x3333333333333333)) | ((x & UINT64_C(0x3333333333333333)) << 2);
    y = ((y >> 2) & UINT64_C(0x3333333333333333)) | ((y & UINT64_C(0x3333333333333333)) << 2);
    x = ((x >> 1) & UINT64_C(0x5555555555555555)) | ((x & UINT64_C(0x5555555555555555)) << 1);
    y = ((y >> 1) & UINT64_C(0x5555555555555555)) | ((y & UINT64_C(0x5555555555555555)) << 1);
}

TLX_ATTRIBUTE_ALWAYS_INLINE
constexpr void scramble_two(int64_t& v0, int64_t& w0, int lgN, uint64_t val0, uint64_t val1) {
    uint64_t v = (uint64_t)v0;
    uint64_t w = (uint64_t)w0;
    v += val0 + val1;
    w += val0 + val1;
    v *= (val0 | UINT64_C(0x4519840211493211));
    w *= (val0 | UINT64_C(0x4519840211493211));
    /*
    v = (bitreverse(v) >> (64 - lgN));
    w = (bitreverse(w) >> (64 - lgN));
    */
    bitreverse_two(v, w);
    v = (v >> (64 - lgN));
    w = (w >> (64 - lgN));

    assert((v >> lgN) == 0);
    assert((w >> lgN) == 0);
    v *= (val1 | UINT64_C(0x3050852102C843A5));
    w *= (val1 | UINT64_C(0x3050852102C843A5));
    v = (bitreverse(v) >> (64 - lgN));
    w = (bitreverse(w) >> (64 - lgN));
    assert((v >> lgN) == 0);
    assert((w >> lgN) == 0);
    v0 = (int64_t)v;
    w0 = (int64_t)w;
}

TLX_ATTRIBUTE_ALWAYS_INLINE
constexpr void scramble_two(int64_t& v0, int64_t& w0, const scramble_state& state) {
    scramble_two(v0, w0, state.lgN, state.val0, state.val1);
}
} // namespace graph500
} // namespace rmat
#endif // RMAT_GRAPH500_HEADER
