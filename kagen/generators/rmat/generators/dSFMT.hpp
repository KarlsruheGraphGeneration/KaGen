/*******************************************************************************
 * rmat/generators/dSFMT.hpp
 *
 * A double-precision SIMD-oriented Fast Mersenne Twister
 *
 * Generates double precision floating point pseudorandom numbers which
 * distribute in the range of [0, 1) with period 19937
 *
 * Adapted from Mutsuo Saito and Makoto Matsumoto's dSFMT 2.2.3, available at
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/ under the following
 * license:
 *
 * Copyright (c) 2007, 2008, 2009 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University.
 * Copyright (c) 2011, 2002 Mutsuo Saito, Makoto Matsumoto, Hiroshima University
 * and The University of Tokyo.  All rights reserved.
 * Copyright (c) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above
 *    copyright notice, this list of conditions and the following
 *    disclaimer in the documentation and/or other materials provided
 *    with the distribution.
 *  * Neither the name of the Hiroshima University nor the names of
 *    its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written
 *    permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ******************************************************************************/


/*
 * inlined: dSFMT.h
 */

#pragma once
#ifndef RMAT_GENERATORS_DSFMT_HEADER
#define RMAT_GENERATORS_DSFMT_HEADER

#include <tlx/logger.hpp>
#include <tlx/define.hpp>

#include <cassert>
#include <cmath>
#include <limits>
#include <vector>


#include <inttypes.h>
#if defined(HAVE_ALTIVEC) && !defined(__APPLE__)
#  include <altivec.h>
#elif defined(SAMPLING_HAVE_SSE2)
#  include <emmintrin.h>
#endif

namespace rmat {
namespace generators {
namespace _dSFMT {

extern "C" {

// this adaptation only supports MT19937
#define DSFMT_MEXP 19937

/*-----------------
  BASIC DEFINITIONS
  -----------------*/
/* Mersenne Exponent. The period of the sequence
 *  is a multiple of 2^DSFMT_MEXP-1.
 * #define DSFMT_MEXP 19937 */
/** DSFMT generator has an internal state array of 128-bit integers,
 * and N is its size. */
#define DSFMT_N ((DSFMT_MEXP - 128) / 104 + 1)
/** N64 is the size of internal state array when regarded as an array
 * of 64-bit integers.*/
#define DSFMT_N64 (DSFMT_N * 2)

#if !defined(DSFMT_BIG_ENDIAN)
#  if defined(__BYTE_ORDER) && defined(__BIG_ENDIAN)
#    if __BYTE_ORDER == __BIG_ENDIAN
#      define DSFMT_BIG_ENDIAN 1
#    endif
#  elif defined(_BYTE_ORDER) && defined(_BIG_ENDIAN)
#    if _BYTE_ORDER == _BIG_ENDIAN
#      define DSFMT_BIG_ENDIAN 1
#    endif
#  elif defined(__BYTE_ORDER__) && defined(__BIG_ENDIAN__)
#    if __BYTE_ORDER__ == __BIG_ENDIAN__
#      define DSFMT_BIG_ENDIAN 1
#    endif
#  elif defined(BYTE_ORDER) && defined(BIG_ENDIAN)
#    if BYTE_ORDER == BIG_ENDIAN
#      define DSFMT_BIG_ENDIAN 1
#    endif
#  elif defined(__BIG_ENDIAN) || defined(_BIG_ENDIAN) \
    || defined(__BIG_ENDIAN__) || defined(BIG_ENDIAN)
#      define DSFMT_BIG_ENDIAN 1
#  endif
#endif

#if defined(DSFMT_BIG_ENDIAN) && defined(__amd64)
#  undef DSFMT_BIG_ENDIAN
#endif

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
//#  include <inttypes.h>
#elif defined(_MSC_VER) || defined(__BORLANDC__)
#  if !defined(DSFMT_UINT32_DEFINED) && !defined(SFMT_UINT32_DEFINED)
typedef unsigned int uint32_t;
typedef unsigned __int64 uint64_t;
#    ifndef UINT64_C
#      define UINT64_C(v) (v ## ui64)
#    endif
#    define DSFMT_UINT32_DEFINED
#    if !defined(inline) && !defined(__cplusplus)
#      define inline __inline
#    endif
#  endif
#else
//#  include <inttypes.h>
#  if !defined(inline) && !defined(__cplusplus)
#    if defined(__GNUC__)
#      define inline __inline__
#    else
#      define inline
#    endif
#  endif
#endif

#ifndef UINT64_C
#  define UINT64_C(v) (v ## ULL)
#endif

/*------------------------------------------
  128-bit SIMD like data type for standard C
  ------------------------------------------*/
#if defined(HAVE_ALTIVEC)
#  if !defined(__APPLE__)
//#    include <altivec.h>
#  endif
/** 128-bit data structure */
union W128_T {
    vector unsigned int s;
    uint64_t u[2];
    uint32_t u32[4];
    double d[2];
};

#elif defined(SAMPLING_HAVE_SSE2)
//#  include <emmintrin.h>

/** 128-bit data structure */
union W128_T {
    __m128i si;
    __m128d sd;
    uint64_t u[2];
    uint32_t u32[4];
    double d[2];
};
#else  /* standard C */
/** 128-bit data structure */
union W128_T {
    uint64_t u[2];
    uint32_t u32[4];
    double d[2];
};
#endif

/** 128-bit data type */
typedef union W128_T w128_t;

/** the 128-bit internal state array */
struct DSFMT_T {
    w128_t status[DSFMT_N + 1];
    int idx;
};
typedef struct DSFMT_T dsfmt_t;

void dsfmt_fill_array_close_open(dsfmt_t *dsfmt, double array[], int size);
void dsfmt_chk_init_gen_rand(dsfmt_t *dsfmt, uint32_t seed, int mexp);
int dsfmt_get_min_array_size(void);

/**
 * This function initializes the internal state array with a 32-bit
 * integer seed.
 * @param dsfmt dsfmt state vector.
 * @param seed a 32-bit integer used as the seed.
 */
inline static void dsfmt_init_gen_rand(dsfmt_t *dsfmt, uint32_t seed) {
    dsfmt_chk_init_gen_rand(dsfmt, seed, DSFMT_MEXP);
}

/*
 * inlined: dSFMT-params.h
 */
#define DSFMT_LOW_MASK  UINT64_C(0x000FFFFFFFFFFFFF)
#define DSFMT_HIGH_CONST UINT64_C(0x3FF0000000000000)
#define DSFMT_SR	12

/* for sse2 */
#if defined(SAMPLING_HAVE_SSE2)
  #define SSE2_SHUFF 0x1b
#elif defined(HAVE_ALTIVEC)
  #if defined(__APPLE__)  /* For OSX */
    #define ALTI_SR (vector unsigned char)(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)
    #define ALTI_SR_PERM \
        (vector unsigned char)(15,0,1,2,3,4,5,6,15,8,9,10,11,12,13,14)
    #define ALTI_SR_MSK \
        (vector unsigned int)(0x000fffffU,0xffffffffU,0x000fffffU,0xffffffffU)
    #define ALTI_PERM \
        (vector unsigned char)(12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3)
  #else
    #define ALTI_SR      {4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4}
    #define ALTI_SR_PERM {15,0,1,2,3,4,5,6,15,8,9,10,11,12,13,14}
    #define ALTI_SR_MSK  {0x000fffffU,0xffffffffU,0x000fffffU,0xffffffffU}
    #define ALTI_PERM    {12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3}
  #endif
#endif

/*
 * inlined: dSFMT-params19937.h
 */

/* #define DSFMT_N	191 */
/* #define DSFMT_MAXDEGREE	19992 */
#define DSFMT_POS1	117
#define DSFMT_SL1	19
#define DSFMT_MSK1	UINT64_C(0x000ffafffffffb3f)
#define DSFMT_MSK2	UINT64_C(0x000ffdfffc90fffd)
#define DSFMT_MSK32_1	0x000ffaffU
#define DSFMT_MSK32_2	0xfffffb3fU
#define DSFMT_MSK32_3	0x000ffdffU
#define DSFMT_MSK32_4	0xfc90fffdU
#define DSFMT_FIX1	UINT64_C(0x90014964b32f4329)
#define DSFMT_FIX2	UINT64_C(0x3b8d12ac548a7c7a)
#define DSFMT_PCV1	UINT64_C(0x3d84e1ac0dc82880)
#define DSFMT_PCV2	UINT64_C(0x0000000000000001)


/* PARAMETERS FOR ALTIVEC */
#if defined(__APPLE__)	/* For OSX */
    #define ALTI_SL1 	(vector unsigned char)(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)
    #define ALTI_SL1_PERM \
	(vector unsigned char)(2,3,4,5,6,7,30,30,10,11,12,13,14,15,0,1)
    #define ALTI_SL1_MSK \
	(vector unsigned int)(0xffffffffU,0xfff80000U,0xffffffffU,0xfff80000U)
    #define ALTI_MSK	(vector unsigned int)(DSFMT_MSK32_1, \
			DSFMT_MSK32_2, DSFMT_MSK32_3, DSFMT_MSK32_4)
#else	/* For OTHER OSs(Linux?) */
    #define ALTI_SL1 	{3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3}
    #define ALTI_SL1_PERM \
	{2,3,4,5,6,7,30,30,10,11,12,13,14,15,0,1}
    #define ALTI_SL1_MSK \
	{0xffffffffU,0xfff80000U,0xffffffffU,0xfff80000U}
    #define ALTI_MSK \
	{DSFMT_MSK32_1, DSFMT_MSK32_2, DSFMT_MSK32_3, DSFMT_MSK32_4}
#endif

/*
 * inlined: dSFMT-common.h
 */

#if defined(SAMPLING_HAVE_SSE2)
//#  include <emmintrin.h>
union X128I_T {
    uint64_t u[2];
    __m128i  i128;
};
union X128D_T {
    double d[2];
    __m128d d128;
};
/** mask data for sse2 */
static const union X128I_T sse2_param_mask = {{DSFMT_MSK1, DSFMT_MSK2}};
#endif

#if defined(HAVE_ALTIVEC)
inline static void do_recursion(w128_t *r, w128_t *a, w128_t * b,
				w128_t *lung) {
    const vector unsigned char sl1 = ALTI_SL1;
    const vector unsigned char sl1_perm = ALTI_SL1_PERM;
    const vector unsigned int sl1_msk = ALTI_SL1_MSK;
    const vector unsigned char sr1 = ALTI_SR;
    const vector unsigned char sr1_perm = ALTI_SR_PERM;
    const vector unsigned int sr1_msk = ALTI_SR_MSK;
    const vector unsigned char perm = ALTI_PERM;
    const vector unsigned int msk1 = ALTI_MSK;
    vector unsigned int w, x, y, z;

    z = a->s;
    w = lung->s;
    x = vec_perm(w, (vector unsigned int)perm, perm);
    y = vec_perm(z, (vector unsigned int)sl1_perm, sl1_perm);
    y = vec_sll(y, sl1);
    y = vec_and(y, sl1_msk);
    w = vec_xor(x, b->s);
    w = vec_xor(w, y);
    x = vec_perm(w, (vector unsigned int)sr1_perm, sr1_perm);
    x = vec_srl(x, sr1);
    x = vec_and(x, sr1_msk);
    y = vec_and(w, msk1);
    z = vec_xor(z, y);
    r->s = vec_xor(z, x);
    lung->s = w;
}
#elif defined(SAMPLING_HAVE_SSE2)
/**
 * This function represents the recursion formula.
 * @param r output 128-bit
 * @param a a 128-bit part of the internal state array
 * @param b a 128-bit part of the internal state array
 * @param d a 128-bit part of the internal state array (I/O)
 */
inline static void do_recursion(w128_t *r, w128_t *a, w128_t *b, w128_t *u) {
    __m128i v, w, x, y, z;

    x = a->si;
    z = _mm_slli_epi64(x, DSFMT_SL1);
    y = _mm_shuffle_epi32(u->si, SSE2_SHUFF);
    z = _mm_xor_si128(z, b->si);
    y = _mm_xor_si128(y, z);

    v = _mm_srli_epi64(y, DSFMT_SR);
    w = _mm_and_si128(y, sse2_param_mask.i128);
    v = _mm_xor_si128(v, x);
    v = _mm_xor_si128(v, w);
    r->si = v;
    u->si = y;
}
#else
/**
 * This function represents the recursion formula.
 * @param r output 128-bit
 * @param a a 128-bit part of the internal state array
 * @param b a 128-bit part of the internal state array
 * @param lung a 128-bit part of the internal state array (I/O)
 */
inline static void do_recursion(w128_t *r, w128_t *a, w128_t * b,
				w128_t *lung) {
    uint64_t t0, t1, L0, L1;

    t0 = a->u[0];
    t1 = a->u[1];
    L0 = lung->u[0];
    L1 = lung->u[1];
    lung->u[0] = (t0 << DSFMT_SL1) ^ (L1 >> 32) ^ (L1 << 32) ^ b->u[0];
    lung->u[1] = (t1 << DSFMT_SL1) ^ (L0 >> 32) ^ (L0 << 32) ^ b->u[1];
    r->u[0] = (lung->u[0] >> DSFMT_SR) ^ (lung->u[0] & DSFMT_MSK1) ^ t0;
    r->u[1] = (lung->u[1] >> DSFMT_SR) ^ (lung->u[1] & DSFMT_MSK2) ^ t1;
}
#endif

} // extern "C"

} // namespace _dSFMT


/*!
 * A wrapper around dSFMT
 */
class dSFMT {
public:
    static constexpr bool debug = true;
    static const char* name;

    dSFMT(size_t seed) : index_(0), block_size_(0), block_id_(0) {
        _dSFMT::dsfmt_init_gen_rand(&dsfmt_, seed);
    }

    //! non-copyable: delete copy-constructor
    dSFMT(const dSFMT &) = delete;
    //! non-copyable: delete assignment operator
    dSFMT & operator = (const dSFMT &) = delete;
    //! move-constructor: default
    dSFMT(dSFMT &&) = default;
    //! move-assignment operator: default
    dSFMT & operator = (dSFMT &&) = default;

    //! Re-seed the dsfmt
    void seed(size_t seed) {
        _dSFMT::dsfmt_init_gen_rand(&dsfmt_, seed);
        // Reset all counters, too
        block_id_ = 0;
        block_size_ = 0;
        index_ = 0;
    }

    //! Minimum number of elements that needs to be generated at a time
    size_t minimum_block_size() const {
        return _dSFMT::dsfmt_get_min_array_size();
    }

    //! Minimum number of elements that needs to be generated at a time for
    //! reasonable performance
    size_t minimum_reasonable_block_size() const {
        return minimum_block_size();
    }

    //! Generate `size` [0,1) doubles in `output`
    void generate_block(std::vector<double> &output, size_t size)
    {
        // Ensure minimum block size (normally 382)
        const size_t min_size = _dSFMT::dsfmt_get_min_array_size();
        if (size < min_size) {
            sLOG << "dSFMT: requested fewer than" << min_size
                 << "deviates, namely" << size;
            size = min_size;
        }
        // resize if the output vector is too small
        if (size > output.size()) {
            output.resize(size);
        }
        _dSFMT::dsfmt_fill_array_close_open(&dsfmt_, output.data(), size);
    }

    //! Generate `size` [0,1) doubles in `output`
    void generate_block(double *arr, size_t size) {
        constexpr bool debug = false;

        // Ensure minimum block size (normally 382)
        const size_t min_size = _dSFMT::dsfmt_get_min_array_size();
        if (size < min_size) {
            sLOG << "dSFMT: requested fewer than" << min_size
                 << "deviates, namely" << size << " -- fallback to next()";
            for (size_t i = 0; i < size; i++) {
                arr[i] = next();
            }
            return;
        }

        _dSFMT::dsfmt_fill_array_close_open(&dsfmt_, arr, size);
    }

    //! Get a single [0,1) double. Computes increasingly large blocks internally
    //! so that this is fast.  May block while next block is generated.
    TLX_ATTRIBUTE_ALWAYS_INLINE
    double next() {
        if (TLX_UNLIKELY(index_ >= block_size_)) {
            if (block_id_ > 2 && ((block_id_ + 1) & block_id_) == 0) {
                // block_id_ + 1 is a power of two. We appear to need a lot of
                // random numbers, increase the blocksize to reduce RNG overhead
                block_size_ *= 2;
            }
            block_size_ = std::max(block_size_,
                                   minimum_reasonable_block_size());
            // generate_block takes care of resizing the vector for us
            generate_block(randblock_, block_size_);
            index_ = 0;
            block_id_++;
        }
        return randblock_[index_++];
    }

    TLX_ATTRIBUTE_ALWAYS_INLINE
    double next_exponential(double lambda) {
        return -std::log(1.0 - next()) / lambda;
    }

    //! Generate a uniform integer from [min, max] (i.e., both inclusive)
    template <typename int_t>
    TLX_ATTRIBUTE_ALWAYS_INLINE
    int_t next_int(int_t min, int_t max) {
        return next() * (max - min + 1) + min;
    }

    //! Bernoulli trial with success probability p
    TLX_ATTRIBUTE_ALWAYS_INLINE
    bool next_bernoulli(double p) {
        assert(0 <= p && p <= 1);
        return next() < p;
    }

    //! Bernoulli trial with success probability cutoff/max
    TLX_ATTRIBUTE_ALWAYS_INLINE
    bool next_bernoulli(double cutoff, double max) {
        assert(0 <= cutoff && cutoff <= max);
        return next() * max < cutoff;
    }

    //! Generate a normally distributed value. This needs two uniform deviates,
    //! if you need more than one, look at next_two_gaussians.
    TLX_ATTRIBUTE_ALWAYS_INLINE
    double next_gaussian(double mean, double stdev) {
        double U = next(), V = next();
        double a = stdev * std::sqrt(-2*std::log(U));
        double b = 2 * M_PI * V;

        return mean + a * std::cos(b);
    }

    //! Generate two independent normally distributed values
    TLX_ATTRIBUTE_ALWAYS_INLINE
    std::pair<double, double> next_two_gaussians(double mean, double stdev) {
        double U = next(), V = next();
        double a = stdev * std::sqrt(-2*std::log(U));
        double b = 2 * M_PI * V;

        U = mean + a * std::cos(b);
        V = mean + a * std::sin(b);
        return std::make_pair(U, V);
    }

    //! Generate `size` uniform integers from [min, max] (i.e., both inclusive)
    template <typename int_t>
    void generate_int_block(int_t min, int_t max, int_t *arr, size_t size)
    {
        for (size_t i = 0; i < size; ++i) {
            arr[i] = next_int(min, max);
        }
    }

    //! Generate `size` uniform integers from [min, max] (i.e., both inclusive)
    template <typename int_t>
    void generate_int_block(int_t min, int_t max, std::vector<int_t> &output,
                            size_t size)
    {
        if (size > output.size()) {
            output.resize(size);
        }
        generate_int_block(min, max, output.data(), size);
    }

    //! Generate `size` geometrically integers with parameter p
    template <typename int_t>
    void generate_geometric_block(double p, int_t* arr,
                                  size_t size)
    {
        const double denominator = std::log(1.0 - p);
        for (size_t i = 0; i < size; ++i) {
            arr[i] = std::log(next()) / denominator;
        }
    }

    //! Generate `size` geometrically integers with parameter p
    template <typename int_t>
    void generate_geometric_block(double p, std::vector<int_t> &output,
                                  size_t size)
    {
        if (size > output.size()) {
            output.resize(size);
        }
        generate_geometric_block(p, output.data(), size);
    }

    //! Generate `size` exponentially distributed integers with rate `lambda`
    //! and displacement `displacement`
    void generate_exponential_block(double lambda, double *arr, size_t size) {
        generate_block(arr, size);

        for (size_t i = 0; i < size; ++i) {
            arr[i] = -std::log(arr[i]) / lambda;
        }
    }

    //! Generate `size` exponentially distributed integers with rate `lambda`
    //! and displacement `displacement`
    void generate_exponential_block(double lambda, std::vector<double> &output,
                                    size_t size) {
        if (size > output.size()) {
            output.resize(size);
        }
        generate_exponential_block(lambda, output.data(), size);
    }

    //! Generate `size` normally distributed integers with mean `mean` and
    //! standard deviation `stdev` using the Box-Muller (2) method.  Faster
    //! methods exist (e.g. g++'s std::normal_distribution)
    void generate_gaussian_block(double mean, double stdev, double *arr,
                                 size_t size) {
        // this method generates two at a time, handle the last element
        // differently if size is odd
        const bool odd = (size % 2 == 1);
        if (odd) --size;

        // first generate uniform values
        generate_block(arr, size);

        for (size_t i = 0; i < size; i += 2) {
            double U = arr[i], V = arr[i+1];
            double a = stdev * std::sqrt(-2*std::log(U));
            double b = 2 * M_PI * V;

            arr[i]   = mean + a * std::cos(b);
            arr[i+1] = mean + a * std::sin(b);
        }

        if (odd) {
            arr[size - 1] = next_gaussian(mean, stdev);
        }
    }

    //! Generate `size` normally distributed integers with mean `mean` and
    //! standard deviation `stdev` using the Box-Muller (2) method.  Faster
    //! methods exist (e.g. g++'s std::normal_distribution)
    void generate_gaussian_block(double mean, double stdev,
                                 std::vector<double> &output, size_t size) {
        // only even sizes are supported
        if (size % 2 == 1) ++size;

        if (size > output.size()) {
            output.resize(size);
        }
        generate_gaussian_block(mean, stdev, output.data(), size);
    }

    //! Alias for next()
    TLX_ATTRIBUTE_ALWAYS_INLINE
    double operator()() {
        return next();
    }

private:
    _dSFMT::dsfmt_t dsfmt_;
    std::vector<double> randblock_;
    size_t index_, block_size_, block_id_;
};

} // namespace generators
} // namespace rmat

#endif // RMAT_GENERATORS_DSFMT_HEADER
