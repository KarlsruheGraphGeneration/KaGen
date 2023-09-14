/*******************************************************************************
 * rmat/generators/dSFMT.cpp
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

#include "dSFMT.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace rmat {
namespace generators {

// set the name
const char* dSFMT::name = "dSFMT";

namespace _dSFMT {

/*----------------
  STATIC FUNCTIONS
  ----------------*/
inline static void gen_rand_array_c0o1(dsfmt_t* dsfmt, w128_t* array, int size);
inline static int  idxof(int i);
static void        initial_mask(dsfmt_t* dsfmt);
static void        period_certification(dsfmt_t* dsfmt);

#if defined(SAMPLING_HAVE_SSE2)
/** 1 in 64bit for sse2 */
static const union X128I_T sse2_int_one = {{1, 1}};
/** 2.0 double for sse2 */
static const union X128D_T sse2_double_two = {{2.0, 2.0}};
/** -1.0 double for sse2 */
static const union X128D_T sse2_double_m_one = {{-1.0, -1.0}};
#endif

/**
 * This function simulate a 32-bit array index overlapped to 64-bit
 * array of LITTLE ENDIAN in BIG ENDIAN machine.
 */
#if defined(DSFMT_BIG_ENDIAN)
inline static int idxof(int i) {
    return i ^ 1;
}
#else
inline static int idxof(int i) {
    return i;
}
#endif

#if defined(SAMPLING_HAVE_SSE2)
/**
 * This function converts the double precision floating point numbers which
 * distribute uniformly in the range [1, 2) to those which distribute uniformly
 * in the range [0, 1).
 * @param w 128bit stracture of double precision floating point numbers (I/O)
 */
inline static void convert_c0o1(w128_t* w) {
    w->sd = _mm_add_pd(w->sd, sse2_double_m_one.d128);
}
#else /* standard C and altivec */
/**
 * This function converts the double precision floating point numbers which
 * distribute uniformly in the range [1, 2) to those which distribute uniformly
 * in the range [0, 1).
 * @param w 128bit stracture of double precision floating point numbers (I/O)
 */
inline static void convert_c0o1(w128_t* w) {
    w->d[0] -= 1.0;
    w->d[1] -= 1.0;
}
#endif

/**
 * This function fills the user-specified array with double precision
 * floating point pseudorandom numbers of the IEEE 754 format.
 * @param dsfmt dsfmt state vector.
 * @param array an 128-bit array to be filled by pseudorandom numbers.
 * @param size number of 128-bit pseudorandom numbers to be generated.
 */
inline static void gen_rand_array_c0o1(dsfmt_t* dsfmt, w128_t* array, int size) {
    int    i, j;
    w128_t lung;

    lung = dsfmt->status[DSFMT_N];
    do_recursion(&array[0], &dsfmt->status[0], &dsfmt->status[DSFMT_POS1], &lung);
    for (i = 1; i < DSFMT_N - DSFMT_POS1; i++) {
        do_recursion(&array[i], &dsfmt->status[i], &dsfmt->status[i + DSFMT_POS1], &lung);
    }
    for (; i < DSFMT_N; i++) {
        do_recursion(&array[i], &dsfmt->status[i], &array[i + DSFMT_POS1 - DSFMT_N], &lung);
    }
    for (; i < size - DSFMT_N; i++) {
        do_recursion(&array[i], &array[i - DSFMT_N], &array[i + DSFMT_POS1 - DSFMT_N], &lung);
        convert_c0o1(&array[i - DSFMT_N]);
    }
    for (j = 0; j < 2 * DSFMT_N - size; j++) {
        dsfmt->status[j] = array[j + size - DSFMT_N];
    }
    for (; i < size; i++, j++) {
        do_recursion(&array[i], &array[i - DSFMT_N], &array[i + DSFMT_POS1 - DSFMT_N], &lung);
        dsfmt->status[j] = array[i];
        convert_c0o1(&array[i - DSFMT_N]);
    }
    for (i = size - DSFMT_N; i < size; i++) {
        convert_c0o1(&array[i]);
    }
    dsfmt->status[DSFMT_N] = lung;
}

/**
 * This function initializes the internal state array to fit the IEEE
 * 754 format.
 * @param dsfmt dsfmt state vector.
 */
static void initial_mask(dsfmt_t* dsfmt) {
    int       i;
    uint64_t* psfmt;

    psfmt = &dsfmt->status[0].u[0];
    for (i = 0; i < DSFMT_N * 2; i++) {
        psfmt[i] = (psfmt[i] & DSFMT_LOW_MASK) | DSFMT_HIGH_CONST;
    }
}

/**
 * This function certificate the period of 2^{SFMT_MEXP}-1.
 * @param dsfmt dsfmt state vector.
 */
static void period_certification(dsfmt_t* dsfmt) {
    uint64_t pcv[2] = {DSFMT_PCV1, DSFMT_PCV2};
    uint64_t tmp[2];
    uint64_t inner;
    int      i;
#if (DSFMT_PCV2 & 1) != 1
    int      j;
    uint64_t work;
#endif

    tmp[0] = (dsfmt->status[DSFMT_N].u[0] ^ DSFMT_FIX1);
    tmp[1] = (dsfmt->status[DSFMT_N].u[1] ^ DSFMT_FIX2);

    inner = tmp[0] & pcv[0];
    inner ^= tmp[1] & pcv[1];
    for (i = 32; i > 0; i >>= 1) {
        inner ^= inner >> i;
    }
    inner &= 1;
    /* check OK */
    if (inner == 1) {
        return;
    }
    /* check NG, and modification */
#if (DSFMT_PCV2 & 1) == 1
    dsfmt->status[DSFMT_N].u[1] ^= 1;
#else
    for (i = 1; i >= 0; i--) {
        work = 1;
        for (j = 0; j < 64; j++) {
            if ((work & pcv[i]) != 0) {
                dsfmt->status[DSFMT_N].u[i] ^= work;
                return;
            }
            work = work << 1;
        }
    }
#endif
    return;
}

/*----------------
  PUBLIC FUNCTIONS
  ----------------*/
/**
 * This function returns the minimum size of array used for \b
 * fill_array functions.
 * @return minimum size of array used for fill_array functions.
 */
int dsfmt_get_min_array_size(void) {
    return DSFMT_N64;
}

/**
 * This function generates double precision floating point
 * pseudorandom numbers which distribute in the range [0, 1) to the
 * specified array[] by one call. This function is the same as
 * fill_array_close1_open2() except the distribution range.
 *
 * @param array an array where pseudorandom numbers are filled
 * by this function.
 * @param dsfmt dsfmt state vector.
 * @param size the number of pseudorandom numbers to be generated.
 * see also \sa fill_array_close1_open2()
 */
void dsfmt_fill_array_close_open(dsfmt_t* dsfmt, double array[], int size) {
    assert(size % 2 == 0);
    assert(size >= DSFMT_N64);
    gen_rand_array_c0o1(dsfmt, (w128_t*)array, size / 2);
}

#if defined(__INTEL_COMPILER)
    #pragma warning(disable : 981)
#endif
/**
 * This function initializes the internal state array with a 32-bit
 * integer seed.
 * @param dsfmt dsfmt state vector.
 * @param seed a 32-bit integer used as the seed.
 * @param mexp caller's mersenne expornent
 */
void dsfmt_chk_init_gen_rand(dsfmt_t* dsfmt, uint32_t seed, int mexp) {
    int       i;
    uint32_t* psfmt;

    /* make sure caller program is compiled with the same MEXP */
    if (mexp != DSFMT_MEXP) {
        fprintf(stderr, "DSFMT_MEXP doesn't match with dSFMT.c\n");
        exit(1);
    }
    psfmt           = &dsfmt->status[0].u32[0];
    psfmt[idxof(0)] = seed;
    for (i = 1; i < (DSFMT_N + 1) * 4; i++) {
        psfmt[idxof(i)] = 1812433253UL * (psfmt[idxof(i - 1)] ^ (psfmt[idxof(i - 1)] >> 30)) + i;
    }
    initial_mask(dsfmt);
    period_certification(dsfmt);
    dsfmt->idx = DSFMT_N64;
}
#if defined(__INTEL_COMPILER)
    #pragma warning(default : 981)
#endif

} // namespace _dSFMT
} // namespace generators
} // namespace rmat
