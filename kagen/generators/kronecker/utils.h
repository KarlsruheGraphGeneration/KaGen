/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef UTILS_H
#define UTILS_H

#ifndef __STDC_CONSTANT_MACROS
    #define __STDC_CONSTANT_MACROS
#endif
#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef GRAPH_GENERATOR_MPI
    #include <mpi.h>
#endif
#ifdef GRAPH_GENERATOR_OMP
    #include <omp.h>
#endif

#include "kagen/generators/kronecker/splittable_mrg.h"

#if defined(_OPENMP)
    #define OMP(x_) _Pragma(x_)
#else
    #define OMP(x_)
#endif

#ifdef __cplusplus
extern "C" {
#endif

// void* xrealloc(void* p, size_t nbytes); /* In utils.c */
// uint_fast64_t random_up_to(mrg_state* st, uint_fast64_t n);

void make_mrg_seed(uint64_t userseed1, uint64_t userseed2, uint_fast32_t* seed) {
    seed[0] = (uint32_t)(userseed1 & UINT32_C(0x3FFFFFFF)) + 1;
    seed[1] = (uint32_t)((userseed1 >> 30) & UINT32_C(0x3FFFFFFF)) + 1;
    seed[2] = (uint32_t)(userseed2 & UINT32_C(0x3FFFFFFF)) + 1;
    seed[3] = (uint32_t)((userseed2 >> 30) & UINT32_C(0x3FFFFFFF)) + 1;
    seed[4] = (uint32_t)((userseed2 >> 60) << 4) + (uint32_t)(userseed1 >> 60) + 1;
}

#ifdef __cplusplus
}
#endif

#endif /* UTILS_H */
