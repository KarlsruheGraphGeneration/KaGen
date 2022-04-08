/*******************************************************************************
 * include/tools/macros_assertions.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef ASSERT_H
#define ASSERT_H

#include <cassert>
#include <cmath>    // fabs
#include <cstdio>   // fprintf(), stderr
#include <cstdlib>  // abort()
#include <iostream> // cerr

#include "macros_common.h"

namespace kagen {

// A custom assertion macro that does not kill the program but prints to
// stderr instead.
#if (defined(NDEBUG) || defined(SPEEDPROFILING))
    #define ASSERT_TRUE(x) \
        do {               \
        } while (false);
#else
    #define ASSERT_TRUE(expression)                                                                                    \
        do {                                                                                                           \
            if (not(expression)) {                                                                                     \
                std::cerr << "ASSERTION FAILED [" << __FILE__ << ":" << __LINE__ << "]. Asserted: " << STR(expression) \
                          << std::endl;                                                                                \
                abort();                                                                                               \
            }                                                                                                          \
        } while (false)
#endif

// Assert: left != right.
//#ifdef NDEBUG
#if (defined(NDEBUG) || defined(SPEEDPROFILING))
    #define ASSERT_NEQ(left, right) \
        do {                        \
        } while (false);
#else
    #define ASSERT_NEQ(left, right)                                                                              \
        do {                                                                                                     \
            if ((left) == (right)) {                                                                             \
                std::cerr << "ASSERTION FAILED [" << __FILE__ << ":" << __LINE__ << "]. Asserted: " << STR(left) \
                          << " != " << STR(right) << " but was " << left << " == " << right << std::endl;        \
                abort();                                                                                         \
            }                                                                                                    \
        } while (false)
#endif

// Assert: left == right.
//#ifdef NDEBUG
#if (defined(NDEBUG) || defined(SPEEDPROFILING))
    #define ASSERT_EQ(left, right) \
        do {                       \
        } while (false);
#else
    #define ASSERT_EQ(left, right)                                                                               \
        do {                                                                                                     \
            if ((left) != (right)) {                                                                             \
                std::cerr << "ASSERTION FAILED [" << __FILE__ << ":" << __LINE__ << "]. Asserted: " << STR(left) \
                          << " == " << STR(right) << " but was " << left << " != " << right << std::endl;        \
                abort();                                                                                         \
            }                                                                                                    \
        } while (false)
#endif

// Assert: left < right.
//#ifdef NDEBUG
#if (defined(NDEBUG) || defined(SPEEDPROFILING))
    #define ASSERT_LT(left, right) \
        do {                       \
        } while (false);
#else
    #define ASSERT_LT(left, right)                                                                               \
        do {                                                                                                     \
            if ((left) >= (right)) {                                                                             \
                std::cerr << "ASSERTION FAILED [" << __FILE__ << ":" << __LINE__ << "]. Asserted: " << STR(left) \
                          << " < " << STR(right) << " but was " << left << " >= " << right << std::endl;         \
                abort();                                                                                         \
            }                                                                                                    \
        } while (false)
#endif

// Assert: left > right.
//#ifdef NDEBUG
#if (defined(NDEBUG) || defined(SPEEDPROFILING))
    #define ASSERT_GT(left, right) \
        do {                       \
        } while (false);
#else
    #define ASSERT_GT(left, right)                                                                               \
        do {                                                                                                     \
            if ((left) <= (right)) {                                                                             \
                std::cerr << "ASSERTION FAILED [" << __FILE__ << ":" << __LINE__ << "]. Asserted: " << STR(left) \
                          << " > " << STR(right) << " but was " << left << " <= " << right << std::endl;         \
                abort();                                                                                         \
            }                                                                                                    \
        } while (false)
#endif

// Assert: left <= right.
//#ifdef NDEBUG
#if (defined(NDEBUG) || defined(SPEEDPROFILING))
    #define ASSERT_LEQ(left, right) \
        do {                        \
        } while (false);
#else
    #define ASSERT_LEQ(left, right)                                                                              \
        do {                                                                                                     \
            if ((left) > (right)) {                                                                              \
                std::cerr << "ASSERTION FAILED [" << __FILE__ << ":" << __LINE__ << "]. Asserted: " << STR(left) \
                          << " <= " << STR(right) << " but was " << left << " > " << right << std::endl;         \
                abort();                                                                                         \
            }                                                                                                    \
        } while (false)
#endif

// Assert: left >= right.
//#ifdef NDEBUG
#if (defined(NDEBUG) || defined(SPEEDPROFILING))
    #define ASSERT_GEQ(left, right) \
        do {                        \
        } while (false);
#else
    #define ASSERT_GEQ(left, right)                                                                              \
        do {                                                                                                     \
            if ((left) < (right)) {                                                                              \
                std::cerr << "ASSERTION FAILED [" << __FILE__ << ":" << __LINE__ << "]. Asserted: " << STR(left) \
                          << " >= " << STR(right) << " but was " << left << " < " << right << std::endl;         \
                abort();                                                                                         \
            }                                                                                                    \
        } while (false)
#endif

// Assert: x <= y <= z.
//#ifdef NDEBUG
#if (defined(NDEBUG) || defined(SPEEDPROFILING))
    #define ASSERT_BETWEEN(x, y, z) \
        do {                        \
        } while (false);
#else
    #define ASSERT_BETWEEN(left, x, right)                                                                       \
        do {                                                                                                     \
            if (((left) > (x)) or ((right) < (x))) {                                                             \
                std::cerr << "ASSERTION FAILED [" << __FILE__ << ":" << __LINE__ << "]. Asserted: " << STR(x)    \
                          << " in {" << STR(left) << ", ..., " << STR(right) << "} but was " << x << " not in {" \
                          << left << ", ..., " << right << "}." << std::endl;                                    \
                abort();                                                                                         \
            }                                                                                                    \
        } while (false)
#endif

// Assert: \forall begin <= i < end: sequence[i] > x.
//#ifdef NDEBUG
#if (defined(NDEBUG) || defined(SPEEDPROFILING))
    #define ASSERT_RANGE_GT(sequence, begin, end, x, i) \
        do {                                            \
        } while (false);
#else
    #define ASSERT_RANGE_GT(sequence, begin, end, x, i) \
        for (int i = begin; i < end; ++i) {             \
            ASSERT_GT(sequence[i], x);                  \
        }
#endif

// Assert: \forall begin <= i < end: sequence[i] >= x.
//#ifdef NDEBUG
#if (defined(NDEBUG) || defined(SPEEDPROFILING))
    #define ASSERT_RANGE_GEQ(sequence, begin, end, x, i) \
        do {                                             \
        } while (false);
#else
    #define ASSERT_RANGE_GEQ(sequence, begin, end, x, i) \
        for (int i = begin; i < end; ++i) {              \
            ASSERT_GEQ(sequence[i], x);                  \
        }
#endif

//#ifdef NDEBUG
#if (defined(NDEBUG) || defined(SPEEDPROFILING))
    #define ASSERT_RANGE_EQ(sequence, begin, end, x) \
        do {                                         \
        } while (false);
#else
    #define ASSERT_RANGE_EQ(sequence, begin, end, x) \
        for (unsigned int i = begin; i < end; ++i) { \
            ASSERT_EQ(sequence[i], x);               \
        }
#endif

} // namespace kagen
#endif
