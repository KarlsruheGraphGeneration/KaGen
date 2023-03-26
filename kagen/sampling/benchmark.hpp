/*******************************************************************************
 * sampling/benchmark.hpp
 *
 * Copyright (C) 2016-2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef SAMPLING_BENCHMARK_HEADER
#define SAMPLING_BENCHMARK_HEADER

#include <cmath>

namespace sampling {

// template-based loop unrolling
template <std::size_t N> struct faux_unroll {
    template <typename F> static void call(F &&f) {
        faux_unroll<N-1>::call(f);
        f(N-1);
    }
};

template <> struct faux_unroll<0> {
    template <typename F> static void call(F &&) {}
};


struct statistics {
    // Single-pass standard deviation calculation as described in Donald Knuth:
    // The Art of Computer Programming, Volume 2, Chapter 4.2.2, Equations 15&16
    double mean;
    double nvar; // approx n * variance; stddev = sqrt(nvar / (count-1))
    std::size_t count;

    statistics() : mean(0.0), nvar(0.0), count(0) {}

    void push(double t) {
        ++count;
        if (count == 1) {
            mean = t;
        } else {
            double oldmean = mean;
            mean += (t - oldmean) / count;
            nvar += (t - oldmean) * (t - mean);
        }
    }

    double avg() {
        return mean;
    }
    double stddev(std::size_t ddof = 1) {
        if (count <= 1) return 0.0;
        // ddof = delta degrees of freedom
        // Set to 0 if you have the entire distribution
        // Set to 1 if you have a sample (to correct for bias)
        return sqrt(nvar / (count - ddof));
    }
};

}

#endif // SAMPLING_BENCHMARK_HEADER

