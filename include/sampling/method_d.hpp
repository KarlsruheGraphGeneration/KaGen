/*******************************************************************************
 * include/methodD.hpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef METHOD_D_HEADER
#define METHOD_D_HEADER

#include "definitions.hpp"
#include "dSFMT.hpp"

#include <cmath>

namespace sampling {

class Vitter {
    public:
        Vitter(ULONG seed)
            : gen(seed)
        { }

        // Sampling method A from Vitter et al.
        //
        // \param N Size of population.
        // \param n Number of samples.
        // \param gen Uniform random variate generator.
        // \param callback Function to process sample.
        //
        template <typename F>
        void methodA(ULONG N,
                     ULONG n,
                     F &&callback) {
            assert(n <= N);

            // Initialization
            ULONG sample = 0;
            double Nreal = (double) N;
            double top = Nreal - n;

            // Main loop
            while (n >= 2) {
                ULONG S = 0;
                double V = gen();
                double quot = top / Nreal;
                while (quot > V) {
                    S++;
                    top -= 1.0;
                    Nreal -= 1.0;
                    quot = (quot * top) / Nreal;
                }
                // Skip over next S records and select the following one
                sample += S + 1;
                // samples[num_sample++] = sample;
                callback(sample);
                Nreal -= 1.0;
                n--;
            }
            if (n == 1) {
                ULONG S = round(Nreal) * gen();
                sample += S + 1;
                // samples[num_sample++] = sample;
                callback(sample);
            }
        }

        // Sampling method D from Vitter et al.
        //
        // \param N Size of population.
        // \param n Number of samples.
        // \param gen Uniform random variate generator.
        // \param samples Function to process sample.
        //
        template <typename F>
        void sample(ULONG N,
                    ULONG n,
                    F &&callback) {
            assert(n <= N);

            ULONG initialN = N;
            // Initialization
            ULONG sample = 0;
            // ULONG num_sample = 0;
            double nreal = (double) n;
            double ninv = 1.0 / nreal;
            double Nreal = (double) N;
            double Vprime = exp(log(gen()) * ninv);
            ULONG qu1 = N + 1 - n;
            double qu1real = Nreal + 1.0 - nreal;
            ULONG negalphainv = -13;
            ULONG threshold = n * (-negalphainv);
            ULONG S = 0;

            // Main loop
            while (n > 1 && threshold < N) {
                double nmin1inv = 1.0 / (nreal - 1.0);
                double negSreal = 0.0;

                while (true) {
                    // Step D2: Generate U and X
                    double X;
                    while (true) {
                        X = Nreal * (1.0 - Vprime);
                        S = X;
                        if (S < qu1) break;
                        Vprime = exp(log(gen()) * ninv);
                    }

                    double U = gen();
                    negSreal = -(double)S;

                    // Step D3: Accept?
                    double y1 = exp(log(U * Nreal / qu1real) * nmin1inv);
                    Vprime = y1 * (-X / Nreal + 1.0) * (qu1real / (negSreal + qu1real));
                    if (Vprime <= 1.0) break; // Accept!

                    // Step D4: Accept?
                    double y2 = 1.0; double top = Nreal - 1.0;
                    double bottom;
                    double limit;
                    if (n - 1 > S) {
                        bottom = Nreal - nreal;
                        limit = N - S;
                    } else {
                        bottom = negSreal + Nreal - 1.0;
                        limit = qu1;
                    }

                    for (ULONG t = N; t > limit; t--) {
                        y2 = (y2 * top) / bottom;
                        top -= 1.0;
                        bottom -= 1.0;
                    }

                    if (Nreal / (Nreal - X) >= y1 * exp(log(y2) * nmin1inv)) {
                        // Accept!
                        Vprime = exp(log(gen()) * nmin1inv);
                        break;
                    }
                    Vprime = exp(log(gen()) * ninv);
                }
                // Skip over next S records and select the following one
                sample += S + 1;
                // samples[num_sample++] = sample;
                callback(sample);
                N = (N - 1) - S;
                Nreal = (Nreal - 1.0) + negSreal;
                n--;
                nreal -= 1.0;
                ninv = nmin1inv;
                qu1 -= S;
                qu1real += negSreal;
                threshold += negalphainv;
            }

            if (n > 1) {
                ULONG currentN = N;
                methodA(N, n, [&](ULONG sample) {
                        callback(sample + initialN - currentN);
                    });
            } else if (n == 1) {
                S = N * Vprime;
                // Skip over next S records and select the following one
                sample += S + 1;
                // samples[num_sample++] = sample;
                callback(sample);
            }
        }

    private:
        dSFMT gen;
};

} // namespace sampling

#endif // METHOD_D_HEADER
