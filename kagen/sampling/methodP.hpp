/*******************************************************************************
 * sampling/methodP.hpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef SAMPLING_METHOD_P_HEADER
#define SAMPLING_METHOD_P_HEADER

#include "kagen/sampling/definitions.hpp"
#include "kagen/sampling/hash.hpp"
#include "kagen/sampling/methodR.hpp"
#include "kagen/sampling/sampling_config.hpp"

#include <random>
#include <vector>

namespace sampling {

template <typename LocalSampler = SeqDivideSampling<>, typename H = CRCHash /*Spooky*/>
class ParDivideSampling {
public:
    ParDivideSampling(ULONG universe_size, ULONG base_size, ULONG seed, PEID size)
        : hyp(seed),
          hash_seed(seed),
          base_size(base_size) {
        // Compute input distribution
        rem = universe_size % size;
        div = universe_size / size;
    }

    template <typename F>
    void sample(ULONG n, PEID j, PEID k, PEID i, F&& callback, ULONG offset = 0) {
        if (j - k == 0) {
            ULONG                            h = H::hash(hash_seed + i);
            typename LocalSampler::base_type base_sampler(h, base_size);
            // Allocate hash table for base case

            LocalSampler local_sampler(base_sampler, base_size, h);
            local_sampler.sample(N(i + 1) - N(i), n, callback, offset);
            return;
        }

        PEID  m = (j + k) / 2;
        ULONG h = H::hash(hash_seed + j + k);
        // stocc.RandomInit(h);
        ULONG N_split = N(m) - N(j - 1);
        ULONG x;
        if (use_binomial) {
            std::mt19937                 rng(h);
            std::binomial_distribution<> binom(n, (double)N_split / (N(k) - N(j - 1)));
            x = binom(rng);
            // x = stocc.Binomial(n, (double)N_split/(N(k) - N(j-1)));
        } else {
            hyp.seed(h);
            // ULONG x = stocc.Hypergeometric(N_split, n, N(k) - N(j-1));
            x = hyp(n, N(k) - N(j - 1) - n, N_split);
        }
        if (i < m)
            sample(x, j, m, i, callback, offset);
        else
            sample(n - x, m + 1, k, i, callback, offset + N_split);
    }

private:
    hypergeometric_distribution<ULONG> hyp;
    ULONG                              hash_seed;
    ULONG                              div;
    PEID                               rem;
    ULONG                              base_size;
    bool                               use_binomial;

    inline ULONG N(PEID i) {
        return i * div + std::min(i, rem);
    }
};

} // namespace sampling

#endif // SAMPLING_METHOD_P_HEADER
