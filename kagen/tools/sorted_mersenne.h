/*******************************************************************************
 * include/tools/sorted_mersenne.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include "kagen/kagen.h"

#include <random>

namespace kagen {
class SortedMersenne {
public:
    SortedMersenne() : SortedMersenne(0) {}
    SortedMersenne(SInt seed) : gen_(seed), dis_(0.0, 1.0), num_samples_(100), ln_cur_max_(0.0) {}

    void RandomInit(SInt seed, SInt samples) {
        gen_.seed(seed);
        num_samples_ = samples;
        ln_cur_max_  = 0.0;
    }

    void RandomInitByArray(SInt seeds[], SInt NumSeeds) {
        std::seed_seq sseq(seeds, seeds + NumSeeds);
        gen_.seed(sseq);
    }

    SInt BRandom() {
        return gen_();
    }

    double Random() {
        double rand = dis_(gen_);
        ln_cur_max_ += std::log(rand) / (double)num_samples_;
        num_samples_--;
        return std::exp(ln_cur_max_);
    }

    SInt IRandom(SInt min, SInt max) {
        if (max == min)
            return min;
        SInt r = (SInt)((double)(SInt)(max - min + 1) * Random() + min);
        if (r > max)
            r = max;
        return r;
    }

private:
    std::mt19937_64                        gen_;
    std::uniform_real_distribution<double> dis_;

    SInt   num_samples_;
    double ln_cur_max_;
};
} // namespace kagen
