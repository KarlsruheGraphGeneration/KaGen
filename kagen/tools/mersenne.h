/*******************************************************************************
 * include/tools/mersenne.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include "kagen/kagen.h"

#include <random>

namespace kagen {
class Mersenne {
public:
    Mersenne() : Mersenne(0) {}
    Mersenne(const SInt seed) : gen_(seed), dis_(0.0, 1.0) {}

    void RandomInit(SInt seed) {
        gen_.seed(seed);
    }

    void RandomInitByArray(SInt seeds[], SInt NumSeeds) {
        std::seed_seq sseq(seeds, seeds + NumSeeds);
        gen_.seed(sseq);
    }

    SInt BRandom() {
        return gen_();
    }

    double Random() {
        return dis_(gen_);
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
};
} // namespace kagen
