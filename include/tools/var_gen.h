/*******************************************************************************
 * include/tools/var_gen.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _VAR_GEN_H_
#define _VAR_GEN_H_

#include <cstddef>

#include "definitions.h"

template <typename Float = HPFloat>
class VarGen {
 public:
  VarGen(SInt seed) : rng_(seed), hyp_(seed) { }

  void RandomInit(SInt seed) {
    rng_.seed(seed);
    hyp_.seed(seed);
  }

  Float Binomial(Float n, double p) {
    std::binomial_distribution<> bin(n, p);
    return bin(rng_);
  }

  Float Hypergeometric(Float n, Float m, Float N) {
    if (m < 1) return 0;
    return hyp_(n, N-n, m);
  }

 private:
  std::mt19937 rng_;
  sampling::hypergeometric_distribution<> hyp_;
};

#endif
