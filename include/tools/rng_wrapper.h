/*******************************************************************************
 * include/tools/rng_wrapper.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _RNG_WRAPPER_H_
#define _RNG_WRAPPER_H_

#include <random>

#include "methodR.hpp"

namespace kagen {

class RNGWrapper {
 public:
  RNGWrapper(const PGeneratorConfig &config)
      : config_(config),
        rng_(0),
        hyp_(0) {};

  SInt GenerateHypergeometric(SInt seed, SInt n, SInt m, SInt N) {
    SInt variate = 0;
    if (config_.use_binom)
      variate = GenerateBinomial(seed, n, (LPFloat)(m / N));
    else
      hyp_.seed(seed);
      if (m < 1) return 0;
      variate = hyp_(n, N-n, m);
    return variate;
  }

  SInt GenerateBinomial(SInt seed, SInt n, LPFloat p) {
    rng_.seed(seed);
    std::binomial_distribution<SInt> bin(n, p);
    return bin(rng_);
  }

  template <typename F>
  void GenerateSample(SInt seed, SInt N, SInt n, F &&callback) {
    sampling::HashSampling<> hs(seed, config_.base_size);
    sampling::SeqDivideSampling<> sds(hs, config_.base_size, seed, config_.use_binom);
    sds.sample(N, n, callback);
  }

 private:
  const PGeneratorConfig &config_;

  std::mt19937 rng_;
  sampling::hypergeometric_distribution<> hyp_;
};

}
#endif
