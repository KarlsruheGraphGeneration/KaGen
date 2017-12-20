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
#include "var_gen.h"

template <typename LPStocc = VarGen<LPFloat>, typename HPStocc = VarGen<>>
class RNGWrapper {
 public:
  RNGWrapper(const PGeneratorConfig &config)
      : config_(config),
        lp_gen_(0),
        hp_gen_(0) {};

  HPFloat GenerateHypergeometric(SInt seed, HPFloat n, HPFloat m, HPFloat N) {
    hp_gen_.RandomInit(seed);
    HPFloat variate = 0;
    if (config_.use_binom)
      variate = hp_gen_.Binomial(n, (double)(m / N));
    else
      variate = hp_gen_.Hypergeometric(n, m, N);
    return variate;
  }

  HPFloat GenerateBinomial(SInt seed, HPFloat n, double p) {
    hp_gen_.RandomInit(seed);
    HPFloat variate = hp_gen_.Binomial(n, p);
    return variate;
  }

  SInt GenerateBinomial(SInt seed, SInt n, double p) {
    lp_gen_.RandomInit(seed);
    HPFloat variate = lp_gen_.Binomial(n, p);
    return variate;
  }

  template <typename F>
  void GenerateSample(SInt seed, HPFloat N, SInt n, F &&callback) {
    sampling::HashSampling<> hs(seed, config_.base_size);
    sampling::SeqDivideSampling<> sds(hs, config_.base_size, seed, config_.use_binom);
    sds.sample(N, n, callback);
  }

 private:
  const PGeneratorConfig &config_;

  LPStocc lp_gen_;
  HPStocc hp_gen_;
};

#endif
