/******************************************************************************
 * rng_wrapper.h
 *
 * Source of the sampling routine
 ******************************************************************************
 * Copyright (C) 2016 Sebastian Lamm <lamm@ira.uka.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

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
