/******************************************************************************
 * sorted_mersenne.h
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

#ifndef _SORTED_MERSENNE_H_
#define _SORTED_MERSENNE_H_

#include <random>

void EndOfProgram(void);  // System-specific exit code (userintf.cpp)

void FatalError(
    const char *ErrorText);  // System-specific error reporting (userintf.cpp)

class SortedMersenne {
 public:
  SortedMersenne() { SortedMersenne(0); };

  SortedMersenne(SInt seed)
      : gen_(seed), dis_(0.0, 1.0), num_samples_(100), ln_cur_max_(0.0){};

  void RandomInit(SInt seed, SInt samples) {
    gen_.seed(seed);
    num_samples_ = samples;
    ln_cur_max_ = 0.0;
  }

  void RandomInitByArray(SInt seeds[], SInt NumSeeds) {
    std::seed_seq sseq(seeds, seeds + NumSeeds);
    gen_.seed(sseq);
  }

  SInt BRandom() { return gen_(); }

  double Random() {
    double rand = dis_(gen_);
    ln_cur_max_ += std::log(rand) / (double)num_samples_;
    num_samples_--;
    return std::exp(ln_cur_max_);
  }

  SInt IRandom(SInt min, SInt max) {
    if (max == min) return min;
    SInt r = (SInt)((double)(SInt)(max - min + 1) * Random() + min);
    if (r > max) r = max;
    return r;
  }

 private:
  std::mt19937_64 gen_;
  std::uniform_real_distribution<double> dis_;

  SInt num_samples_;
  double ln_cur_max_;
};

#endif
