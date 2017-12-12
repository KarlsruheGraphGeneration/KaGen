/******************************************************************************
 * mersenne.h
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

#ifndef _MERSENNE_H_
#define _MERSENNE_H_

#include <random>

void EndOfProgram(void);  // System-specific exit code (userintf.cpp)

void FatalError(
    const char *ErrorText);  // System-specific error reporting (userintf.cpp)

class Mersenne {
 public:
  Mersenne() { Mersenne(0); };

  Mersenne(SInt seed) : gen_(seed), dis_(0.0, 1.0){};

  void RandomInit(SInt seed) { gen_.seed(seed); }

  void RandomInitByArray(SInt seeds[], SInt NumSeeds) {
    std::seed_seq sseq(seeds, seeds + NumSeeds);
    gen_.seed(sseq);
  }

  SInt BRandom() { return gen_(); }

  double Random() { return dis_(gen_); }

  SInt IRandom(SInt min, SInt max) {
    if (max == min) return min;
    SInt r = (SInt)((double)(SInt)(max - min + 1) * Random() + min);
    if (r > max) r = max;
    return r;
  }

 private:
  std::mt19937_64 gen_;
  std::uniform_real_distribution<double> dis_;
};

#endif
