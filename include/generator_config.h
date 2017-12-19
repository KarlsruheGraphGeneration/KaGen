/******************************************************************************
 * generator_config.h
 *
 * Source of the graph generator
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

#ifndef _GENERATOR_CONFIG_H_
#define _GENERATOR_CONFIG_H_

#include <string>
#include "definitions.h"

// Configuration for the generator.
struct PGeneratorConfig {
  PGeneratorConfig() {}

  // Type of generator
  std::string generator;
  // Seed for the PRNG
  int seed;
  // Number of nodes/edges
  ULONG n, m;
  // Chunk size
  ULONG k;
  // Edge probability
  double p;
  // Edge radius
  double r;
  // Output filename
  std::string output_file;
  // Debug output
  std::string debug_output;
  // Use hash tryagain sampling
  bool hash_sample;
  // Allow self loops
  bool self_loops;
  // Power-law exponent
  double plexp;
  // Avg. degree
  double avg_degree;
  // RHG clique threshold
  double thres;
  // RHG query strategy
  double query_both;
  // BA minimum degree
  double min_degree;
  // Size of histogramm
  ULONG dist_size;
  // Use binomial approximation to hypergeometric
  bool use_binom;
  // Floating-point precision
  ULONG precision;
  // Sampler base size
  ULONG base_size;
  ULONG hyp_base;
  // Benchmarks
  ULONG iterations;
};

#endif
