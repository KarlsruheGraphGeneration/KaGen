/*******************************************************************************
 * include/generator_config.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _GENERATOR_CONFIG_H_
#define _GENERATOR_CONFIG_H_

#include <string>
#include "definitions.h"

namespace kagen {

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

}
#endif
