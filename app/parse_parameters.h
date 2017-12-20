/******************************************************************************
 * parse_parameters.h
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

#ifndef _PARSE_PARAMETERS_H_
#define _PARSE_PARAMETERS_H_

#include <string.h>

#include "generator_config.h"
#include "tools/arg_parser.h"

#include "definitions.h"

void ParseParameters(int argn, char **argv,
                     PGeneratorConfig &generator_config) {
  ArgParser args(argn, argv);

  // Generator
  generator_config.generator = args.Get<std::string>("gen", "");

  // Nodes
  bool exact_n = args.IsSet("exact_n");
  if (exact_n)
    generator_config.n = args.Get<ULONG>("n", 100);
  else
    generator_config.n = (ULONG)1 << args.Get<ULONG>("n", 3);

  // Blocks
  generator_config.k = args.Get<ULONG>("k", 4);

  // RNG
  generator_config.seed = args.Get<ULONG>("seed", 1);
  generator_config.hash_sample = args.Get<bool>("hash_sample", false);
  generator_config.use_binom = args.IsSet("binom");

  // I/O
  generator_config.output_file = args.Get<std::string>("output", "out");
  generator_config.debug_output = args.Get<std::string>("debug", "dbg");
  generator_config.dist_size = args.Get<ULONG>("dist", 10);

  // Edges
  bool exact_m = args.IsSet("exact_m");
  if (exact_m)
    generator_config.m = args.Get<ULONG>("m", 0);
  else
    generator_config.m = (ULONG)1 << args.Get<ULONG>("m", 0);
  generator_config.p = args.Get<double>("p", 0.0);
  generator_config.self_loops = args.Get<bool>("self_loops", false);

  // Radius/Edges
  generator_config.r = args.Get<double>("r", 0.125);

  // Average degree
  generator_config.avg_degree = args.Get<double>("d", 5.0);
  generator_config.plexp = args.Get<double>("gamma", 2.6);

  // RHG 
  generator_config.thres = args.Get<ULONG>("t", 0);
  generator_config.query_both = args.Get<bool>("qb", false);

  // BA
  generator_config.min_degree = args.Get<ULONG>("md", 4);

  // Floating-point precision
  generator_config.precision = args.Get<ULONG>("prec", 32);

  // Sampling algorithm
  generator_config.base_size = (ULONG)1 << args.Get<ULONG>("sk", 8);
  generator_config.hyp_base = (ULONG)1 << args.Get<ULONG>("hk", 8);

  // Benchmarks
  generator_config.iterations = args.Get<ULONG>("i", 1);
}

#endif
