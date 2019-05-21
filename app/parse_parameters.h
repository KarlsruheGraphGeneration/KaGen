/*******************************************************************************
 * app/parse_parameters.cpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _PARSE_PARAMETERS_H_
#define _PARSE_PARAMETERS_H_

#include <string.h>

#include "generator_config.h"
#include "tools/arg_parser.h"

#include "definitions.h"

namespace kagen {

void ParseParameters(int argn, char **argv,
                     PEID, PEID size,
                     PGeneratorConfig &generator_config) {
  ArgParser args(argn, argv);

  // Generator
  generator_config.generator = args.Get<std::string>("gen", "");

  if (args.IsSet("help") || argn < 2) {
    if (generator_config.generator == "") {
      std::cout << "================================================" << std::endl;
      std::cout << "==================== KaGen =====================" << std::endl;
      std::cout << "================================================" << std::endl;
      std::cout << "Usage:\t\t\tmpirun -n <num_proc> ./kagen -gen <generator> [additional parameters]" << std::endl;
      std::cout << "Generators:\t\tgnm_directed|gnm_undirected|gnp_directed|gnp_undirected|rgg_2d|rgg_3d|rdg_2d|rdg_3d|ba|rhg" << std::endl;
      std::cout << "Additional help:\t./kagen -gen <generator> -help" << std::endl;
    }
    
    if (generator_config.generator == "gnm_undirected" || generator_config.generator == "gnm_directed") {
      std::cout << "================================================" << std::endl;
      std::cout << "========== Erdos-Renyi Graphs G(n,m) ===========" << std::endl;
      std::cout << "================================================" << std::endl;
      std::cout << "Parameters:" << std::endl;
      std::cout << "-n\t\t<number of vertices as a power of two>" << std::endl;
      std::cout << "-m\t\t<number of edges as a power of two>" << std::endl;
      std::cout << "-k\t\t<number of chunks>" << std::endl;
      std::cout << "-seed\t\t<seed for PRNGs>" << std::endl;
      std::cout << "-output\t\t<output file>" << std::endl;
      std::cout << "-self_loops" << std::endl;
      std::cout << "\nExample:" << std::endl;
      std::cout << "mpirun -n 16 ./build/app/kagen -gen gnm_directed -n 20 -m 22 -self_loops -output tmp" << std::endl;
    } else if (generator_config.generator == "gnp_undirected" || generator_config.generator == "gnp_directed") {
      std::cout << "================================================" << std::endl;
      std::cout << "========== Erdos-Renyi Graphs G(n,p) ===========" << std::endl;
      std::cout << "================================================" << std::endl;
      std::cout << "Parameters:" << std::endl;
      std::cout << "-n\t\t<number of vertices as a power of two>" << std::endl;
      std::cout << "-p\t\t<edge probability>" << std::endl;
      std::cout << "-k\t\t<number of chunks>" << std::endl;
      std::cout << "-seed\t\t<seed for PRNGs>" << std::endl;
      std::cout << "-output\t\t<output file>" << std::endl;
      std::cout << "-self_loops" << std::endl;
      std::cout << "\nExample:" << std::endl;
      std::cout << "mpirun -n 16 ./build/app/kagen -gen gnp_directed -n 20 -p 0.001 -self_loops -output tmp" << std::endl;
    } else if (generator_config.generator == "rgg_2d" || generator_config.generator == "rgg_3d") {
      std::cout << "================================================" << std::endl;
      std::cout << "======= Random Geometric Graphs RGG(n,d) ========" << std::endl;
      std::cout << "================================================" << std::endl;
      std::cout << "Parameters:" << std::endl;
      std::cout << "-n\t\t<number of vertices as a power of two>" << std::endl;
      std::cout << "-r\t\t<radius for vertices to be connected> (r <= 1.0)>" << std::endl;
      std::cout << "-k\t\t<number of chunks>" << std::endl;
      std::cout << "-seed\t\t<seed for PRNGs>" << std::endl;
      std::cout << "-output\t\t<output file>" << std::endl;
      std::cout << "\nExample:" << std::endl;
      std::cout << "mpirun -n 16 ./build/app/kagen -gen rgg_3d -n 20 -r 0.00275 -output tmp" << std::endl;
    } else if (generator_config.generator == "rdg_2d" || generator_config.generator == "rdg_3d") {
      std::cout << "================================================" << std::endl;
      std::cout << "======== Random Delaunay Graphs RDG(n) =========" << std::endl;
      std::cout << "================================================" << std::endl;
      std::cout << "Parameters:" << std::endl;
      std::cout << "-n\t\t<number of vertices as a power of two>" << std::endl;
      std::cout << "-k\t\t<number of chunks>" << std::endl;
      std::cout << "-seed\t\t<seed for PRNGs>" << std::endl;
      std::cout << "-output\t\t<output file>" << std::endl;
      std::cout << "\nExample:" << std::endl;
      std::cout << "mpirun -n 16 ./build/app/kagen -gen rdg_3d -n 20 -output tmp" << std::endl;
    } else if (generator_config.generator == "ba") {
      std::cout << "================================================" << std::endl;
      std::cout << "======= Barabassi-Albert Graphs BA(n,d) ========" << std::endl;
      std::cout << "================================================" << std::endl;
      std::cout << "Parameters:" << std::endl;
      std::cout << "-n\t\t<number of vertices as a power of two>" << std::endl;
      std::cout << "-md\t\t<min degree for each vertex>" << std::endl;
      std::cout << "-k\t\t<number of chunks>" << std::endl;
      std::cout << "-seed\t\t<seed for PRNGs>" << std::endl;
      std::cout << "-output\t\t<output file>" << std::endl;
      std::cout << "\nExample:" << std::endl;
      std::cout << "mpirun -n 16 ./build/app/kagen -gen ba -n 20 -md 4 -output tmp" << std::endl;
    } else if (generator_config.generator == "rhg") {
      std::cout << "Parameters for Random Hyperbolic Graphs RHG(n,gamma,d)" << std::endl;
      std::cout << "================================================" << std::endl;
      std::cout << "=== Random Hyperbolic Graphs RHG(n,gamma,d) ====" << std::endl;
      std::cout << "================================================" << std::endl;
      std::cout << "Parameters:" << std::endl;
      std::cout << "-n\t\t<number of vertices as a power of two>" << std::endl;
      std::cout << "-gamma\t\t<power-law exponent>" << std::endl;
      std::cout << "-d\t\t<average degree>" << std::endl;
      std::cout << "-k\t\t<number of chunks>" << std::endl;
      std::cout << "-seed\t\t<seed for PRNGs>" << std::endl;
      std::cout << "-output\t\t<output file>" << std::endl;
      std::cout << "\nExample:" << std::endl;
      std::cout << "mpirun -n 16 ./build/app/kagen -gen rhg -n 20 -d 8 -gamma 2.2 -output tmp" << std::endl;
    } else if (generator_config.generator == "rmat") {
      std::cout << "Parameters for Kronecker Graphs RMAT(n,m)" << std::endl;
      std::cout << "================================================" << std::endl;
      std::cout << "=========== Kronecker Graphs RMAT(n,m) ==========" << std::endl;
      std::cout << "================================================" << std::endl;
      std::cout << "Parameters:" << std::endl;
      std::cout << "-n\t\t<number of vertices as a power of two>" << std::endl;
      std::cout << "-m\t\t<number of edges as a power of two>" << std::endl;
      std::cout << "-k\t\t<number of chunks>" << std::endl;
      std::cout << "-seed\t\t<seed for PRNGs>" << std::endl;
      std::cout << "-output\t\t<output file>" << std::endl;
      std::cout << "\nExample:" << std::endl;
      std::cout << "mpirun -n 16 ./build/app/kagen -gen rmat -n 20 -m 22 -output tmp" << std::endl;
    }
    exit(0);
  }

  // Nodes
  bool exact_n = args.IsSet("exact_n");
  if (exact_n)
    generator_config.n = args.Get<ULONG>("n", 100);
  else
    generator_config.n = (ULONG)1 << args.Get<ULONG>("n", 3);

  // Blocks
  generator_config.k = args.Get<ULONG>("k", size);

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

}
#endif
