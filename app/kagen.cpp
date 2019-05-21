/*******************************************************************************
 * app/generate_kagen.cpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Daniel Funke <funke@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
//#define DEL_STATS 1

#include <mpi.h>

#include "benchmark.h"
#include "generator_config.h"
#include "io/generator_io.h"
#include "parse_parameters.h"
#include "timer.h"

#include "geometric/delaunay/delaunay_2d.h"
#include "geometric/delaunay/delaunay_3d.h"
#include "geometric/rgg/rgg_2d.h"
#include "geometric/rgg/rgg_3d.h"
#include "gnm/gnm_directed.h"
#include "gnm/gnm_undirected.h"
#include "gnp/gnp_directed.h"
#include "gnp/gnp_undirected.h"
#include "hyperbolic/hyperbolic.h"
#include "barabassi/barabassi.h"
#include "kronecker/kronecker.h"

using namespace kagen;

void OutputParameters(const PGeneratorConfig &config, const PEID /* rank */,
                      const PEID size) {
  if (config.generator == "gnm_directed" ||
      config.generator == "gnm_undirected" ||
      config.generator == "gnp_directed" ||
      config.generator == "gnp_undirected")
    std::cout << "generate graph (n=" << config.n << ", m=" << config.m
              << " (p=" << config.p << "), k=" << config.k
              << ", s=" << config.seed << ", P=" << size << ")" << std::endl;

  else if (config.generator == "rgg_2d" || config.generator == "rgg_3d")
    std::cout << "generate graph (n=" << config.n << ", r=" << config.r
              << ", k=" << config.k << ", s=" << config.seed << ", P=" << size
              << ")" << std::endl;

  else if (config.generator == "rdg_2d" || config.generator == "rdg_3d")
    std::cout << "generate graph (n=" << config.n << ", k=" << config.k
              << ", s=" << config.seed << ", P=" << size << ")" << std::endl;

  else if (config.generator == "rhg")
    std::cout << "generate graph (n=" << config.n << ", d=" << config.avg_degree
              << ", gamma=" << config.plexp << ", k=" << config.k
              << ", s=" << config.seed << ", P=" << size << ")" << std::endl;

  else if (config.generator == "ba")
    std::cout << "generate graph (n=" << config.n << ", d=" << config.min_degree
              << ", k=" << config.k << ", s=" << config.seed << ", P=" << size
              << ")" << std::endl;

  else if (config.generator == "rmat")
    std::cout << "generate graph (n=" << config.n << ", m=" << config.m
              << ", k=" << config.k << ", s=" << config.seed << ", P=" << size
              << ")" << std::endl;
}

template <typename Generator, typename EdgeCallback>
void RunGenerator(const PGeneratorConfig &config, const PEID rank,
                  const PEID /* size */, Statistics &stats, Statistics &edge_stats,
                  Statistics &edges, const EdgeCallback &cb) {
  // Start timers
  Timer t;
  double local_time = 0.0;
  double total_time = 0.0;
  t.Restart();

  // Chunk distribution
  Generator gen(config, rank, cb);
  gen.Generate();

  // Output
  local_time = t.Elapsed();
  MPI_Reduce(&local_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, ROOT,
             MPI_COMM_WORLD);
  if (rank == ROOT) {
    stats.Push(total_time);
    edge_stats.Push(total_time / gen.NumberOfEdges());
    edges.Push(gen.NumberOfEdges());
  }

  if (rank == ROOT) std::cout << "write output..." << std::endl;
  gen.Output();
}

int main(int argn, char **argv) {
  // Init MPI
  MPI_Init(&argn, &argv);
  PEID rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Read command-line args
  PGeneratorConfig generator_config;
  ParseParameters(argn, argv, rank, size, generator_config);

  if (rank == ROOT) OutputParameters(generator_config, rank, size);

  // Statistics
  Statistics stats;
  Statistics edge_stats;
  Statistics edges;

  auto edge_cb = [](SInt, SInt){};
  ULONG user_seed = generator_config.seed;
  for (ULONG i = 0; i < generator_config.iterations; ++i) {
    MPI_Barrier(MPI_COMM_WORLD);
    generator_config.seed = user_seed + i;
    if (generator_config.generator == "gnm_directed")
      RunGenerator<GNMDirected<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "gnm_undirected")
      RunGenerator<GNMUndirected<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "gnp_directed")
      RunGenerator<GNPDirected<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "gnp_undirected")
      RunGenerator<GNPUndirected<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "rgg_2d")
      RunGenerator<RGG2D<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "rgg_3d")
      RunGenerator<RGG3D<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "rdg_2d")
      RunGenerator<Delaunay2D<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "rdg_3d")
      RunGenerator<Delaunay3D<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "rhg")
      RunGenerator<Hyperbolic<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "ba")
      RunGenerator<Barabassi<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "rmat")
      RunGenerator<Kronecker<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else 
      if (rank == ROOT) std::cout << "generator not supported" << std::endl;
  }

  if (rank == ROOT) {
    std::cout << "RESULT runner=" << generator_config.generator
              << " time=" << stats.Avg() << " stddev=" << stats.Stddev()
              << " iterations=" << generator_config.iterations
              << " edges=" << edges.Avg()
              << " time_per_edge=" << edge_stats.Avg() << std::endl;
  }

  MPI_Finalize();
  return 0;
}
