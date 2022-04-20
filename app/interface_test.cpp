/*******************************************************************************
 * app/interface_test.cpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Daniel Funke <funke@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "kagen_interface.h"

using namespace kagen;

int main(int argn, char **argv) {
  // Init MPI
  MPI_Init(&argn, &argv);
  PEID rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Params
  SInt num_vertices = (ULONG)1 << 16;
  SInt num_edges = (ULONG)1 << 20;

  // Run generator
  KaGen gen(rank, size);
  gen.GenerateDirectedGNM(num_vertices, num_edges);
  gen.GenerateUndirectedGNM(num_vertices, num_edges);
  gen.Generate2DRGG(num_vertices, 0.0072);
  gen.Generate3DRGG(num_vertices, 0.01);
  gen.Generate2DRDG(num_vertices);
  gen.Generate3DRDG(num_vertices);
  gen.GenerateRHG(num_vertices, 3.0, 16);
  gen.GenerateBA(num_vertices, 16);
  const SInt sqrt_num_vertices = std::sqrt(num_vertices);
  gen.Generate2DGrid(sqrt_num_vertices, sqrt_num_vertices, 0.9, true);
  
  MPI_Finalize();
  return 0;
}
