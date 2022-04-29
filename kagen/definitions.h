/*******************************************************************************
 * include/definitions.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <tuple>
#include <vector>

#include <mpi.h>

namespace kagen {
// Constants
typedef long long          LONG;
typedef unsigned long long ULONG;
typedef int                INT;
typedef unsigned int       UINT;
typedef int                PEID;

const PEID ROOT = 0;

// High/low prec
typedef long double HPFloat;
typedef double      LPFloat;
typedef ULONG       SInt;
typedef LONG        SSInt;

#define KAGEN_MPI_LONG    MPI_LONG_LONG
#define KAGEN_MPI_ULONG   MPI_UNSIGNED_LONG_LONG
#define KAGEN_MPI_INT     MPI_INT
#define KAGEN_MPI_UINT    MPI_UNSIGNED
#define KAGEN_MPI_PEID    MPI_INT
#define KAGEN_MPI_HPFLOAT MPI_LONG_DOUBLE
#define KAGEN_MPI_LPFLOAT MPI_DOUBLE
#define KAGEN_MPI_SINT    MPI_UNSIGNED_LONG_LONG
#define KAGEN_MPI_SSINT   MPI_LONG_LONG

enum Direction { Up, Down, Left, Right, Front, Back };

using EdgeList    = std::vector<std::tuple<SInt, SInt>>;
using VertexRange = std::pair<SInt, SInt>;
} // namespace kagen
