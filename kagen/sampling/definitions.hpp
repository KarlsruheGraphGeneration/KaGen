/*******************************************************************************
 * sampling/config.hpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef SAMPLING_DEFINITIONS_HEADER
#define SAMPLING_DEFINITIONS_HEADER

namespace sampling {

// Constants
typedef long long LONG;
typedef unsigned long long ULONG;
typedef int INT;
typedef unsigned int UINT;
typedef unsigned long long NodeID;
typedef unsigned long long EdgeID;
typedef unsigned long long NodeWeight;
typedef unsigned long long EdgeWeight;
typedef int PEID;

const PEID ROOT = 0;

#define SAMPLING_SEEDA 28475421
// #define SAMPLING_SEEDB 52150599

} // namespace sampling

#endif // SAMPLING_DEFINITIONS_HEADER
