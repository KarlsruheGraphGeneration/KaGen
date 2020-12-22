/*******************************************************************************
 * include/definitions.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _DEFINITIONS_H_
#define _DEFINITIONS_H_

namespace kagen {

// Constants
typedef long long LONG;
typedef unsigned long long ULONG;
typedef int INT;
typedef unsigned int UINT;
typedef int PEID;

const PEID ROOT = 0;

// High/low prec
typedef long double HPFloat;
typedef double LPFloat;
typedef ULONG SInt;
typedef LONG SSInt;

enum Direction {
  Up, Down, Left, Right, Front, Back
};

}
#endif
