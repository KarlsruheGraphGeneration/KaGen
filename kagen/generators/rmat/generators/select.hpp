/*******************************************************************************
 * rmat/rng/stl.hpp
 *
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef RMAT_GENERATORS_SELECT_HEADER
#define RMAT_GENERATORS_SELECT_HEADER

#include "dSFMT.hpp"
#include "stl.hpp"
#ifdef RMAT_HAVE_MKL
#include "mkl.hpp"
#endif

namespace rmat {
namespace generators {

struct select {
#ifdef RMAT_HAVE_MKL
    // MKL is much faster than anything else, by a factor that's not even funny
    // any more for large block sizes
    using type = mkl;
#else
    // dSFMT is at least twice as fast as std::mt19937_64 for large blocks with
    // gcc, and more using clang. It's practically never slower, so prefer it.
    using type = dSFMT;
#endif // RMAT_HAVE_MKL
};

using select_t = select::type;

} // namespace generators
} // namespace rmat

#endif // RMAT_GENERATORS_SELECT_HEADER
