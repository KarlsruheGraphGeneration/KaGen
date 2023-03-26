/*******************************************************************************
 * sampling/rng/stl.hpp
 *
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef SAMPLING_RNG_SELECT_HEADER
    #define SAMPLING_RNG_SELECT_HEADER

    #include "kagen/sampling/config.hpp"
    #include "kagen/sampling/rng/dSFMT.hpp"
    #include "kagen/sampling/rng/stl.hpp"
    #ifdef SAMPLING_HAVE_MKL
        #include "kagen/sampling/rng/mkl.hpp"
    #endif

namespace sampling {
namespace rng {

struct select {
    #ifdef SAMPLING_HAVE_MKL
    // MKL is much faster than anything else, by a factor that's not even funny
    // any more for large block sizes
    using type = mkl;
    #else
    // dSFMT is at least twice as fast as std::mt19937_64 for large blocks with
    // gcc, and more using clang. It's practically never slower, so prefer it.
    using type = dSFMT;
    #endif
};

using select_t = select::type;

} // namespace rng
} // namespace sampling

#endif // SAMPLING_RNG_SELECT_HEADER
