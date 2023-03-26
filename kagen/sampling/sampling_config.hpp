/*******************************************************************************
 * sampling/sampling_config.hpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef SAMPLING_SAMPLING_CONFIG_HEADER
#define SAMPLING_SAMPLING_CONFIG_HEADER

#include "kagen/sampling/definitions.hpp"

#include <string>

namespace sampling {

// Configuration for the generator.
struct SamplingConfig {
    SamplingConfig() {}

    // Seed for the PRNG
    int seed;
    // Number of samples
    ULONG n;
    // Size of population
    ULONG N;
    // Base case size
    ULONG k;
    // Sample probability (Bernoulli)
    double p;
    // Output filename
    std::string output_file;
    // Write edges directly to disk
    bool write_to_disk;
    // Number of iterations
    ULONG iterations;
    // Use binomial approximation
    bool use_binom;

    // Shared Mem OMP
    void LogDump(FILE* /* out */) const {}
};

} // namespace sampling

#endif // SAMPLING_SAMPLING_CONFIG_HEADER
