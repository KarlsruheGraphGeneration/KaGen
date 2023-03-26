/*******************************************************************************
 * sampling/rng/stl.hpp
 *
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef SAMPLING_RAND_STL_HEADER
#define SAMPLING_RAND_STL_HEADER

#include <random>

namespace sampling {
namespace rng {

//! Wrapper around STL random generators and distributions
template <typename Generator = std::mt19937_64>
class stl {
public:
    static constexpr const char* name = "stl";

    using generator_t = Generator;

    //! Initialize new stl random generator
    stl(size_t seed) : generator_(seed) {}

    //! Re-seed the random generator
    void seed(size_t seed) {
        generator_.seed(seed);
    }

    //! Minimum number of elements that needs to be generated at a time
    size_t minimum_block_size() const {
        return 1;
    }

    //! Minimum number of elements that needs to be generated at a time for
    //! reasonable performance
    size_t minimum_reasonable_block_size() const {
        return 1;
    }

    //! Generate a uniform double from [0, 1)
    double next() {
        return uniform_double_(generator_);
    }

    //! Generate a uniform integer from [min, max] (i.e., both inclusive)
    template <typename int_t>
    int_t next_int(int_t min, int_t max) {
        std::uniform_int_distribution<int_t> dist(min, max);
        return dist(generator_);
    }

    //! Generate `size` uniform doubles from [0, 1)
    void generate_block(std::vector<double> &output, size_t size) {
        if (size > output.size()) {
            output.resize(size);
        }
        for (size_t i = 0; i < size; ++i) {
            output[i] = next();
        }
    }

    //! Generate `size` uniform integers from [min, max] (i.e., both inclusive)
    template <typename int_t>
    void generate_int_block(int_t min, int_t max, std::vector<int_t> &output,
                            size_t size)
    {
        if (size > output.size()) {
            output.resize(size);
        }
        std::uniform_int_distribution<int_t> dist(min, max);
        for (size_t i = 0; i < size; ++i) {
            output[i] = dist(generator_);
        }
    }

    //! Generate `size` geometrically integers with parameter p
    template <typename int_t>
    void generate_geometric_block(double p, std::vector<int_t> &output,
                                  size_t size)
    {
        if (size > output.size()) {
            output.resize(size);
        }
        std::geometric_distribution<int_t> dist(p);
        for (size_t i = 0; i < size; ++i) {
            output[i] = dist(generator_);
        }
    }

    //! Alias for next()
    double operator()() {
        return next();
    }

private:
    generator_t generator_;
    std::uniform_real_distribution<double> uniform_double_;
};

} // namespace rng
} // namespace sampling

#endif // SAMPLING_RAND_STL_HEADER
