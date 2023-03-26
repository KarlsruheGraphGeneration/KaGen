/*******************************************************************************
 * sampling/rng/mkl.hpp
 *
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef SAMPLING_RNG_MKL_HEADER
    #define SAMPLING_RNG_MKL_HEADER

    #ifdef SAMPLING_HAVE_MKL

        #include "kagen/sampling/rng/errcheck.inc"

        #include <mkl.h>
        #include <mkl_vsl.h>

        #include <limits>

namespace sampling {
namespace rng {

/*!
 * MKL generator wrapper
 *
 * If you need more elements at a time than an int can hold, you need to
 * `#define MKL_INT size_t` before including this file
 */
class mkl {
public:
    static constexpr const char* name = "mkl";

    //! Initialize new MKL random generator
    mkl(size_t seed) {
        vslNewStream(&stream, VSL_BRNG_SFMT19937, seed);
    }

    ~mkl() {
        vslDeleteStream(&stream);
    }

    //! Re-seed the random generator
    void seed(size_t seed) {
        vslDeleteStream(&stream);
        vslNewStream(&stream, VSL_BRNG_SFMT19937, seed);
    }

    //! Minimum number of elements that needs to be generated at a time
    size_t minimum_block_size() const {
        return 1; // it's unknown how many MKL generates internally
    }

    //! Minimum number of elements that needs to be generated at a time for
    //! reasonable performance
    size_t minimum_reasonable_block_size() const {
        return 512; // chosen by fair guess?
    }

    //! Generate `size` uniform doubles from [0, 1)
    void generate_block(std::vector<double>& output, size_t size) {
        check_size(size);
        if (size > output.size()) {
            output.resize(size);
        }

        MKL_INT count = static_cast<MKL_INT>(size);
        // VSL_RNG_METHOD_UNIFORM_STD_ACCURATE?
        int status = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, count, output.data(), 0.0, 1.0);
        CheckVslError(status);
    }

    //! Generate `size` uniform integers from [min, max] (i.e., both inclusive).
    //! Unfortunately, MKL does not support generating 64-bit random integers.
    void generate_int_block(int min, int max, std::vector<int>& output, size_t size) {
        check_size(size);
        if (size > output.size()) {
            output.resize(size);
        }

        MKL_INT count = static_cast<MKL_INT>(size);
        // VSL_RNG_METHOD_UNIFORM_STD_ACCURATE?
        int status = viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, count, output.data(), min, max + 1);
        CheckVslError(status);
    }

    //! Generate `size` geometrically distributed integers with parameter p
    void generate_geometric_block(double p, std::vector<int>& output, size_t size) {
        check_size(size);

        MKL_INT count  = static_cast<MKL_INT>(size);
        int     status = viRngGeometric(VSL_RNG_METHOD_GEOMETRIC_ICDF, stream, count, output.data(), p);
        CheckVslError(status);
    }

    //! Get a single [0,1) double. Computes increasingly large blocks internally
    //! so that this is fast.  May block while next block is generated.
    TLX_ATTRIBUTE_ALWAYS_INLINE
    double next() {
        if (TLX_UNLIKELY(index_ >= block_size_)) {
            if (block_id_ > 2 && ((block_id_ + 1) & block_id_) == 0) {
                // block_id_ + 1 is a power of two. We appear to need a lot of
                // random numbers, increase the blocksize to reduce RNG overhead
                block_size_ *= 2;
            }
            block_size_ = std::max(block_size_, minimum_reasonable_block_size());
            // generate_block takes care of resizing the vector for us
            generate_block(randblock_, block_size_);
            index_ = 0;
            block_id_++;
        }
        return randblock_[index_++];
    }

    //! Generate a uniform integer from [min, max] (i.e., both inclusive)
    template <typename int_t>
    TLX_ATTRIBUTE_ALWAYS_INLINE int_t next_int(int_t min, int_t max) {
        return next() * (max - min + 1) + min;
    }

    //! Alias for next()
    TLX_ATTRIBUTE_ALWAYS_INLINE
    double operator()() {
        return next();
    }

private:
    //! Check that `size` fits into an MKL_INT
    bool check_size(size_t size) {
        if (size >= std::numeric_limits<MKL_INT>::max()) {
            /*
            LOG1
                << "Error: MKL generator block size exceeds value range of MKL_INT:"
                << size << " >= " << std::numeric_limits<MKL_INT>::max();
            */
            return false;
        }
        return true;
    }

    //! MKL state object
    VSLStreamStatePtr   stream;
    std::vector<double> randblock_;
    size_t              index_ = 0, block_size_ = 0, block_id_ = 0;
};

} // namespace rng
} // namespace sampling

    #endif // SAMPLING_HAVE_MKL

#endif // SAMPLING_RNG_MKL_HEADER
