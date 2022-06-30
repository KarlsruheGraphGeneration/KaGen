/*******************************************************************************
 * rmat/generators/mkl.hpp
 *
 * Generate random deviates using Intel MKL
 *
 * Copyright (C) 2018 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef RMAT_GENERATORS_MKL_HEADER
#define RMAT_GENERATORS_MKL_HEADER

#ifdef RMAT_HAVE_MKL

#include <tlx/define.hpp>
#include <tlx/logger.hpp>

#include <mkl.h>
#include <mkl_vsl.h>

#include <cassert>
#include <cmath>
#include <limits>

namespace rmat {
namespace generators {

void CheckVslError(int);

/*!
 * MKL generator wrapper
 *
 * If you need more elements at a time than an int can hold, you need to
 * `#define MKL_INT size_t` before including this file
 */
class mkl {
public:
    static const char* name;

    //! Initialize new MKL random generator
    mkl(size_t seed) {
        vslNewStream(&stream, VSL_BRNG_SFMT19937, seed);
    }

    ~mkl() {
        vslDeleteStream(&stream);
    }

    //! non-copyable: delete copy-constructor
    mkl(const mkl &) = delete;
    //! non-copyable: delete assignment operator
    mkl & operator = (const mkl &) = delete;
    //! move-constructor: default
    mkl(mkl &&) = default;
    //! move-assignment operator: default
    mkl & operator = (mkl &&) = default;

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
    void generate_block(std::vector<double> &output, size_t size) {
        check_size(size);
        if (size > output.size()) {
            output.resize(size);
        }
        generate_block(output.data(), size);
    }

    //! Generate `size` uniform doubles from [0, 1)
    void generate_block(double *ptr, size_t size) {
        check_size(size);

        MKL_INT count = static_cast<MKL_INT>(size);
        // VSL_RNG_METHOD_UNIFORM_STD_ACCURATE?
        int status = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream,
                                  count, ptr, 0.0, 1.0);
        CheckVslError(status);
    }

    //! Generate `size` log(uniform double) values
    void generate_log_block(double *output, size_t size) {
        generate_block(output, size);

        MKL_INT count = static_cast<MKL_INT>(size);
        // TODO check whether MKL supports inplace log
        vdLn(count, output, output);
    }

    //! Generate `size` log(uniform double) values
    void generate_log_block(std::vector<double> &output, size_t size) {
        generate_block(output, size); // resizes if needed

        MKL_INT count = static_cast<MKL_INT>(size);
        // TODO check whether MKL supports inplace log
        vdLn(count, output.data(), output.data());
    }

    //! Generate `size` uniform integers from [min, max] (i.e., both inclusive).
    //! Unfortunately, MKL does not support generating 64-bit random integers.
    void generate_int_block(int min, int max, int *output, size_t size) {
        check_size(size);
        MKL_INT count = static_cast<MKL_INT>(size);
        // VSL_RNG_METHOD_UNIFORM_STD_ACCURATE?
        int status = viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream,
                                  count, output, min, max + 1);
        CheckVslError(status);
    }

    //! Generate `size` uniform integers from [min, max] (i.e., both inclusive).
    //! Unfortunately, MKL does not support generating 64-bit random integers.
    void generate_int_block(int min, int max, std::vector<int> &output,
                            size_t size) {
        check_size(size);
        if (size > output.size()) {
            output.resize(size);
        }
        generate_int_block(min, max, output.data(), size);
    }

    //! Generate `size` geometrically distributed integers with parameter p
    void generate_geometric_block(double p, int *output, size_t size) {
        check_size(size);
        MKL_INT count = static_cast<MKL_INT>(size);
        int status = viRngGeometric(VSL_RNG_METHOD_GEOMETRIC_ICDF, stream,
                                    count, output, p);
        CheckVslError(status);
    }

    //! Generate `size` geometrically distributed integers with parameter p
    void generate_geometric_block(double p, std::vector<int> &output,
                                  size_t size) {
        check_size(size);
        if (size > output.size()) {
            output.resize(size);
        }
        generate_geometric_block(p, output.data(), size);
    }

    //! Generate `size` exponentially distributed integers with rate `lambda`
    //! and displacement `displacement`
    void generate_exponential_block(double lambda, double *output,
                                    size_t size, double displacement = 0) {
        check_size(size);

        double scale = 1.0 / lambda;

        MKL_INT count = static_cast<MKL_INT>(size);
        int status = vdRngExponential(VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream,
                                      count, output, displacement, scale);
        CheckVslError(status);
    }

    //! Generate `size` exponentially distributed integers with rate `lambda`
    //! and displacement `displacement`
    void generate_exponential_block(double lambda, std::vector<double> &output,
                                    size_t size, double displacement = 0) {
        check_size(size);
        if (size > output.size()) {
            output.resize(size);
        }
        generate_exponential_block(lambda, output.data(), size, displacement);
    }

    //! Generate `size` normally distributed integers with mean `mean` and
    //! standard deviation `stdev`
    void generate_gaussian_block(double mean, double stdev,
                                 double *output, size_t size) {
        check_size(size);
        MKL_INT count = static_cast<MKL_INT>(size);
        int status = vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream,
                                   count, output, mean, stdev);
        CheckVslError(status);
    }

    //! Generate `size` normally distributed integers with mean `mean` and
    //! standard deviation `stdev`
    void generate_gaussian_block(double mean, double stdev,
                                 std::vector<double> &output, size_t size) {
        check_size(size);
        if (size > output.size()) {
            output.resize(size);
        }
        generate_gaussian_block(mean, stdev, output.data(), size);
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
            block_size_ = std::max(block_size_,
                                   minimum_reasonable_block_size());
            // generate_block takes care of resizing the vector for us
            generate_block(randblock_, block_size_);
            index_ = 0;
            block_id_++;
        }
        return randblock_[index_++];
    }

    TLX_ATTRIBUTE_ALWAYS_INLINE
    double next_exponential(double lambda) {
        if (TLX_UNLIKELY(index_ >= block_size_)) {
            if (block_id_ > 2 && ((block_id_ + 1) & block_id_) == 0) {
                // block_id_ + 1 is a power of two. We appear to need a lot of
                // random numbers, increase the blocksize to reduce RNG overhead
                block_size_ *= 2;
            }
            block_size_ = std::max(block_size_,
                                   minimum_reasonable_block_size());
            // generate_log_block takes care of resizing the vector for us
            generate_log_block(randblock_, block_size_);
            index_ = 0;
            block_id_++;
        }
        return -randblock_[index_++] / lambda;

    }

    //! Generate a uniform integer from [min, max] (i.e., both inclusive)
    template <typename int_t>
    TLX_ATTRIBUTE_ALWAYS_INLINE
    int_t next_int(int_t min, int_t max) {
        return next() * (max - min + 1) + min;
    }

    //! Bernoulli trial with success probability p
    TLX_ATTRIBUTE_ALWAYS_INLINE
    bool next_bernoulli(double p) {
        assert(0 <= p && p <= 1);
        return next() < p;
    }

    //! Bernoulli trial with success probability cutoff/max
    TLX_ATTRIBUTE_ALWAYS_INLINE
    bool next_bernoulli(double cutoff, double max) {
        assert(0 <= cutoff && cutoff <= max);
        return next() * max < cutoff;
    }

    //! Generate a normally distributed value. This needs two uniform deviates,
    //! if you need more than one, look at next_two_gaussians.
    TLX_ATTRIBUTE_ALWAYS_INLINE
    double next_gaussian(double mean, double stdev) {
        double U = next(), V = next();
        double a = stdev * std::sqrt(-2*std::log(U));
        double b = 2 * M_PI * V;

        return mean + a * std::cos(b);
    }

    //! Generate two independent normally distributed values
    TLX_ATTRIBUTE_ALWAYS_INLINE
    std::pair<double, double> next_two_gaussians(double mean, double stdev) {
        double U = next(), V = next();
        double a = stdev * std::sqrt(-2*std::log(U));
        double b = 2 * M_PI * V;

        U = mean + a * std::cos(b);
        V = mean + a * std::sin(b);
        return std::make_pair(U, V);
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
            LOG1
                << "Error: MKL generator block size exceeds value range of MKL_INT:"
                << size << " >= " << std::numeric_limits<MKL_INT>::max();
            return false;
        }
        return true;
    }

    //! MKL state object
    VSLStreamStatePtr stream;
    std::vector<double> randblock_;
    size_t index_ = 0, block_size_ = 0, block_id_ = 0;
};

} // namespace generators
} // namespace rmat

#endif // RMAT_HAVE_MKL

#endif // RMAT_GENERATORS_MKL_HEADER
