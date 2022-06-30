/*******************************************************************************
 * rmat/generators/stl.hpp
 *
 * Generate random deviates using the STL
 *
 * Copyright (C) 2018 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef RMAT_GENERATORS_STL_HEADER
#define RMAT_GENERATORS_STL_HEADER

#include <cassert>
#include <random>
#include <type_traits>
#include <vector>

namespace rmat {
namespace generators {

//! Wrapper around STL random generators and distributions
template <typename Generator = std::mt19937_64>
class stl {
public:
    static constexpr const char* name = "stl";

    using generator_t = Generator;

    //! Initialize new stl random generator
    stl(size_t seed) : generator_(seed == 0 ? std::random_device{}() : seed) {}

    stl(const stl & other) {
        generator_ = other.generator_;
    }
    stl & operator = (const stl &other) {
        generator_ = other.generator_;
        return *this;
    }
    //! move-constructor: default
    stl(stl &&) = default;
    //! move-assignment operator: default
    stl & operator = (stl &&) = default;

    //! Re-seed the random generator
    void seed(size_t seed) {
        if (seed == 0) seed = std::random_device{}();
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

    double next_exponential(double lambda) {
        std::exponential_distribution<double> distribution(lambda);
        return distribution(generator_);
    }

    //! Generate a uniform integer from [min, max] (i.e., both inclusive)
    template <typename int_t>
    int_t next_int(int_t min, int_t max) {
        std::uniform_int_distribution<int_t> dist(min, max);
        return dist(generator_);
    }

    //! Bernoulli trial with success probability p
    bool next_bernoulli(double p) {
        assert(0 <= p && p <= 1);
        std::bernoulli_distribution dist(p);
        return dist(generator_);
    }

    //! Bernoulli trial with success probability cutoff/max
    bool next_bernoulli(double cutoff, double max) {
        assert(0 <= cutoff && cutoff <= max);
        return next() * max < cutoff;
    }

    //! Binomial distribution
    size_t next_binomial(size_t n, double p) {
        std::binomial_distribution<> dist(n, p);
        return dist(generator_);
    }

    //! Generate a normally distributed value. This needs two uniform deviates,
    //! if you need more than one, look at next_two_gaussians.
    double next_gaussian(double mean, double stdev) {
        std::normal_distribution<double> dist(mean, stdev);
        return dist(generator_);
    }

    //! Generate two independent normally distributed values
    std::pair<double, double> next_two_gaussians(double mean, double stdev) {
        std::normal_distribution<double> dist(mean, stdev);
        double U = dist(generator_), V = dist(generator_);
        return std::make_pair(U, V);
    }


    template <typename distribution_type, typename Iterator>
    void generate_block(Iterator begin, Iterator end,
                        distribution_type distribution) {
        using value_type = typename std::iterator_traits<Iterator>::value_type;

        for (auto it = begin; it < end; ++it) {
            *it = static_cast<value_type>(distribution(generator_));
        }
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

    //! Generate `size` uniform doubles from [0, 1)
    void generate_block(double *arr, size_t size) {
        for (size_t i = 0; i < size; ++i) {
            arr[i] = next();
        }
    }

    //! Generate `size` uniform integers from [min, max] (i.e., both inclusive)
    template <typename int_t>
    void generate_int_block(int_t min, int_t max, int_t *arr, size_t size)
    {
        std::uniform_int_distribution<int_t> dist(min, max);
        for (size_t i = 0; i < size; ++i) {
            arr[i] = dist(generator_);
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
        generate_int_block(min, max, output.data(), size);
    }

    //! Generate `size` geometrically integers with parameter p
    template <typename int_t>
    void generate_geometric_block(double p, int_t *arr, size_t size)
    {
        std::geometric_distribution<int_t> dist(p);
        for (size_t i = 0; i < size; ++i) {
            arr[i] = dist(generator_);
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
        generate_geometric_block(p, output.data(), size);
    }

    //! Generate `size` exponentially distributed integers with rate `lambda`
    //! and displacement `displacement`
    void generate_exponential_block(double lambda, double *arr, size_t size) {
        std::exponential_distribution<double> dist(lambda);
        for (size_t i = 0; i < size; ++i) {
            arr[i] = dist(generator_);
        }
    }

    //! Generate `size` exponentially distributed integers with rate `lambda`
    //! and displacement `displacement`
    void generate_exponential_block(double lambda, std::vector<double> &output,
                                    size_t size) {
        if (size > output.size()) {
            output.resize(size);
        }
        generate_exponential_block(lambda, output.data(), size);
    }

    //! Generate `size` normally distributed integers with mean `mean` and
    //! standard deviation `stdev`
    void generate_gaussian_block(double mean, double stdev,
                                 double *arr, size_t size) {
        std::normal_distribution<double> dist(mean, stdev);
        for (size_t i = 0; i < size; ++i) {
            arr[i] = dist(generator_);
        }
    }

    //! Generate `size` normally distributed integers with mean `mean` and
    //! standard deviation `stdev`
    void generate_gaussian_block(double mean, double stdev,
                                 std::vector<double> &output, size_t size) {
        if (size > output.size()) {
            output.resize(size);
        }
        generate_gaussian_block(mean, stdev, output.data(), size);
    }

    //! Alias for next()
    double operator()() {
        return next();
    }

private:
    generator_t generator_;
    std::uniform_real_distribution<double> uniform_double_;
};

} // namespace generators
} // namespace rmat

#endif // RMAT_GENERATORS_STL_HEADER
