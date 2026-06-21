/*******************************************************************************
 * include/tools/rng_wrapper.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/kagen.h"
#include "kagen/sampling/methodR.hpp"

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <limits>
#include <random>

namespace kagen {
template <typename int_t = std::int64_t>
class RNGWrapper {
    using CountInt   = boost::multiprecision::cpp_int;
    using CountFloat = boost::multiprecision::cpp_dec_float_50;

public:
    RNGWrapper(const PGeneratorConfig& config) : config_(config), rng_(0), hyp_(0) {}

    int_t GenerateHypergeometric(SInt seed, int_t n, int_t m, int_t N) {
        SInt variate = 0;
        if (config_.use_binom)
            variate = GenerateBinomial(seed, n, (LPFloat)m / N);
        else {
            hyp_.seed(seed);
            if (m < 1)
                return 0;
            variate = hyp_(n, N - n, m);
        }
        return variate;
    }

    SInt GenerateBinomial(SInt seed, SInt n, LPFloat p) {
        rng_.seed(seed);
        std::binomial_distribution<SInt> bin(n, p);
        return bin(rng_);
    }

    template <typename F>
    void GenerateSample(SInt seed, SInt N, SInt n, F&& callback) {
        sampling::HashSampling<>      hs(seed, config_.base_size);
        sampling::SeqDivideSampling<> sds(hs, config_.base_size, seed, config_.use_binom);
        sds.sample(N, n, callback);
    }

    template <typename real_t = double>
    real_t GenerateUniform(SInt seed, real_t min = 0.0, real_t max = 1.0) {
        rng_.seed(seed);

        std::uniform_real_distribution<real_t> dist(min, max);

        return dist(rng_);
    }
    void SeedUniformStream(const SInt seed) {
        rng_.seed(seed);
    }

    double GenerateUniformDoubleStream() {
        return uniform_double_(rng_);
    }

    long double GenerateUniformLongDoubleStream() {
        return uniform_long_double_(rng_);
    }

    double GenerateCanonicalDoubleStream() {
        constexpr double scale = 1.0 / 9007199254740992.0; // 2^-53

        const std::uint64_t a = static_cast<std::uint64_t>(rng_()) >> 5; // 27 bits
        const std::uint64_t b = static_cast<std::uint64_t>(rng_()) >> 6; // 26 bits

        return static_cast<double>((a << 26) | b) * scale;
    }

    SInt GenerateIntegerStream(const SInt upper_exclusive) {
        if (upper_exclusive <= 0) {
            throw ConfigurationError("Invalid integer stream range");
        }

        std::uniform_int_distribution<SInt> dist(0, upper_exclusive - 1);
        return dist(rng_);
    }

    SInt GeneratePoisson(SInt seed, double lambda) {
        rng_.seed(seed);

        std::poisson_distribution<SInt> poisson(lambda);

        return poisson(rng_);
    }

    SInt GenerateBinomialHuge(SInt seed, CountInt n, LPFloat p) {
        if (p <= 0.0) {
            return 0;
        }

        if (p >= 1.0) {
            if (n > std::numeric_limits<SInt>::max()) {
                throw ConfigurationError("Binomial variate exceeds SInt range");
            }

            return n.convert_to<SInt>();
        }

        if (n <= std::numeric_limits<SInt>::max()) {
            return GenerateBinomial(seed, n.convert_to<SInt>(), p);
        }

        const CountFloat mean = CountFloat(n) * CountFloat(p);

        if (mean > CountFloat(std::numeric_limits<SInt>::max())) {
            throw ConfigurationError("Expected binomial variate exceeds SInt range");
        }

        // Temporary exact-path stopgap:
        // mathematically exact huge-n binomial still missing.
        throw ConfigurationError("Huge-population exact binomial sampler not implemented yet");
    }

private:
    const PGeneratorConfig& config_;

    std::mt19937                                         rng_;
    std::uniform_real_distribution<double>               uniform_double_{0.0, 1.0};
    std::uniform_real_distribution<long double>          uniform_long_double_{0.0L, 1.0L};
    sampling::hypergeometric_distribution<int_t, double> hyp_;
};

} // namespace kagen
