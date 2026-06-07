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

    bool FitsNative(const CountInt& value) const {
        return value <= static_cast<CountInt>(std::numeric_limits<int_t>::max());
    }

    int_t CheckedCast(const CountInt& value) const {
        if (!FitsNative(value)) {
            throw ConfigurationError("count exceeds RNG native integer range");
        }

        return static_cast<int_t>(value);
    }
    /*
    SInt GenerateHypergeometricLarge(
        const SInt seed, const CountInt& left_population, const SInt draws, const CountInt& total_population) {
        if (draws == 0 || left_population == 0 || total_population == 0) {
            return 0;
        }

        if (left_population >= total_population) {
            return draws;
        }

        // Exact path: same behavior as existing GNM when counts fit int_t.
        if (FitsNative(left_population) && FitsNative(total_population)) {
            return static_cast<SInt>(GenerateHypergeometric(
                seed, CheckedCast(left_population), static_cast<int_t>(draws), CheckedCast(total_population)));
        }

        // Large-count path:
        // Hypergeometric(N, K, draws) ≈ Binomial(draws, K / N)
        // This is accurate when total_population is huge compared to draws.
        const CountFloat p_float = CountFloat(left_population) / CountFloat(total_population);

        double p = p_float.convert_to<double>();
        p        = std::clamp(p, 0.0, 1.0);

        return GenerateBinomial(seed, draws, static_cast<LPFloat>(p));
    }

    ######################################
        TODO: Temp exact Hypergeometric
        Distribution for large populations
    ######################################*/
    SInt GenerateHypergeometricLarge(
        const SInt seed, const CountInt& left_population, const SInt draws, const CountInt& total_population) {
        if (draws == 0 || left_population == 0 || total_population == 0) {
            return 0;
        }

        if (left_population >= total_population) {
            return draws;
        }

        if (FitsNative(left_population) && FitsNative(total_population)) {
            return static_cast<SInt>(GenerateHypergeometric(
                seed, CheckedCast(left_population), static_cast<int_t>(draws), CheckedCast(total_population)));
        }

        rng_.seed(seed);

        CountInt remaining_left  = left_population;
        CountInt remaining_total = total_population;

        SInt successes = 0;

        for (SInt draw = 0; draw < draws; ++draw) {
            if (remaining_left == 0) {
                break;
            }

            if (remaining_left == remaining_total) {
                successes += draws - draw;
                break;
            }

            const CountInt sample = GenerateUniformBelow(remaining_total);

            if (sample < remaining_left) {
                ++successes;
                --remaining_left;
            }

            --remaining_total;
        }

        return successes;
    }

    std::size_t BitLength(CountInt value) const {
        std::size_t bits = 0;

        while (value > 0) {
            value >>= 1;
            ++bits;
        }

        return bits;
    }

    CountInt GenerateUniformBelow(const CountInt& upper_bound) {
        if (upper_bound <= 0) {
            return 0;
        }

        const std::size_t bits = BitLength(upper_bound);

        while (true) {
            CountInt value = 0;

            for (std::size_t generated = 0; generated < bits; generated += 32) {
                value <<= 32;
                value += static_cast<std::uint32_t>(rng_());
            }

            const std::size_t excess_bits = ((bits + 31) / 32) * 32 - bits;

            if (excess_bits > 0) {
                value >>= excess_bits;
            }

            if (value < upper_bound) {
                return value;
            }
        }
    }

private:
    const PGeneratorConfig& config_;

    std::mt19937                                         rng_;
    sampling::hypergeometric_distribution<int_t, double> hyp_;
};

} // namespace kagen
