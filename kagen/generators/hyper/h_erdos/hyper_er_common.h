#pragma once

#include "kagen/generators/generator.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"
#include "kagen/tools/mersenne.h"
#include "kagen/tools/rng_wrapper.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/unordered/unordered_flat_map.hpp>

#include <algorithm>
#include <fstream>
#include <limits>
#include <unordered_set>
#include <vector>

namespace kagen {

using CountInt = boost::multiprecision::cpp_int;

constexpr unsigned long long kCountSeedMultiplier    = 999999ULL;
constexpr unsigned long long kRankSeedMultiplier     = 104729ULL;
constexpr unsigned long long kEdgeSeedMultiplier     = 1234567ULL;
constexpr unsigned long long kEdgeRankSeedMultiplier = 932456ULL;
constexpr unsigned long long kGNMCountSeedMultiplier = 7919ULL;
constexpr std::size_t        kFNVOffsetBasis         = 1469598103934665603ULL;
constexpr std::size_t        kFNVPrime               = 1099511628211ULL;

struct MinOwnerWeights {
    std::vector<long double> weights;

    // actual_count(s) = weights[s] * exp(log_actual_scale)
    long double log_actual_scale = 0.0L;
};

struct MinOwnerPartition {
    std::vector<SInt>        begin;
    std::vector<SInt>        end;
    std::vector<long double> mass;
    long double              total_mass = 0.0L;
};

struct VectorHash {
    std::size_t operator()(const std::vector<SInt>& vector) const {
        std::size_t hash_value = kFNVOffsetBasis;

        for (const SInt value: vector) {
            hash_value ^= std::hash<SInt>{}(value);
            hash_value *= kFNVPrime;
        }

        return hash_value;
    }
};

long double LogBinomialApprox(SInt n, SInt k);

struct LogBinomCacheStats {
    std::uint64_t calls    = 0;
    std::uint64_t hits     = 0;
    std::uint64_t misses   = 0;
    std::uint64_t inserts  = 0;
    std::uint64_t max_size = 0;

    SInt min_key = std::numeric_limits<SInt>::max();
    SInt max_key = 0;

    std::uint64_t backward_steps = 0;
    std::uint64_t forward_steps  = 0;
    std::uint64_t same_as_last   = 0;
};

struct LogBinomCache {
    SInt                                  fixed_k = -1;
    std::unordered_map<SInt, long double> values;
    double                                log_k_factorial = 0.0;

    explicit LogBinomCache(SInt k = -1, std::size_t expected_size = 4096) : fixed_k(k) {
        values.max_load_factor(0.7);
        values.reserve(expected_size);

        if (fixed_k >= 0 && fixed_k <= 64) {
            for (SInt i = 2; i <= fixed_k; ++i) {
                log_k_factorial += std::log(static_cast<double>(i));
            }
        }
    }

    long double LogBinomialSmallKFast(SInt n, SInt k) const {
        if (k > n) {
            return -std::numeric_limits<long double>::infinity();
        }

        k = std::min(k, n - k);

        double result = 0.0;
        for (SInt i = 0; i < k; ++i) {
            result += std::log(static_cast<double>(n - i));
        }

        result -= log_k_factorial;
        return static_cast<long double>(result);
    }

    long double Get(SInt x, SInt k) {
        auto it = values.find(x);
        if (it != values.end()) {
            return it->second;
        }

        const long double value = (k <= 32 && k == fixed_k) ? LogBinomialSmallKFast(x, k) : LogBinomialApprox(x, k);

        values.emplace(x, value);
        return value;
    }
};

bool ExpectedCountIsNegligible(double expected_m_k);

std::vector<SInt> DeterministicRankCounts(SInt m, const std::vector<long double>& mass, long double total_mass);

long double LogDifferenceOfExponentials(long double log_larger, long double log_smaller);

long double LogAdd(long double log_a, long double log_b);

long double MinRangeMassApprox(SInt begin, SInt end, SInt n, SInt k);

long double MinRangeMassApproxCached(SInt begin, SInt end, SInt n, SInt k, LogBinomCache& cache);

SInt FindMinBoundaryByMass(SInt n, SInt k, SInt rank, SInt size);

CountInt BinomialExact(SInt n, SInt k);

CountInt MinRangeMassExact(SInt begin, SInt end, SInt n, SInt k);

bool MinRangeMassDefinitelyExceedsSInt(SInt begin, SInt end, SInt n, SInt k);

long double MersenneUniform01(Mersenne& mersenne);

long double BinomSmallExactLD(SInt n, SInt q);

std::vector<SInt> FloydSample(SInt universe_offset, SInt universe_size, SInt sample_size, Mersenne& mersenne);

template <typename RNG>
std::vector<SInt> MultinomialRankCounts(
    const SInt m, const std::vector<long double>& mass, const long double total_mass, RNG& rng, SInt seed) {
    const SInt p = static_cast<SInt>(mass.size());

    std::vector<long double> cdf(p);
    long double              prefix = 0.0L;

    for (SInt rank = 0; rank < p; ++rank) {
        prefix += mass[rank] / total_mass;
        cdf[rank] = prefix;
    }

    cdf[p - 1] = 1;

    std::vector<SInt> counts(p, 0);

    for (SInt i = 0; i < m; ++i) {
        const SInt        draw_seed = sampling::Spooky::hash(seed++);
        const long double x         = static_cast<long double>(rng.GenerateUniform(draw_seed, 0.0L, 1.0L));

        const auto it   = std::lower_bound(cdf.begin(), cdf.end(), x);
        const SInt rank = static_cast<SInt>(it - cdf.begin());

        ++counts[rank];
    }

    return counts;
}

template <typename RNG>
SInt SampleMinimumFromCDF(const SInt min_begin, const std::vector<long double>& cdf, RNG& rng, SInt& seed) {
    const SInt        draw_seed = sampling::Spooky::hash(seed++);
    const long double x         = static_cast<long double>(rng.GenerateUniform(draw_seed, 0.0L, 1.0L));

    const auto it = std::lower_bound(cdf.begin(), cdf.end(), x);

    return min_begin + static_cast<SInt>(it - cdf.begin());
}

template <typename RNG>
SInt PoissonLocalCountFromScaledMass(
    const long double local_mass, const long double log_actual_scale, const double p, RNG& rng, const SInt seed) {
    if (local_mass <= 0.0L || p <= 0.0) {
        return 0;
    }

    const long double log_lambda = std::log(static_cast<long double>(p)) + std::log(local_mass) + log_actual_scale;

    if (log_lambda > std::log(static_cast<long double>(std::numeric_limits<SInt>::max()))) {
        throw ConfigurationError("HGNP expected local hyperedge count exceeds SInt range");
    }

    const double lambda = static_cast<double>(expl(log_lambda));

    if (!std::isfinite(lambda) || lambda <= 0.0) {
        return 0;
    }

    return rng.GeneratePoisson(seed, lambda);
}

template <typename RNG>
SInt SampleMinimumImplicitK2(const SInt local_begin, const SInt local_end, const SInt n, RNG& rng, SInt& seed) {
    const long double u = std::min<long double>(
        static_cast<long double>(rng.GenerateCanonicalDoubleStream()), std::nextafter(1.0L, 0.0L));
    ++seed;

    const auto choose2 = [](const long double x) {
        return (x * (x - 1.0L)) / 2.0L;
    };

    const long double begin_tail = choose2(static_cast<long double>(n - local_begin));
    const long double end_tail   = (n - local_end >= 2) ? choose2(static_cast<long double>(n - local_end)) : 0.0L;

    const long double target_tail = begin_tail - u * (begin_tail - end_tail);

    long double x_real = (1.0L + std::sqrt(1.0L + 8.0L * target_tail)) / 2.0L;

    SInt x = static_cast<SInt>(std::floor(x_real));

    SInt s = n - x - 1;
    s      = std::clamp<SInt>(s, local_begin, local_end - 1);

    // Correct possible rounding error.
    while (s > local_begin) {
        const long double prev_tail = choose2(static_cast<long double>(n - s));

        if (prev_tail > target_tail) {
            break;
        }

        --s;
    }

    while (s < local_end - 1) {
        const long double tail = choose2(static_cast<long double>(n - (s + 1)));

        if (tail <= target_tail) {
            break;
        }

        ++s;
    }

    return s;
}

template <typename RNG>
SInt SampleMinimumImplicit(
    const SInt local_begin, const SInt local_end, const SInt n, const SInt k, RNG& rng, SInt& seed,
    LogBinomCache& cache, std::uint64_t* steps = nullptr, std::uint64_t* cache_gets = nullptr) {
    if (k == 2) {
        return SampleMinimumImplicitK2(local_begin, local_end, n, rng, seed);
    }

    const long double uniform01 = std::min<long double>(
        static_cast<long double>(rng.GenerateCanonicalDoubleStream()), std::nextafter(1.0L, 0.0L));
    ++seed;

    if (cache_gets) {
        ++(*cache_gets);
    }
    const long double log_begin_tail = cache.Get(n - local_begin, k);

    if (cache_gets) {
        ++(*cache_gets);
    }
    const long double log_end_tail =
        (n - local_end >= k) ? cache.Get(n - local_end, k) : -std::numeric_limits<long double>::infinity();

    const long double log_local_mass = LogDifferenceOfExponentials(log_begin_tail, log_end_tail);

    const long double log_target_mass = log_local_mass + std::log(uniform01);

    const long double log_target_tail = LogDifferenceOfExponentials(log_begin_tail, log_target_mass);

    SInt left  = local_begin;
    SInt right = local_end - 1;

    while (left < right) {
        if (steps) {
            ++(*steps);
        }

        const SInt mid = left + ((right - left) / 2);

        if (cache_gets) {
            ++(*cache_gets);
        }
        const long double log_mid_tail =
            (n - (mid + 1) >= k) ? cache.Get(n - (mid + 1), k) : -std::numeric_limits<long double>::infinity();

        if (log_mid_tail <= log_target_tail) {
            right = mid;
        } else {
            left = mid + 1;
        }
    }

    return left;
}

template <typename RNG>
SInt SampleBlockCountFromLogSize(
    long double log_block_size, long double log_p, RNG& rng, SInt seed, const char* error_context) {
    const long double log_expected = log_block_size + log_p;

    if (log_expected <= std::log(std::numeric_limits<double>::denorm_min())) {
        return 0;
    }

    if (log_expected > std::log(static_cast<long double>(std::numeric_limits<SInt>::max()))) {
        throw ConfigurationError(std::string(error_context) + " expected local block count exceeds SInt");
    }

    if (log_block_size < std::log(static_cast<long double>(std::numeric_limits<SInt>::max()))) {
        const SInt trials = static_cast<SInt>(std::llround(expl(log_block_size)));
        return rng.GenerateBinomial(seed, trials, static_cast<double>(expl(log_p)));
    }

    return rng.GeneratePoisson(seed, static_cast<double>(expl(log_expected)));
}

template <typename RNG>
SInt HypergeometricCountSequential(CountInt population, CountInt successes, const SInt draws, RNG& rng, SInt& seed) {
    SInt count = 0;

    for (SInt i = 0; i < draws; ++i) {
        if (population == 0 || successes == 0) {
            break;
        }
        const SInt        draw_seed = sampling::Spooky::hash(seed++);
        const long double uniform01 = static_cast<long double>(rng.GenerateUniform(draw_seed, 0.0L, 1.0L));

        const long double p = successes.convert_to<long double>() / population.convert_to<long double>();

        if (uniform01 < p) {
            ++count;
            --successes;
        }

        --population;
    }
    return count;
}

std::vector<double> ReadNumericLines(const std::string& filename, const std::string& error_message);
// Fast approximate Boltzmann-style sampler for HGNM pin-budget size counts.
class BoltzmannPinBudgetSizeCountSampler {
public:
    BoltzmannPinBudgetSizeCountSampler(SInt m, SInt lower, SInt upper, SInt pin_budget)
        : m_(m),
          lower_(lower),
          d_(upper - lower),
          r_(pin_budget - (lower * m)) {}

    std::unordered_map<SInt, SInt> Sample() {
        if (r_ == 0) {
            return {{lower_, m_}};
        }

        if (r_ == d_ * m_) {
            return {{lower_ + d_, m_}};
        }
        const auto probs = BuildProbabilities();

        if (r_ == d_ * m_) {
            return {{lower_ + d_, m_}};
        }

        const SInt        active_d = static_cast<SInt>(probs.size()) - 1;
        std::vector<SInt> counts(probs.size(), 0);

        SInt edges_sum = 0;
        SInt pins_sum  = 0;

        std::vector<std::pair<long double, SInt>> fractional;
        fractional.reserve(probs.size());

        for (SInt t = 0; t <= active_d; ++t) {
            const long double exact = static_cast<long double>(m_) * probs[t];
            const SInt        base  = static_cast<SInt>(std::floor(exact));

            counts[t] = base;
            edges_sum += base;
            pins_sum += base * t;

            fractional.emplace_back(exact - static_cast<long double>(base), t);
        }

        std::sort(fractional.begin(), fractional.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

        for (SInt i = 0; edges_sum < m_; ++i) {
            const SInt t = fractional[i % fractional.size()].second;
            ++counts[t];
            ++edges_sum;
            pins_sum += t;
        }

        Repair(counts, edges_sum, pins_sum, active_d);

        std::unordered_map<SInt, SInt> size_counts;

        for (SInt t = 0; t <= active_d; ++t) {
            if (counts[t] > 0) {
                size_counts[lower_ + t] = counts[t];
            }
        }

        return size_counts;
    }

private:
    std::vector<long double> BuildProbabilities() const {
        const long double target = static_cast<long double>(r_) / static_cast<long double>(m_);

        if (target <= 0.0L) {
            return {1.0L};
        }

        if (target >= static_cast<long double>(d_)) {
            return {1.0L};
        }

        const long double q     = target / (target + 1.0L);
        const long double alpha = 1.0L - q;

        std::vector<long double> probs;
        probs.reserve(1024);

        long double z = 0.0L;
        long double p = alpha;

        for (SInt t = 0; t <= d_; ++t) {
            const long double expected = static_cast<long double>(m_) * p;

            if (t > 0 && expected < 1e-9L) {
                break;
            }

            probs.push_back(p);
            z += p;
            p *= q;
        }

        for (auto& x: probs) {
            x /= z;
        }

        return probs;
    }

    void Repair(std::vector<SInt>& counts, SInt& edges_sum, SInt& pins_sum, SInt active_d) {
        if (edges_sum != m_) {
            throw ConfigurationError("Boltzmann pin-budget sampler expected exact edge count before repair");
        }

        if (pins_sum < r_) {
            MoveUpBatched(counts, pins_sum, active_d);
        } else if (pins_sum > r_) {
            MoveDownBatched(counts, pins_sum, active_d);
        }
    }
    void MoveUpBatched(std::vector<SInt>& counts, SInt& pins_sum, SInt active_d) const {
        SInt delta = r_ - pins_sum;

        for (SInt from = 0; from < active_d && delta > 0; ++from) {
            if (counts[from] == 0) {
                continue;
            }

            const SInt to   = std::min<SInt>(active_d, from + delta);
            const SInt gain = to - from;

            if (gain <= 0) {
                continue;
            }

            const SInt move = std::min<SInt>(counts[from], delta / gain);

            if (move > 0) {
                counts[from] -= move;
                counts[to] += move;
                pins_sum += move * gain;
                delta -= move * gain;
            }

            if (delta > 0 && counts[from] > 0 && from + 1 <= active_d) {
                --counts[from];
                ++counts[from + 1];
                ++pins_sum;
                --delta;
            }
        }

        if (pins_sum != r_) {
            throw ConfigurationError("Failed to repair HGNM pin budget upward");
        }
    }

    void MoveDownBatched(std::vector<SInt>& counts, SInt& pins_sum, SInt active_d) const {
        SInt delta = pins_sum - r_;

        for (SInt from = active_d; from > 0 && delta > 0; --from) {
            if (counts[from] == 0) {
                continue;
            }

            const SInt to   = std::max<SInt>(0, from - delta);
            const SInt loss = from - to;

            if (loss <= 0) {
                continue;
            }

            const SInt move = std::min<SInt>(counts[from], delta / loss);

            if (move > 0) {
                counts[from] -= move;
                counts[to] += move;
                pins_sum -= move * loss;
                delta -= move * loss;
            }

            if (delta > 0 && counts[from] > 0 && from - 1 >= 0) {
                --counts[from];
                ++counts[from - 1];
                --pins_sum;
                --delta;
            }
        }

        if (pins_sum != r_) {
            throw ConfigurationError("Failed to repair HGNM pin budget downward");
        }
    }

    SInt m_;
    SInt lower_;
    SInt d_;
    SInt r_;
};
} // namespace kagen