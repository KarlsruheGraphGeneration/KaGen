#include "kagen/generators/hyper/h_erdos/hyper_gnp.h"

#include "kagen/sampling/hash.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace kagen {

std::unique_ptr<Generator>
HyperGNPFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<HyperGNPBig>(config, rank, size);
}

template <typename BigInt>
HyperGNP<BigInt>::HyperGNP(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size),
      rng_(config) {
    if (config_.n == 0) {
        throw ConfigurationError("HGNP requires n > 0");
    }

    if (config_.k == 0) {
        throw ConfigurationError("HGNP requires at least one chunk");
    }
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateCSR() {
    graph_.hyperedge_offsets.push_back(0);

    const SInt lower_bound = config_.size_dist_lower_bound;
    const SInt upper_bound =
        config_.size_dist_upper_bound != 0 ? std::min(config_.size_dist_upper_bound, config_.n) : config_.n;

    for (SInt hyperedge_size = lower_bound; hyperedge_size <= upper_bound; ++hyperedge_size) {
        const double expected_edges = 1000.0 / std::pow(2.0, hyperedge_size - lower_bound);

        const CountInt universe = Binomial(config_.n, hyperedge_size);
        if (universe == 0) {
            continue;
        }

        const long double log_binom = std::lgammal(static_cast<long double>(config_.n) + 1.0L)
                                      - std::lgammal(static_cast<long double>(hyperedge_size) + 1.0L)
                                      - std::lgammal(static_cast<long double>(config_.n - hyperedge_size) + 1.0L);

        const long double p_ld = static_cast<long double>(expected_edges) / expl(log_binom);
        const double      p    = static_cast<double>(std::min<long double>(1.0L, p_ld));

        if (rank_ == 0) {
            std::cout << "k=" << hyperedge_size << " expected=" << expected_edges << " p=" << p << '\n';
        }

        if (p <= 0.0) {
            continue;
        }

        GenerateHyperedgesOfSizeGNP(hyperedge_size, universe, p);
    }

    const SInt vertices_per_pe = config_.n / size_;
    const SInt remainder       = config_.n % size_;

    const SInt start = rank_ * vertices_per_pe + std::min<SInt>(rank_, remainder);
    const SInt end   = start + vertices_per_pe + (rank_ < remainder);

    SetVertexRange(start, end);
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateHyperedgesOfSizeGNP(
    const SInt hyperedge_size, const CountInt& /*universe*/, const double p) {
    if (hyperedge_size <= 0 || hyperedge_size > config_.n || p <= 0.0) {
        return;
    }

    if (p > 1.0) {
        throw ConfigurationError("HGNP probabilities must be in [0, 1]");
    }

    const SInt global_first_begin = 0;
    const SInt global_first_end   = config_.n - hyperedge_size + 1;

    const SInt local_first_begin = FirstVertexBegin(hyperedge_size);
    const SInt local_first_end   = FirstVertexEnd(hyperedge_size);

    QueryFirstVertexRangeGNP(
        hyperedge_size, global_first_begin, global_first_end, local_first_begin, local_first_end, p, 1);
}

template <typename BigInt>
void HyperGNP<BigInt>::QueryFirstVertexRangeGNP(
    const SInt hyperedge_size, const SInt lo, const SInt hi, const SInt query_lo, const SInt query_hi, const double p,
    const SInt level) {
    if (lo >= hi || p <= 0.0 || query_lo >= hi || query_hi <= lo) {
        return;
    }

    const CountInt universe = CountFirstVertexRange(lo, hi, config_.n, hyperedge_size);

    if (universe == 0) {
        return;
    }

    if (query_lo <= lo && hi <= query_hi && universe <= static_cast<CountInt>(std::numeric_limits<BigInt>::max())) {
        const BigInt universe_rng = CheckedCastCount(universe);

        const SInt seed = sampling::Spooky::hash(config_.seed + 104729 * hyperedge_size + 9176 * level + lo);

        const SInt m = rng_.GenerateBinomial(seed, universe_rng, p);

        if (m == 0) {
            return;
        }

        std::vector<SInt> prefix;
        QueryPrefix(prefix, hyperedge_size, lo, hi, m, level);
        return;
    }

    if (hi - lo == 1) {
        const SInt        first = lo;
        std::vector<SInt> prefix;
        prefix.reserve(hyperedge_size);
        prefix.push_back(first);

        QueryPrefixGNP(prefix, hyperedge_size - 1, first + 1, config_.n, p, level + 1);
        return;
    }

    const SInt mid = lo + ((hi - lo) / 2);

    QueryFirstVertexRangeGNP(hyperedge_size, lo, mid, query_lo, query_hi, p, level + 1);
    QueryFirstVertexRangeGNP(hyperedge_size, mid, hi, query_lo, query_hi, p, level + 1);
}
template <typename BigInt>
void HyperGNP<BigInt>::QueryPrefixGNP(
    std::vector<SInt>& prefix, const SInt remaining_k, const SInt lo, const SInt hi, const double p, const SInt level) {
    if (remaining_k == 0) {
        const SInt seed =
            sampling::Spooky::hash(config_.seed + (99991 * level) + (17 * static_cast<SInt>(prefix.front())));

        const SInt m = rng_.GenerateBinomial(seed, static_cast<BigInt>(1), p);
        if (m == 1) {
            PushHyperedge(prefix);
        }
        return;
    }

    if (lo >= hi || hi - lo < remaining_k || p <= 0.0) {
        return;
    }

    const SInt max_next = hi - remaining_k;

    for (SInt next = lo; next <= max_next; ++next) {
        const CountInt universe = Binomial(hi - next - 1, remaining_k - 1);

        if (universe == 0) {
            continue;
        }

        prefix.push_back(next);

        if (universe <= static_cast<CountInt>(std::numeric_limits<BigInt>::max())) {
            const BigInt universe_rng = CheckedCastCount(universe);

            const SInt seed = sampling::Spooky::hash(
                config_.seed + 4256249 * remaining_k + 741457 * level + 17 * next
                + 131 * static_cast<SInt>(prefix.front()));

            const SInt m = rng_.GenerateBinomial(seed, universe_rng, p);

            if (m != 0) {
                rng_.GenerateSample(seed, universe_rng, m, [&](const BigInt sample) {
                    std::vector<SInt> pins = prefix;

                    std::vector<SInt> suffix = UnrankCombination(sample - 1, hi - next - 1, remaining_k - 1);

                    for (const SInt vertex: suffix) {
                        pins.push_back(next + 1 + vertex);
                    }

                    PushHyperedge(pins);
                });
            }
        } else {
            QueryPrefixGNP(prefix, remaining_k - 1, next + 1, hi, p, level + 1);
        }

        prefix.pop_back();
    }
}
template <typename BigInt>
void HyperGNP<BigInt>::QueryPrefix(
    std::vector<SInt>& prefix, const SInt remaining_k, const SInt lo, const SInt hi, const SInt m, const SInt level) {
    if (m == 0 || remaining_k == 0 || lo >= hi) {
        return;
    }

    const CountInt universe = CountFirstVertexRange(lo, hi, config_.n, remaining_k);

    if (universe == 0) {
        return;
    }

    if (remaining_k == 1) {
        const BigInt universe_rng = CheckedCastCount(static_cast<CountInt>(hi - lo));

        const SInt seed = sampling::Spooky::hash(config_.seed + 99173 * level + lo + 17 * prefix.size());

        rng_.GenerateSample(seed, universe_rng, m, [&](const BigInt sample) {
            std::vector<SInt> pins = prefix;
            pins.push_back(lo + static_cast<SInt>(sample - 1));
            PushHyperedge(pins);
        });

        return;
    }

    if (universe <= static_cast<CountInt>(std::numeric_limits<BigInt>::max())) {
        const BigInt universe_rng = CheckedCastCount(universe);

        const SInt seed = sampling::Spooky::hash(config_.seed + 130363 * remaining_k + 2750159 * level + lo);

        rng_.GenerateSample(seed, universe_rng, m, [&](const BigInt sample) {
            std::vector<SInt> pins = prefix;

            std::vector<SInt> suffix = UnrankFirstVertexRange(sample - 1, lo, hi, config_.n, remaining_k);

            pins.insert(pins.end(), suffix.begin(), suffix.end());
            PushHyperedge(pins);
        });

        return;
    }

    const SInt mid = lo + ((hi - lo) / 2);

    const CountInt left_weight  = CountFirstVertexRange(lo, mid, config_.n, remaining_k);
    const CountInt right_weight = universe - left_weight;

    if (left_weight == 0) {
        QueryPrefix(prefix, remaining_k, mid, hi, m, level + 1);
        return;
    }

    if (right_weight == 0) {
        QueryPrefix(prefix, remaining_k, lo, mid, m, level + 1);
        return;
    }

    const SInt seed = sampling::Spooky::hash(config_.seed + 271828 * remaining_k + 1618033 * level + lo);

    const SInt left_m = rng_.GenerateHypergeometricLarge(seed, left_weight, m, universe);

    QueryPrefix(prefix, remaining_k, lo, mid, left_m, level + 1);
    QueryPrefix(prefix, remaining_k, mid, hi, m - left_m, level + 1);
}

template <typename BigInt>
void HyperGNP<BigInt>::FinalizeCSR(MPI_Comm /*comm*/) {}

template <typename BigInt>
CountInt HyperGNP<BigInt>::Binomial(SInt n, SInt k) const {
    if (k < 0 || k > n) {
        return 0;
    }

    k = std::min(k, n - k);

    CountInt result = 1;

    for (SInt i = 1; i <= k; ++i) {
        result *= static_cast<CountInt>(n - k + i);
        result /= static_cast<CountInt>(i);
    }

    return result;
}

template <typename BigInt>
SInt HyperGNP<BigInt>::FirstVertexBegin(const SInt hyperedge_size) const {
    const SInt num_first_vertices = config_.n - hyperedge_size + 1;

    return (num_first_vertices * static_cast<SInt>(rank_)) / static_cast<SInt>(size_);
}

template <typename BigInt>
SInt HyperGNP<BigInt>::FirstVertexEnd(const SInt hyperedge_size) const {
    const SInt num_first_vertices = config_.n - hyperedge_size + 1;

    return (num_first_vertices * static_cast<SInt>(rank_ + 1)) / static_cast<SInt>(size_);
}

template <typename BigInt>
BigInt HyperGNP<BigInt>::CheckedCastCount(const CountInt& value) const {
    if (value > static_cast<CountInt>(std::numeric_limits<BigInt>::max())) {
        throw ConfigurationError("hyperedge universe exceeds RNG integer range");
    }

    return static_cast<BigInt>(value);
}

template <typename BigInt>
BigInt HyperGNP<BigInt>::BinomialNative(SInt n, SInt k) const {
    if (k < 0 || k > n) {
        return 0;
    }

    k = std::min(k, n - k);

    BigInt result = 1;

    for (SInt i = 1; i <= k; ++i) {
        result *= static_cast<BigInt>(n - k + i);
        result /= static_cast<BigInt>(i);
    }

    return result;
}

template <typename BigInt>
std::vector<SInt> HyperGNP<BigInt>::UnrankCombination(BigInt index, const SInt n, const SInt k) const {
    std::vector<SInt> pins;
    pins.reserve(k);

    SInt x = 0;

    for (SInt remaining = k; remaining > 0; --remaining) {
        BigInt count = BinomialNative(n - x - 1, remaining - 1);

        for (; x < n; ++x) {
            if (index < count) {
                pins.push_back(x);
                ++x;
                break;
            }

            index -= count;

            const SInt a = n - x - 1;
            const SInt b = remaining - 1;

            if (a > b && a > 0) {
                count *= static_cast<BigInt>(a - b);
                count /= static_cast<BigInt>(a);
            } else {
                count = 0;
            }
        }
    }

    return pins;
}

template <typename BigInt>
void HyperGNP<BigInt>::PushHyperedge(const std::vector<SInt>& pins) {
    graph_.hyperedge_pins.insert(graph_.hyperedge_pins.end(), pins.begin(), pins.end());
    graph_.hyperedge_offsets.push_back(static_cast<SInt>(graph_.hyperedge_pins.size()));
}

template <typename BigInt>
CountInt
HyperGNP<BigInt>::CountFirstVertexRange(const SInt lo, const SInt hi, const SInt n, const SInt hyperedge_size) const {
    if (lo >= hi || hyperedge_size == 0 || lo >= n) {
        return 0;
    }

    const SInt clipped_hi = std::min<SInt>(hi, n - hyperedge_size + 1);

    const CountInt left  = Binomial(n - lo, hyperedge_size);
    const CountInt right = Binomial(n - clipped_hi, hyperedge_size);

    return left > right ? left - right : CountInt{0};
}

template <typename BigInt>
std::vector<SInt>
HyperGNP<BigInt>::UnrankFirstVertexRange(BigInt index, const SInt lo, const SInt hi, const SInt n, const SInt k) const {
    SInt left  = lo;
    SInt right = hi;

    while (left + 1 < right) {
        const SInt mid = left + ((right - left) / 2);

        const BigInt left_count = CheckedCastCount(CountFirstVertexRange(lo, mid, n, k));

        if (index < left_count) {
            right = mid;
        } else {
            left = mid;
        }
    }

    const SInt first = left;

    const BigInt before = CheckedCastCount(CountFirstVertexRange(lo, first, n, k));

    index -= before;

    std::vector<SInt> result;
    result.reserve(k);

    result.push_back(first);

    std::vector<SInt> suffix = UnrankCombination(index, n - first - 1, k - 1);

    for (const SInt vertex: suffix) {
        result.push_back(first + 1 + vertex);
    }

    return result;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
template class HyperGNP<__int128>;
#pragma GCC diagnostic pop

} // namespace kagen