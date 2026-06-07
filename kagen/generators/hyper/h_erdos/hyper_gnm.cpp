#include "kagen/generators/hyper/h_erdos/hyper_gnm.h"

#include "kagen/sampling/hash.hpp"

#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace kagen {

std::unique_ptr<Generator>
HyperGNMFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<HyperGNMBig>(config, rank, size);
}

template <typename BigInt>
HyperGNM<BigInt>::HyperGNM(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size),
      rng_(config) {
    if (config_.n == 0) {
        throw ConfigurationError("HGNM requires n > 0");
    }

    if (config_.m == 0) {
        throw ConfigurationError("HGNM requires m > 0");
    }

    if (config_.k == 0) {
        throw ConfigurationError("HGNM requires at least one chunk");
    }
}

template <typename BigInt>
void HyperGNM<BigInt>::GenerateCSR() {
    graph_.hyperedge_offsets.push_back(0);

    const SInt lower_bound = config_.size_dist_lower_bound;
    const SInt upper_bound =
        config_.size_dist_upper_bound != 0 ? std::min(config_.size_dist_upper_bound, config_.n) : config_.n;

    for (SInt hyperedge_size = lower_bound; hyperedge_size <= upper_bound; ++hyperedge_size) {
        const SInt m_k = HyperedgesOfSize(hyperedge_size);

        if (m_k == 0) {
            if (hyperedge_size > lower_bound) {
                break;
            }
            continue;
        }

        GenerateHyperedgesOfSize(hyperedge_size, m_k);
    }

    const SInt vertices_per_pe = config_.n / size_;
    const SInt remainder       = config_.n % size_;

    const SInt start = (rank_ * vertices_per_pe) + std::min<SInt>(rank_, remainder);

    const SInt end = start + vertices_per_pe + (rank_ < remainder);

    SetVertexRange(start, end);
}

template <typename BigInt>
void HyperGNM<BigInt>::FinalizeCSR(MPI_Comm /*comm*/) {}

template <typename BigInt>
SInt HyperGNM<BigInt>::HyperedgesOfSize(const SInt hyperedge_size) const {
    const SInt lower_bound = config_.size_dist_lower_bound;

    const SInt upper_bound =
        config_.size_dist_upper_bound != 0 ? std::min(config_.size_dist_upper_bound, config_.n) : config_.n;

    const double alpha = config_.size_dist_alpha;

    if (hyperedge_size < lower_bound || hyperedge_size > upper_bound) {
        return 0;
    }

    if (hyperedge_size == lower_bound) {
        SInt assigned = 0;

        for (SInt k = lower_bound + 1; k <= upper_bound; ++k) {
            assigned += HyperedgesOfSize(k);
        }

        return config_.m - assigned;
    }

    const double probability = alpha * std::pow(1.0 - alpha, static_cast<double>(hyperedge_size - lower_bound));

    return static_cast<SInt>(std::llround(static_cast<double>(config_.m) * probability));
}

template <typename BigInt>
CountInt HyperGNM<BigInt>::Binomial(SInt n, SInt k) const {
    if (k > n) {
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
CountInt
HyperGNM<BigInt>::CountFirstVertexRange(const SInt lo, const SInt hi, const SInt n, const SInt hyperedge_size) const {
    if (lo >= hi || hyperedge_size == 0 || lo >= n) {
        return 0;
    }

    const SInt clipped_hi = std::min<SInt>(hi, n - hyperedge_size + 1);

    // Number of k-hyperedges whose smallest vertex lies in [lo, hi):
    //
    //   sum_{v=lo}^{hi-1} C(n - v - 1, k - 1)
    // = C(n - lo, k) - C(n - hi, k)
    //
    const CountInt left  = Binomial(n - lo, hyperedge_size);
    const CountInt right = Binomial(n - clipped_hi, hyperedge_size);

    return left > right ? left - right : CountInt{0};
}

template <typename BigInt>
void HyperGNM<BigInt>::GenerateHyperedgesOfSize(const SInt hyperedge_size, const SInt m_k) {
    if (hyperedge_size > config_.n || m_k == 0) {
        return;
    }

    const SInt global_first_begin = 0;
    const SInt global_first_end   = config_.n - hyperedge_size + 1;

    const SInt local_first_begin = FirstVertexBegin(hyperedge_size);
    const SInt local_first_end   = FirstVertexEnd(hyperedge_size);

    QueryFirstVertexRange(
        hyperedge_size, m_k, global_first_begin, global_first_end, local_first_begin, local_first_end, 1);
}

template <typename BigInt>
BigInt HyperGNM<BigInt>::CheckedCastCount(const CountInt& value) const {
    if (value > static_cast<CountInt>(std::numeric_limits<BigInt>::max())) {
        throw ConfigurationError("local hyperedge universe exceeds RNG integer range");
    }

    return static_cast<BigInt>(value);
}

template <typename BigInt>
SInt HyperGNM<BigInt>::FirstVertexBegin(const SInt hyperedge_size) const {
    const SInt num_first_vertices = config_.n - hyperedge_size + 1;
    return (num_first_vertices * static_cast<SInt>(rank_)) / static_cast<SInt>(size_);
}

template <typename BigInt>
SInt HyperGNM<BigInt>::FirstVertexEnd(const SInt hyperedge_size) const {
    const SInt num_first_vertices = config_.n - hyperedge_size + 1;
    return (num_first_vertices * static_cast<SInt>(rank_ + 1)) / static_cast<SInt>(size_);
}

template <typename BigInt>
void HyperGNM<BigInt>::QueryFirstVertexRange(
    const SInt hyperedge_size, const SInt m, const SInt lo, const SInt hi, const SInt query_lo, const SInt query_hi,
    const SInt level) {
    if (m == 0 || lo >= hi || query_lo >= hi || query_hi <= lo) {
        return;
    }

    const CountInt total_weight = CountFirstVertexRange(lo, hi, config_.n, hyperedge_size);

    if (total_weight == 0) {
        return;
    }

    // Entire current range belongs to this rank
    if (query_lo <= lo && hi <= query_hi) {
        std::vector<SInt> prefix;
        QueryPrefix(prefix, hyperedge_size, lo, hi, m, level);
        return;
    }

    const SInt     mid          = lo + ((hi - lo) / 2);
    const CountInt left_weight  = CountFirstVertexRange(lo, mid, config_.n, hyperedge_size);
    const CountInt right_weight = total_weight - left_weight;

    if (left_weight == 0) {
        QueryFirstVertexRange(hyperedge_size, m, mid, hi, query_lo, query_hi, level + 1);
        return;
    }

    if (right_weight == 0) {
        QueryFirstVertexRange(hyperedge_size, m, lo, mid, query_lo, query_hi, level + 1);
        return;
    }

    const SInt seed    = sampling::Spooky::hash(config_.seed + (104729 * hyperedge_size) + (9176 * level) + lo);
    const SInt left_m  = rng_.GenerateHypergeometricLarge(seed, left_weight, m, total_weight);
    const SInt right_m = m - left_m;

    QueryFirstVertexRange(hyperedge_size, left_m, lo, mid, query_lo, query_hi, level + 1);
    QueryFirstVertexRange(hyperedge_size, right_m, mid, hi, query_lo, query_hi, level + 1);
}

template <typename BigInt>
void HyperGNM<BigInt>::QueryPrefix(
    std::vector<SInt>& prefix, const SInt remaining_k, const SInt lo, const SInt hi, const SInt m, const SInt level) {
    if (m == 0 || remaining_k == 0 || lo >= hi) {
        return;
    }

    const CountInt universe = CountFirstVertexRange(lo, hi, config_.n, remaining_k);

    if (universe == 0) {
        return;
    }

    // Suffix of size 1: sample from universe size
    if (remaining_k == 1) {
        const BigInt universe_rng = CheckedCastCount(static_cast<CountInt>(hi - lo));
        const SInt   seed         = sampling::Spooky::hash(config_.seed + (99173 * level) + lo + (17 * prefix.size()));

        rng_.GenerateSample(seed, universe_rng, m, [&](const BigInt sample) {
            std::vector<SInt> pins = prefix;
            pins.push_back(lo + static_cast<SInt>(sample - 1));
            PushHyperedge(pins);
        });
        return;
    }

    // remaining universe fits in range of GenerateSample
    if (universe <= static_cast<CountInt>(std::numeric_limits<BigInt>::max())) {
        const BigInt universe_rng = CheckedCastCount(universe);

        const SInt seed = sampling::Spooky::hash(config_.seed + (130363 * remaining_k) + (2750159 * level) + lo);
        rng_.GenerateSample(seed, universe_rng, m, [&](const BigInt sample) {
            std::vector<SInt> pins   = prefix;
            const BigInt      index  = sample - 1;
            std::vector<SInt> suffix = UnrankFirstVertexRange(index, lo, hi, config_.n, remaining_k);

            for (const SInt vertex: suffix) {
                pins.push_back(vertex);
            }
            PushHyperedge(pins);
        });

        return;
    }

    const SInt     mid         = lo + ((hi - lo) / 2);
    const CountInt left_weight = CountFirstVertexRange(lo, mid, config_.n, remaining_k);

    const CountInt right_weight = universe - left_weight;

    if (left_weight == 0) {
        QueryPrefix(prefix, remaining_k, mid, hi, m, level + 1);
        return;
    }

    if (right_weight == 0) {
        QueryPrefix(prefix, remaining_k, lo, mid, m, level + 1);
        return;
    }
    const SInt seed   = sampling::Spooky::hash(config_.seed + (271828 * remaining_k) + (1618033 * level) + lo);
    const SInt left_m = rng_.GenerateHypergeometricLarge(seed, left_weight, m, universe);

    QueryPrefix(prefix, remaining_k, lo, mid, left_m, level + 1);
    QueryPrefix(prefix, remaining_k, mid, hi, m - left_m, level + 1);
}

template <typename BigInt>
std::vector<SInt>
HyperGNM<BigInt>::UnrankFirstVertexRange(BigInt index, const SInt lo, const SInt hi, const SInt n, const SInt k) const {
    SInt left  = lo;
    SInt right = hi;

    while (left + 1 < right) {
        const SInt mid = left + (right - left) / 2;

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

    auto suffix = UnrankCombination(index, n - first - 1, k - 1);

    for (const SInt vertex: suffix) {
        result.push_back(first + 1 + vertex);
    }

    return result;
}

template <typename BigInt>
BigInt HyperGNM<BigInt>::BinomialNative(SInt n, SInt k) const {
    if (k > n) {
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
std::vector<SInt> HyperGNM<BigInt>::UnrankCombination(BigInt index, const SInt n, const SInt k) const {
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

            // Rolling update of binomial
            const SInt a = n - x - 1;
            const SInt b = remaining - 1;
            if (a > b && a > 0) {
                count *= (a - b);
                count /= a;
            } else {
                count = 0;
            }
        }
    }

    return pins;
}

template <typename BigInt>
void HyperGNM<BigInt>::PushHyperedge(const std::vector<SInt>& pins) {
    graph_.hyperedge_pins.insert(graph_.hyperedge_pins.end(), pins.begin(), pins.end());
    graph_.hyperedge_offsets.push_back(static_cast<SInt>(graph_.hyperedge_pins.size()));
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
template class HyperGNM<__int128>;
#pragma GCC diagnostic pop

} // namespace kagen