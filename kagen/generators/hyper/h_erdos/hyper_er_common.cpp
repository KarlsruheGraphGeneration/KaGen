#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace kagen {

MinOwnerWeights ComputeMinOwnerWeights(const SInt n, const SInt k) {
    if (k == 0 || k > n) {
        throw std::runtime_error("invalid k");
    }

    const SInt num_minima = n - k + 1;

    std::vector<long double> log_w(num_minima);
    log_w[0] = 0.0L;

    for (SInt s = 0; s + 1 < num_minima; ++s) {
        const long double numerator   = static_cast<long double>(n - s - k);
        const long double denominator = static_cast<long double>(n - s - 1);

        log_w[s + 1] = log_w[s] + std::log(numerator) - std::log(denominator);
    }

    const long double max_log = *std::max_element(log_w.begin(), log_w.end());

    MinOwnerWeights result;
    result.weights.resize(log_w.size());

    for (std::size_t i = 0; i < log_w.size(); ++i) {
        result.weights[i] = expl(log_w[i] - max_log);
    }

    result.log_actual_scale = LogBinomialApprox(n - 1, k - 1) + max_log;

    return result;
}

MinOwnerPartition PartitionMinimaByWeight(const std::vector<long double>& weights, const SInt p) {
    const SInt num_minima = static_cast<SInt>(weights.size());

    long double total = 0.0L;
    for (const long double weight: weights) {
        total += weight;
    }

    MinOwnerPartition part;
    part.begin.resize(p);
    part.end.resize(p);
    part.mass.assign(p, 0.0L);
    part.total_mass = total;

    SInt        s      = 0;
    long double prefix = 0.0L;

    for (SInt rank = 0; rank < p; ++rank) {
        part.begin[rank] = s;

        const long double target = total * static_cast<long double>(rank + 1) / static_cast<long double>(p);

        while (s < num_minima && prefix + weights[s] <= target) {
            prefix += weights[s];
            part.mass[rank] += weights[s];
            ++s;
        }

        part.end[rank] = s;
    }

    while (s < num_minima) {
        part.mass[p - 1] += weights[s];
        ++s;
    }

    part.end[p - 1] = num_minima;

    return part;
}

std::vector<long double>
BuildLocalMinCDF(const SInt min_begin, const SInt min_end, const std::vector<long double>& weights) {
    std::vector<long double> cdf;
    cdf.reserve(static_cast<std::size_t>(min_end - min_begin));

    long double prefix = 0.0L;

    for (SInt s = min_begin; s < min_end; ++s) {
        prefix += weights[s];
        cdf.push_back(prefix);
    }

    const long double total = cdf.back();

    for (auto& x: cdf) {
        x /= total;
    }

    cdf.back() = 1;
    return cdf;
}

std::vector<SInt>
DeterministicRankCounts(const SInt m, const std::vector<long double>& mass, const long double total_mass) {
    const SInt p = static_cast<SInt>(mass.size());

    std::vector<SInt>                         counts(p, 0);
    std::vector<std::pair<long double, SInt>> fractional;

    SInt assigned = 0;

    for (SInt rank = 0; rank < p; ++rank) {
        const long double exact = static_cast<long double>(m) * mass[rank] / total_mass;

        const SInt base = static_cast<SInt>(std::floor(exact));

        counts[rank] = base;
        assigned += base;

        fractional.emplace_back(exact - static_cast<long double>(base), rank);
    }

    std::sort(fractional.begin(), fractional.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

    const SInt remaining = m - assigned;

    for (SInt i = 0; i < remaining; ++i) {
        ++counts[fractional[i].second];
    }

    return counts;
}

long double LogBinomialApprox(const SInt n, const SInt k) {
    if (k > n) {
        return -std::numeric_limits<long double>::infinity();
    }

    return std::lgammal(static_cast<long double>(n) + 1) - std::lgammal(static_cast<long double>(k) + 1)
           - std::lgammal(static_cast<long double>(n - k) + 1);
}

bool ExpectedCountIsNegligible(const double expected_m_k) {
    return expected_m_k > 0.0 && expected_m_k < 1e-12;
}

long double LogDifferenceOfExponentials(const long double log_larger, const long double log_smaller) {
    if (log_smaller == -std::numeric_limits<long double>::infinity()) {
        return log_larger;
    }

    return log_larger + std::log1pl(-expl(log_smaller - log_larger));
}

long double MinRangeMassApprox(const SInt begin, const SInt end, const SInt n, const SInt k) {
    if (begin >= end) {
        return 0.0L;
    }

    const long double log_total = LogBinomialApprox(n, k);
    const long double log_a     = LogBinomialApprox(n - begin, k);

    const long double log_b =
        (n - end >= k) ? LogBinomialApprox(n - end, k) : -std::numeric_limits<long double>::infinity();

    const long double log_mass = LogDifferenceOfExponentials(log_a, log_b);

    return expl(log_mass - log_total);
}

SInt FindMinBoundaryByMass(const SInt n, const SInt k, const SInt rank, const SInt size) {
    const SInt num_minima = n - k + 1;

    if (rank <= 0) {
        return 0;
    }

    if (rank >= size) {
        return num_minima;
    }

    const long double target_tail = 1.0L - (static_cast<long double>(rank) / static_cast<long double>(size));

    const long double log_target_tail = std::log(target_tail);
    const long double log_total       = LogBinomialApprox(n, k);
    SInt              left            = 0;
    SInt              right           = num_minima;

    while (left < right) {
        const SInt        mid = left + ((right - left) / 2);
        const long double log_tail =
            (n - mid >= k) ? LogBinomialApprox(n - mid, k) - log_total : -std::numeric_limits<long double>::infinity();

        if (log_tail <= log_target_tail) {
            right = mid;
        } else {
            left = mid + 1;
        }
    }
    return left;
}

CountInt BinomialExact(SInt n, SInt k) {
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

CountInt MinRangeMassExact(const SInt begin, const SInt end, const SInt n, const SInt k) {
    if (begin >= end) {
        return 0;
    }

    const CountInt a = BinomialExact(n - begin, k);

    const CountInt b = (n - end >= k) ? BinomialExact(n - end, k) : CountInt(0);

    return a - b;
}

long double
MinRangeMassApproxCached(const SInt begin, const SInt end, const SInt n, const SInt k, LogBinomCache& cache) {
    if (begin >= end) {
        return 0.0L;
    }

    const long double log_total = cache.Get(n, k);
    const long double log_a     = cache.Get(n - begin, k);

    const long double log_b = (n - end >= k) ? cache.Get(n - end, k) : -std::numeric_limits<long double>::infinity();

    const long double log_mass = LogDifferenceOfExponentials(log_a, log_b);

    return expl(log_mass - log_total);
}

bool MinRangeMassDefinitelyExceedsSInt(const SInt begin, const SInt end, const SInt n, const SInt k) {
    if (begin >= end) {
        return false;
    }

    LogBinomCache cache(k);

    const long double log_total = std::log(static_cast<long double>(std::numeric_limits<SInt>::max()));

    const long double log_a = cache.Get(n - begin, k);

    const long double log_b = (n - end >= k) ? cache.Get(n - end, k) : -std::numeric_limits<long double>::infinity();

    const long double log_mass = LogDifferenceOfExponentials(log_a, log_b);

    return log_mass > log_total;
}

} // namespace kagen