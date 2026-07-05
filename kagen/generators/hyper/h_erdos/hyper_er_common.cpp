#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace kagen {

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

long double LogAdd(long double a, long double b) {
    if (a == -std::numeric_limits<long double>::infinity()) {
        return b;
    }
    if (b == -std::numeric_limits<long double>::infinity()) {
        return a;
    }
    if (a < b) {
        std::swap(a, b);
    }
    return a + std::log1pl(expl(b - a));
}

long double BinomSmallExactLD(SInt n, SInt q) {
    if (q < 0 || q > n) {
        return 0.0L;
    }

    q = std::min(q, n - q);

    long double result = 1.0L;
    for (SInt t = 1; t <= q; ++t) {
        result *= static_cast<long double>(n - q + t);
        result /= static_cast<long double>(t);
    }
    return result;
}

long double MersenneUniform01(Mersenne& mersenne) {
    return std::min<long double>(static_cast<long double>(mersenne.Random()), std::nextafter(1.0L, 0.0L));
}
std::vector<SInt>
FloydSample(const SInt universe_offset, const SInt universe_size, const SInt sample_size, Mersenne& mersenne) {
    if (sample_size > universe_size) {
        throw ConfigurationError("Cannot sample more pins than available vertices");
    }

    if (sample_size <= 32) {
        std::vector<SInt> result;
        result.reserve(sample_size);

        while (static_cast<SInt>(result.size()) < sample_size) {
            const long double x = MersenneUniform01(mersenne);

            const SInt candidate = universe_offset + static_cast<SInt>(x * static_cast<long double>(universe_size));

            bool duplicate = false;
            for (const SInt v: result) {
                if (v == candidate) {
                    duplicate = true;
                    break;
                }
            }

            if (!duplicate) {
                result.push_back(candidate);
            }
        }

        std::sort(result.begin(), result.end());
        return result;
    }

    std::unordered_set<SInt> selected;
    selected.reserve(static_cast<std::size_t>(sample_size));

    const SInt start = universe_size - sample_size;

    for (SInt j = start; j < universe_size; ++j) {
        const long double x = MersenneUniform01(mersenne);

        const SInt t = static_cast<SInt>(x * static_cast<long double>(j + 1));

        const SInt candidate = universe_offset + t;
        const SInt fallback  = universe_offset + j;

        selected.insert(selected.contains(candidate) ? fallback : candidate);
    }

    std::vector<SInt> result(selected.begin(), selected.end());
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<double> ReadNumericLines(const std::string& filename, const std::string& error_message) {
    std::ifstream in(filename);
    if (!in) {
        throw ConfigurationError(error_message);
    }

    std::vector<double> values;
    std::string         line;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        values.push_back(std::stod(line));
    }

    return values;
}

} // namespace kagen