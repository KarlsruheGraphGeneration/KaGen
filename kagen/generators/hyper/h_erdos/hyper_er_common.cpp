#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
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

void FloydSampleInto(
    const SInt universe_offset, const SInt universe_size, const SInt sample_size, Mersenne& mersenne,
    std::vector<SInt>& out, FloydScratchSet& scratch, const std::size_t offset) {
    if (sample_size > universe_size) {
        throw ConfigurationError("Cannot sample more pins than available vertices");
    }

    out.resize(offset);

    if (sample_size <= 32) {
        out.reserve(offset + static_cast<std::size_t>(sample_size));

        while (static_cast<SInt>(out.size() - offset) < sample_size) {
            const long double x = MersenneUniform01(mersenne);

            const SInt candidate = universe_offset + static_cast<SInt>(x * static_cast<long double>(universe_size));

            bool duplicate = false;

            for (std::size_t i = offset; i < out.size(); ++i) {
                if (out[i] == candidate) {
                    duplicate = true;
                    break;
                }
            }

            if (!duplicate) {
                out.push_back(candidate);
            }
        }

        std::sort(out.begin() + offset, out.end());

        return;
    }

    scratch.clear();

    scratch.reserve(static_cast<std::size_t>(sample_size));

    const SInt start = universe_size - sample_size;

    for (SInt j = start; j < universe_size; ++j) {
        const long double x = MersenneUniform01(mersenne);

        const SInt t = static_cast<SInt>(x * static_cast<long double>(j + 1));

        const SInt candidate = universe_offset + t;

        const SInt fallback = universe_offset + j;

        scratch.insert(scratch.contains(candidate) ? fallback : candidate);
    }

    out.reserve(offset + scratch.size());

    out.insert(out.end(), scratch.begin(), scratch.end());

    std::sort(out.begin() + offset, out.end());
}

std::vector<SInt>
FloydSample(const SInt universe_offset, const SInt universe_size, const SInt sample_size, Mersenne& mersenne) {
    std::vector<SInt> result;
    FloydScratchSet   scratch;
    FloydSampleInto(universe_offset, universe_size, sample_size, mersenne, result, scratch, 0);
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

SInt MersenneUniformSIntBelow(const SInt upper, Mersenne& mersenne) {
    if (upper == 0) {
        throw ConfigurationError("MersenneUniformSIntBelow requires upper > 0");
    }

    const auto bound = static_cast<std::uint64_t>(upper);

    const std::uint64_t threshold = static_cast<std::uint64_t>(-bound) % bound;

    while (true) {
        const std::uint64_t value =
            static_cast<std::uint64_t>(mersenne.BRandom()) | (static_cast<std::uint64_t>(mersenne.BRandom()) << 32);

        if (value >= threshold) {
            return static_cast<SInt>(value % bound);
        }
    }
}

/**
    Returns exact C(n, k) if it fits in SInt, else std::nullopt
*/
std::optional<SInt> TryBinomialSInt(SInt n, SInt k) {
    if (k > n) {
        return SInt{0};
    }

    k = std::min(k, n - k);

    unsigned __int128       result = 1;
    const unsigned __int128 limit  = std::numeric_limits<SInt>::max();

    for (SInt i = 1; i <= k; ++i) {
        unsigned __int128 numerator   = n - k + i;
        unsigned __int128 denominator = i;

        const auto g1 = std::gcd(static_cast<std::uint64_t>(numerator), static_cast<std::uint64_t>(denominator));

        numerator /= g1;
        denominator /= g1;

        const auto g2 =
            std::gcd(static_cast<std::uint64_t>(result % denominator), static_cast<std::uint64_t>(denominator));

        result /= g2;
        denominator /= g2;

        if (numerator != 0 && result > limit / numerator) {
            return std::nullopt;
        }

        result *= numerator;
        result /= denominator;

        if (result > limit) {
            return std::nullopt;
        }
    }

    return static_cast<SInt>(result);
}

template <typename BigInt>
ExactFixedCountHyperedgeGenerator<BigInt>::ExactFixedCountHyperedgeGenerator(
    const PGeneratorConfig& config, const PEID rank, const PEID size, RNGWrapper<>& rng, Mersenne& mersenne,
    Graph& graph, HypergraphMemoryStats& memory_stats, FloydScratchSet& floyd_scratch,
    ErdosHypergraphDebugLogger* debug_logger, SInt* next_debug_hyperedge_id
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    ,
    HGNPInstrumentation* instrumentation
#endif
    )
    : config_(config),
      rank_(rank),
      size_(size),
      rng_(rng),
      mersenne_(mersenne),
      graph_(graph),
      memory_stats_(memory_stats),
      floyd_scratch_(floyd_scratch),
      debug_logger_(debug_logger),
      next_debug_hyperedge_id_(next_debug_hyperedge_id)
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
      ,
      instrumentation_(instrumentation)
#endif
{
}

template <typename BigInt>
void ExactFixedCountHyperedgeGenerator<BigInt>::ValidateDuplicateCheckingIsFeasible(const SInt local_m) const {
    if (!config_.allow_duplicates && local_m > (SInt{1} << 26)) {
        throw ConfigurationError(
            "Duplicate checking for huge hypergraph generation "
            "is infeasible; use --fast");
    }
}

template <typename BigInt>
SInt ExactFixedCountHyperedgeGenerator<BigInt>::EdgeSeed(const SInt hyperedge_size) const {
    return sampling::Spooky::hash(
        static_cast<unsigned long long>(config_.seed)
        + (static_cast<unsigned long long>(rank_) * kEdgeRankSeedMultiplier)
        + (static_cast<unsigned long long>(hyperedge_size) * kEdgeSeedMultiplier));
}

template <typename BigInt>
SInt ExactFixedCountHyperedgeGenerator<BigInt>::RankSplitSeed(
    const SInt hyperedge_size, const SInt rank_begin, const SInt rank_end, const SInt level) const {
    return sampling::Spooky::hash(
        (static_cast<unsigned long long>(config_.seed) * kGNMCountSeedMultiplier)
        + (static_cast<unsigned long long>(hyperedge_size) * kCountSeedMultiplier)
        + (static_cast<unsigned long long>(rank_begin) * kRankSeedMultiplier)
        + (static_cast<unsigned long long>(rank_end) * kEdgeRankSeedMultiplier)
        + static_cast<unsigned long long>(level));
}

template <typename BigInt>
double ExactFixedCountHyperedgeGenerator<BigInt>::SplitProbability(
    const CountInt& child_population, const CountInt& parent_population) const {
    if (child_population <= 0) {
        return 0.0;
    }
    if (child_population >= parent_population) {
        return 1.0;
    }

    // This conversion is performed only once per level of the rank's
    // root-to-leaf path. Using a high-precision floating ratio avoids
    // magnitude-dependent native-integer dispatch and does not require
    // materializing or sampling from a machine-sized population.
    using ProbabilityFloat = boost::multiprecision::cpp_dec_float_50;

    const ProbabilityFloat numerator   = child_population.convert_to<ProbabilityFloat>();
    const ProbabilityFloat denominator = parent_population.convert_to<ProbabilityFloat>();
    const ProbabilityFloat probability = numerator / denominator;

    return std::clamp(probability.convert_to<double>(), 0.0, 1.0);
}

template <typename BigInt>
MinimumPartitionLocalCount
ExactFixedCountHyperedgeGenerator<BigInt>::ComputeLocalCount(const SInt hyperedge_size, const SInt global_m) {
    const CountInt total_population = BinomialExact(config_.n, hyperedge_size);

    if (CountInt(global_m) > total_population) {
        throw ConfigurationError("Fixed-count generation exceeds the hyperedge population");
    }

    if (size_ <= 0 || (size_ & (size_ - 1)) != 0) {
        throw ConfigurationError("Minimum-partition fixed-count generation requires a power-of-two PE count");
    }

    if (size_ == 1) {
        return {
            .local_m          = global_m,
            .local_population = total_population,
        };
    }

    const SInt num_minima = config_.n - hyperedge_size + 1;

    return ComputeLocalCountRecursive(hyperedge_size, 0, size_, rank_, 0, num_minima, total_population, global_m, 0);
}

template <typename BigInt>
MinimumPartitionLocalCount ExactFixedCountHyperedgeGenerator<BigInt>::ComputeLocalCountRecursive(
    const SInt hyperedge_size, const SInt rank_begin, const SInt rank_end, const SInt target_rank, const SInt min_begin,
    const SInt min_end, const CountInt& population, const SInt draws, const SInt level) {
    if (rank_end - rank_begin == 1) {
        return {
            .local_m          = draws,
            .local_population = population,
        };
    }

    const SInt rank_mid = rank_begin + ((rank_end - rank_begin) / 2);
    const SInt min_mid  = FindMinBoundaryByMass(config_.n, hyperedge_size, rank_mid, size_);

    const CountInt left_population  = MinRangeMassExact(min_begin, min_mid, config_.n, hyperedge_size);
    const CountInt right_population = population - left_population;

    if (left_population < 0 || right_population < 0) {
        throw ConfigurationError("Invalid minimum-partition population split");
    }

    SInt left_draws = 0;

    if (draws != 0 && left_population != 0) {
        if (right_population == 0) {
            left_draws = draws;
        } else {
            const double probability = SplitProbability(left_population, population);
            const SInt   seed        = RankSplitSeed(hyperedge_size, rank_begin, rank_end, level);

            // The split is pseudorandom and independently reproducible on every
            // descendant rank. It conserves the exact parent count because the
            // right child receives draws - left_draws.
            left_draws = rng_.GenerateBinomial(seed, draws, probability);

            // A binomial proposal may leave the finite-population support in
            // dense cases. Clamp to the exact feasible support, calculated from
            // the combinatorial interval populations.
            const auto [minimum_left, maximum_left] = HypergeometricSupport(population, left_population, draws);
            left_draws                              = std::clamp(left_draws, minimum_left, maximum_left);
        }
    }

    const SInt right_draws = draws - left_draws;

    if (CountInt(left_draws) > left_population || CountInt(right_draws) > right_population) {
        throw ConfigurationError("Minimum-partition rank split exceeds child population");
    }

    if (target_rank < rank_mid) {
        return ComputeLocalCountRecursive(
            hyperedge_size, rank_begin, rank_mid, target_rank, min_begin, min_mid, left_population, left_draws,
            level + 1);
    }

    return ComputeLocalCountRecursive(
        hyperedge_size, rank_mid, rank_end, target_rank, min_mid, min_end, right_population, right_draws, level + 1);
}

template <typename BigInt>
void ExactFixedCountHyperedgeGenerator<BigInt>::ValidateExactSparseDensity(
    const SInt local_m, const CountInt& local_population) const {
    if (CountInt(local_m) > local_population) {
        throw ConfigurationError(
            "Exact fixed-count local edge count exceeds "
            "the local hyperedge population");
    }

    if (local_population == 0) {
        return;
    }

    if (CountInt(local_m) * 4 > local_population) {
        throw ConfigurationError(
            "Dense exact fixed-count hypergraph generation "
            "is not implemented");
    }
}

template <typename BigInt>
FixedCountLocalRange
ExactFixedCountHyperedgeGenerator<BigInt>::PrepareLocalRange(const SInt hyperedge_size, const SInt global_m) {
    if (hyperedge_size == 0 || hyperedge_size > config_.n || global_m == 0) {
        return {};
    }

    // Minimum-vertex ownership remains the work partition. Each rank needs
    // only its own two boundaries and its O(log P) deterministic split path.
    const SInt min_begin = FindMinBoundaryByMass(config_.n, hyperedge_size, rank_, size_);
    const SInt min_end   = FindMinBoundaryByMass(config_.n, hyperedge_size, rank_ + 1, size_);

    const MinimumPartitionLocalCount local = ComputeLocalCount(hyperedge_size, global_m);

    ValidateDuplicateCheckingIsFeasible(local.local_m);
    ValidateExactSparseDensity(local.local_m, local.local_population);

    return {
        .min_begin = min_begin,
        .min_end   = min_end,
        .local_m   = local.local_m,
    };
}

template <typename BigInt>
void ExactFixedCountHyperedgeGenerator<BigInt>::SampleHyperedgeInto(
    const SInt minimum_vertex, const SInt hyperedge_size, std::vector<SInt>& pins) {
    pins.clear();
    pins.push_back(minimum_vertex);

    FloydSampleInto(
        minimum_vertex + 1, config_.n - minimum_vertex - 1, hyperedge_size - 1, mersenne_, pins, floyd_scratch_, 1);
}

template <typename BigInt>
bool ExactFixedCountHyperedgeGenerator<BigInt>::TryPushHyperedge(
    const std::vector<SInt>& pins, HyperedgeSeenSet& local_seen) {
    bool accepted = true;

    if (!config_.allow_duplicates) {
        const double duplicate_start = MPI_Wtime();

        accepted = local_seen.insert(FingerprintPins(pins)).second;

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION

        if (instrumentation_ != nullptr) {
            instrumentation_->duplicate_check_seconds += MPI_Wtime() - duplicate_start;
        }
    }

    if (!accepted) {
        if (instrumentation_ != nullptr) {
            ++instrumentation_->duplicate_rejects;
        }

        return false;
#endif
    }
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    const double write_start = MPI_Wtime();
#endif
    PushUncompressedHyperedge(graph_, memory_stats_, pins);
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    if (instrumentation_ != nullptr) {
        instrumentation_->csr_write_seconds += MPI_Wtime() - write_start;

        ++instrumentation_->generated_edges;

        instrumentation_->generated_pins += static_cast<std::uint64_t>(pins.size());
    }
#endif
    return true;
}

template <typename BigInt>
void ExactFixedCountHyperedgeGenerator<BigInt>::GenerateSampledHyperedges(
    const SInt hyperedge_size, const SInt local_min_begin, const SInt local_min_end, const SInt local_m,
    SInt& edge_seed, LogBinomCache& log_binom_cache) {
    ValidateDuplicateCheckingIsFeasible(local_m);

    std::vector<SInt> pins;
    pins.reserve(static_cast<std::size_t>(hyperedge_size));

    auto local_seen = MakeLocalSeenSet(config_.allow_duplicates, local_m);

    const std::uint64_t local_m_u64 = static_cast<std::uint64_t>(local_m);

    const std::uint64_t max_attempts = local_m_u64 > std::numeric_limits<std::uint64_t>::max() / 10
                                           ? std::numeric_limits<std::uint64_t>::max()
                                           : std::max<std::uint64_t>(local_m_u64 * 10, 1000);

    std::uint64_t total_attempts = 0;
    SInt          generated      = 0;

    while (generated < local_m) {
        const auto event_start = std::chrono::steady_clock::now();

        std::uint64_t sampling_attempts    = 0;
        std::uint64_t duplicate_rejections = 0;
        std::uint64_t minimum_search_steps = 0;
        std::uint64_t minimum_cache_gets   = 0;

        while (true) {
            std::uint64_t attempt_minimum_steps = 0;
            std::uint64_t attempt_cache_gets    = 0;

            SInt minimum_vertex;

            {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
                ScopedMPITimer timer(instrumentation_->minimum_sample_seconds);
#endif
                minimum_vertex = SampleMinimumImplicit(
                    local_min_begin, local_min_end, config_.n, hyperedge_size, rng_, edge_seed, log_binom_cache,
                    &attempt_minimum_steps, &attempt_cache_gets);
            }

            ++sampling_attempts;
            ++total_attempts;

            minimum_search_steps += attempt_minimum_steps;
            minimum_cache_gets += attempt_cache_gets;

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
            if (instrumentation_ != nullptr) {
                ++instrumentation_->minimum_samples;
                instrumentation_->minimum_search_steps += attempt_minimum_steps;
                instrumentation_->minimum_cache_gets += attempt_cache_gets;
            }
#endif

            {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
                ScopedMPITimer timer(instrumentation_->pin_sample_seconds);
#endif
                SampleHyperedgeInto(minimum_vertex, hyperedge_size, pins);
            }

            if (!TryPushHyperedge(pins, local_seen)) {
                ++duplicate_rejections;

                if (total_attempts > max_attempts) {
                    throw ConfigurationError("Exact fixed-count rejection sampling exceeded the attempt limit");
                }

                continue;
            }

            ++generated;

            if (debug_logger_ != nullptr && next_debug_hyperedge_id_ != nullptr) {
                const auto duration_ns =
                    std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - event_start)
                        .count();

                debug_logger_->LogHyperedge({
                    .hyperedge_id         = (*next_debug_hyperedge_id_)++,
                    .hyperedge_size       = hyperedge_size,
                    .minimum_vertex       = pins.front(),
                    .sampling_attempts    = static_cast<SInt>(sampling_attempts),
                    .duplicate_rejections = static_cast<SInt>(duplicate_rejections),
                    .minimum_search_steps = static_cast<SInt>(minimum_search_steps),
                    .minimum_cache_gets   = static_cast<SInt>(minimum_cache_gets),
                    .duration_ns          = duration_ns,
                });
            }

            break;
        }
    }

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    if (instrumentation_ != nullptr) {
        instrumentation_->attempts += total_attempts;
    }
#endif
}

template <typename BigInt>
void ExactFixedCountHyperedgeGenerator<BigInt>::Generate(const SInt hyperedge_size, const FixedCountLocalRange& range) {
    if (hyperedge_size == 0 || hyperedge_size > config_.n || range.local_m == 0) {
        return;
    }

    if (range.min_begin >= range.min_end) {
        throw ConfigurationError(
            "Nonzero local edge count assigned to an "
            "empty minimum range");
    }

    const double cache_start = MPI_Wtime();

    const std::size_t static_hint = ComputeLogBinomCacheSize(range.local_m, range.min_begin, range.min_end);

    const std::size_t cache_size = std::min(static_hint, log_binom_cache_reserve_hint_);

    LogBinomCache log_binom_cache(hyperedge_size, cache_size);

    log_binom_cache.InitializeRange(config_.n, range.min_begin, range.min_end);
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    if (instrumentation_ != nullptr) {
        instrumentation_->cache_init_seconds += MPI_Wtime() - cache_start;

        instrumentation_->max_cache_initial_reserve = std::max<std::uint64_t>(
            instrumentation_->max_cache_initial_reserve, static_cast<std::uint64_t>(cache_size));
    }
#endif

    SInt edge_seed = EdgeSeed(hyperedge_size);

    mersenne_.RandomInit(edge_seed);

    GenerateSampledHyperedges(
        hyperedge_size, range.min_begin, range.min_end, range.local_m, edge_seed, log_binom_cache);

    const std::size_t observed_size = 0;

    constexpr std::size_t kMinimumHint = 4096;
    constexpr std::size_t kMaximumHint = 1U << 17;

    const std::size_t observed_with_headroom =
        observed_size > std::numeric_limits<std::size_t>::max() / 3 ? kMaximumHint : ((observed_size * 3) + 1) / 2;

    log_binom_cache_reserve_hint_ = std::clamp(observed_with_headroom, kMinimumHint, kMaximumHint);
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    // Put the instrumentation block here.
    if (instrumentation_ != nullptr) {
        const LogBinomCacheStats& stats = log_binom_cache.stats;

        instrumentation_->cache_map_calls += stats.map_calls;

        instrumentation_->cache_map_hits += stats.map_hits;

        instrumentation_->cache_map_misses += stats.map_misses;

        instrumentation_->cache_map_inserts += stats.map_inserts;

        instrumentation_->cache_candidate_calls += stats.candidate_calls;

        instrumentation_->cache_candidate_exact_hits += stats.candidate_exact_hits;

        instrumentation_->cache_candidate_recurrence_hits += stats.candidate_recurrence_hits;

        instrumentation_->cache_candidate_direct_evals += stats.candidate_direct_evals;

        instrumentation_->cache_candidate_below_range += stats.candidate_below_range;

        instrumentation_->cache_candidate_distance_sum += stats.candidate_distance_sum;

        instrumentation_->cache_candidate_max_distance =
            std::max(instrumentation_->cache_candidate_max_distance, stats.candidate_max_distance);

        instrumentation_->cache_backward_steps += stats.backward_steps;

        instrumentation_->cache_forward_steps += stats.forward_steps;

        instrumentation_->cache_max_size = std::max(instrumentation_->cache_max_size, stats.max_size);

        const std::uint64_t total_calls = stats.map_calls + stats.candidate_calls;

        if (total_calls > 0) {
            instrumentation_->cache_min_key = std::min(instrumentation_->cache_min_key, stats.min_key);

            instrumentation_->cache_max_key = std::max(instrumentation_->cache_max_key, stats.max_key);
        }
    }
#endif
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"

template class ExactFixedCountHyperedgeGenerator<SInt>;
template class ExactFixedCountHyperedgeGenerator<__int128>;

#pragma GCC diagnostic pop

} // namespace kagen