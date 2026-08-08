#pragma once

#include "kagen/generators/generator.h"
#include "kagen/hypergraph/debug_logger_erdos.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"
#include "kagen/tools/mersenne.h"
#include "kagen/tools/rng_wrapper.h"

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_flat_set.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <limits>
#include <vector>

namespace kagen {

using CountInt         = boost::multiprecision::cpp_int;
using FloydScratchSet  = boost::unordered_flat_set<SInt>;
using HyperedgeSeenSet = boost::unordered_flat_set<std::uint64_t>;

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

inline std::uint64_t FingerprintPins(const std::vector<SInt>& pins) {
    std::uint64_t hash_value = kFNVOffsetBasis;

    for (const SInt value: pins) {
        hash_value ^= static_cast<std::uint64_t>(value);
        hash_value *= kFNVPrime;
    }

    return hash_value;
}

long double LogBinomialApprox(SInt n, SInt k);

struct LogBinomCacheStats {
    // Map-backed Get().
    std::uint64_t map_calls   = 0;
    std::uint64_t map_hits    = 0;
    std::uint64_t map_misses  = 0;
    std::uint64_t map_inserts = 0;

    // Cursor-backed GetCandidate().
    std::uint64_t candidate_calls           = 0;
    std::uint64_t candidate_exact_hits      = 0;
    std::uint64_t candidate_recurrence_hits = 0;
    std::uint64_t candidate_direct_evals    = 0;
    std::uint64_t candidate_below_range     = 0;

    std::uint64_t candidate_distance_sum          = 0;
    std::uint64_t candidate_distance_observations = 0;
    std::uint64_t candidate_max_distance          = 0;

    std::uint64_t backward_steps = 0;
    std::uint64_t forward_steps  = 0;

#ifdef KAGEN_USE_LOGBINOM_CURSOR_CACHE
    std::uint64_t get_exact_hits            = 0;
    std::uint64_t get_recurrence_hits       = 0;
    std::uint64_t get_direct_evals          = 0;
    std::uint64_t get_distance_sum          = 0;
    std::uint64_t get_distance_observations = 0;
    std::uint64_t get_max_distance          = 0;
    std::uint64_t get_forward_steps         = 0;
    std::uint64_t get_backward_steps        = 0;
#endif

    // Across both access methods.
    SInt min_key = std::numeric_limits<SInt>::max();

    SInt max_key = 0;

    std::uint64_t max_size = 0;
};

struct LogBinomCache {
    SInt fixed_k = -1;

    long double log_k_factorial = 0.0L;
    long double inv_k;

    // Fixed rank-local interval tails.
    bool        range_initialized = false;
    SInt        range_n           = 0;
    SInt        range_begin       = 0;
    SInt        range_end         = 0;
    long double log_begin_tail    = -std::numeric_limits<long double>::infinity();
    long double log_end_tail      = -std::numeric_limits<long double>::infinity();

    // Moving candidate cursor.
    bool        cursor_initialized = false;
    SInt        cursor_x           = 0;
    long double cursor_value       = -std::numeric_limits<long double>::infinity();

#ifdef KAGEN_USE_LOGBINOM_CURSOR_CACHE
    bool        get_cursor_initialized = false;
    SInt        get_cursor_x           = 0;
    long double get_cursor_value       = -std::numeric_limits<long double>::infinity();
#endif

    LogBinomCacheStats stats;

    explicit LogBinomCache(const SInt k = -1, [[maybe_unused]] const std::size_t expected_size = 4096) : fixed_k(k) {
        if (fixed_k >= 0) {
            log_k_factorial = std::lgammal(static_cast<long double>(fixed_k) + 1.0L);

            inv_k = 1.0L / static_cast<long double>(fixed_k);
        }
    }

    long double EvaluateDirect(const SInt x) const {
        if (x < fixed_k) {
            return -std::numeric_limits<long double>::infinity();
        }

        if (fixed_k <= 32) {
            return LogBinomialSmallKFast(x);
        }

        return LogBinomialApprox(x, fixed_k);
    }

    long double LogBinomialSmallKFast(const SInt n) const {
        long double result = -log_k_factorial;

        for (SInt i = 0; i < fixed_k; ++i) {
            result += std::log(static_cast<long double>(n - i));
        }

        return result;
    }

    long double LogBinomialRatioSmallK(const SInt x, const SInt anchor_x) const {
        assert(fixed_k >= 0);
        assert(x >= fixed_k);
        assert(anchor_x >= fixed_k);

        /*
         * log(C(x,k) / C(anchor_x,k))
         *
         *   = sum_i log((x-i)/(anchor_x-i))
         *   = sum_i log(1 + (x-anchor_x)/(anchor_x-i)).
         *
         * In the minimum sampler x <= anchor_x, so each term is
         * nonpositive.
         */
        const long double difference = static_cast<long double>(x - anchor_x);

        long double result = 0.0L;

        for (SInt i = 0; i < fixed_k; ++i) {
            const long double denominator = static_cast<long double>(anchor_x - i);

            result += std::log1pl(difference / denominator);
        }

        return result;
    }

    long double StepBackward(const SInt x, const long double value) const {
        assert(x >= fixed_k);
        assert(x > 0);

        // C(k - 1, k) = 0.
        if (x == fixed_k) {
            return -std::numeric_limits<long double>::infinity();
        }

        // log C(x - 1, k)
        // = log C(x, k) + log((x - k) / x)
        return value + std::log1pl(-static_cast<long double>(fixed_k) / static_cast<long double>(x));
    }

    long double StepForward(const SInt x, const long double value) const {
        assert(x >= fixed_k);

        // log C(x + 1, k)
        // = log C(x, k) + log((x + 1) / (x + 1 - k))
        return value + std::log1pl(static_cast<long double>(fixed_k) / static_cast<long double>(x + 1 - fixed_k));
    }

    long double GetCandidate(const SInt target_x) {
        constexpr SInt kMaxRecurrenceDistance = 16;

        ++stats.candidate_calls;

        stats.min_key = std::min(stats.min_key, target_x);
        stats.max_key = std::max(stats.max_key, target_x);

        if (target_x < fixed_k) {
            ++stats.candidate_below_range;
            return -std::numeric_limits<long double>::infinity();
        }

        if (!cursor_initialized) {
            cursor_value       = EvaluateDirect(target_x);
            cursor_x           = target_x;
            cursor_initialized = true;

            ++stats.candidate_direct_evals;
            return cursor_value;
        }

        if (target_x == cursor_x) {
            ++stats.candidate_exact_hits;
            return cursor_value;
        }

        const SInt distance = target_x > cursor_x ? target_x - cursor_x : cursor_x - target_x;

        ++stats.candidate_distance_observations;

        stats.candidate_distance_sum += static_cast<std::uint64_t>(distance);

        stats.candidate_max_distance =
            std::max<std::uint64_t>(stats.candidate_max_distance, static_cast<std::uint64_t>(distance));

        if (distance > kMaxRecurrenceDistance) {
            cursor_value = EvaluateDirect(target_x);
            cursor_x     = target_x;

            ++stats.candidate_direct_evals;
            return cursor_value;
        }

        ++stats.candidate_recurrence_hits;

        while (cursor_x < target_x) {
            cursor_value = StepForward(cursor_x, cursor_value);

            ++cursor_x;
            ++stats.forward_steps;
        }

        while (cursor_x > target_x) {
            cursor_value = StepBackward(cursor_x, cursor_value);

            --cursor_x;
            ++stats.backward_steps;
        }

        return cursor_value;
    }

    void SetCandidateCursor(const SInt x, const long double value) {
        cursor_x           = x;
        cursor_value       = value;
        cursor_initialized = true;
    }

    long double Get(const SInt x, const SInt k) {
        assert(k == fixed_k);

        ++stats.map_calls;

        stats.min_key = std::min(stats.min_key, x);
        stats.max_key = std::max(stats.max_key, x);

#ifdef KAGEN_USE_LOGBINOM_CURSOR_CACHE
        return GetCached(x);
#else
        ++stats.map_misses;
        return EvaluateDirect(x);
#endif
    }
#ifdef KAGEN_USE_LOGBINOM_CURSOR_CACHE
    long double GetCached(const SInt target_x) {
        constexpr SInt kMaxRecurrenceDistance = 16;

        if (target_x < fixed_k) {
            ++stats.map_misses;
            return -std::numeric_limits<long double>::infinity();
        }

        if (!get_cursor_initialized) {
            get_cursor_value       = EvaluateDirect(target_x);
            get_cursor_x           = target_x;
            get_cursor_initialized = true;

            ++stats.map_misses;
            ++stats.get_direct_evals;
            return get_cursor_value;
        }

        if (target_x == get_cursor_x) {
            ++stats.map_hits;
            ++stats.get_exact_hits;
            return get_cursor_value;
        }

        const SInt distance = target_x > get_cursor_x ? target_x - get_cursor_x : get_cursor_x - target_x;

        ++stats.get_distance_observations;
        stats.get_distance_sum += static_cast<std::uint64_t>(distance);
        stats.get_max_distance = std::max(stats.get_max_distance, static_cast<std::uint64_t>(distance));

        if (distance > kMaxRecurrenceDistance) {
            get_cursor_value = EvaluateDirect(target_x);
            get_cursor_x     = target_x;

            ++stats.map_misses;
            ++stats.get_direct_evals;
            return get_cursor_value;
        }

        ++stats.map_hits;
        ++stats.get_recurrence_hits;

        while (get_cursor_x < target_x) {
            get_cursor_value = StepForward(get_cursor_x, get_cursor_value);
            ++get_cursor_x;
            ++stats.get_forward_steps;
        }

        while (get_cursor_x > target_x) {
            get_cursor_value = StepBackward(get_cursor_x, get_cursor_value);
            --get_cursor_x;
            ++stats.get_backward_steps;
        }

        return get_cursor_value;
    }
#endif
    void InitializeRange(const SInt n, const SInt local_begin, const SInt local_end) {
        if (local_begin >= local_end) {
            throw ConfigurationError("LogBinomCache requires a nonempty local minimum range");
        }

        range_n     = n;
        range_begin = local_begin;
        range_end   = local_end;

        log_begin_tail = Get(n - local_begin, fixed_k);

        log_end_tail =
            n - local_end >= fixed_k ? Get(n - local_end, fixed_k) : -std::numeric_limits<long double>::infinity();

        range_initialized = true;
    }
};

struct HyperERPartitionRange {
    SInt begin;
    SInt end;
};

HyperERPartitionRange AssignPartitionsToPE(SInt num_partitions, PEID rank, PEID size);

inline std::size_t ComputeLogBinomCacheSize(const SInt local_m, const SInt local_min_begin, const SInt local_min_end) {
    constexpr std::size_t kMaximumInitialSize = 1U << 16;

    if (local_m <= 0 || local_min_begin >= local_min_end) {
        return 0;
    }

    const auto search_width = static_cast<std::uint64_t>(local_min_end - local_min_begin);

    return static_cast<std::size_t>(std::min<std::uint64_t>(search_width, kMaximumInitialSize));
}

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

SInt MersenneUniformSIntBelow(SInt upper, Mersenne& mersenne);

long double BinomSmallExactLD(SInt n, SInt q);

std::optional<SInt> TryBinomialSInt(SInt n, SInt k);

std::vector<SInt> FloydSample(SInt universe_offset, SInt universe_size, SInt sample_size, Mersenne& mersenne);
void              FloydSampleInto(
    SInt universe_offset, SInt universe_size, SInt sample_size, Mersenne& mersenne, std::vector<SInt>& out,
    FloydScratchSet& scratch, std::size_t offset);

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

class LogBinomialCorrectionState {
public:
    LogBinomialCorrectionState(
        const SInt candidate, const SInt n, const SInt k, const long double candidate_residual,
        const long double k_tail_residual, const SInt local_begin, const SInt local_end)
        : candidate_(candidate),
          candidate_x_(n - (candidate + 1)),
          n_(n),
          k_(k),
          candidate_residual_(candidate_residual),
          k_tail_residual_(k_tail_residual),
          local_begin_(local_begin),
          local_end_(local_end) {
        AssertInvariant();
    }

    long double CandidateResidual() const {
        return candidate_residual_;
    }

    bool HasPredecessor() const {
        return candidate_ > local_begin_;
    }

    long double CurrentPredecessorResidual() const {
        assert(HasPredecessor());
        return PredecessorResidual();
    }

    SInt Candidate() const {
        return candidate_;
    }

    bool CandidateIsAtOrRight() const {
        return candidate_residual_ <= 0.0L;
    }

    bool PredecessorIsAtOrRight() const {
        if (candidate_ == local_begin_) {
            return false;
        }

        return PredecessorResidual() <= 0.0L;
    }

    void MoveLeft() {
        assert(candidate_ > local_begin_);

        candidate_residual_ = PredecessorResidual();

        --candidate_;
        ++candidate_x_;

        AssertInvariant();
    }

    void MoveRight() {
        assert(candidate_ < local_end_ - 1);

        candidate_residual_ = SuccessorResidual();

        ++candidate_;
        --candidate_x_;

        AssertInvariant();
    }

private:
    void AssertInvariant() const {
        assert(local_begin_ < local_end_);
        assert(candidate_ >= local_begin_);
        assert(candidate_ < local_end_);
        assert(candidate_x_ == n_ - (candidate_ + 1));
    }

    long double PredecessorResidual() const {
        if (candidate_x_ < k_) {
            /*
             * candidate_x_ == k - 1.
             * The predecessor has x == k and C(k,k) == 1.
             *
             * Its normalized target residual was computed directly
             * by the caller.
             */
            return k_tail_residual_;
        }

        /*
         * Moving to the predecessor changes x to x + 1:
         *
         * log C(x+1,k)
         *   = log C(x,k)
         *     + log(1 + k/(x+1-k)).
         *
         * The anchor and target terms cancel, so the same recurrence
         * applies to the residual.
         */
        return candidate_residual_
               + std::log1pl(static_cast<long double>(k_) / static_cast<long double>(candidate_x_ + 1 - k_));
    }

    long double SuccessorResidual() const {
        if (candidate_x_ <= k_) {
            /*
             * Moving from x == k to x == k - 1 gives C(x,k) == 0.
             */
            return -std::numeric_limits<long double>::infinity();
        }

        /*
         * Moving to the successor changes x to x - 1:
         *
         * log C(x-1,k)
         *   = log C(x,k)
         *     + log(1 - k/x).
         */
        return candidate_residual_
               + std::log1pl(-static_cast<long double>(k_) / static_cast<long double>(candidate_x_));
    }

    SInt candidate_;
    SInt candidate_x_;
    SInt n_;
    SInt k_;

    long double candidate_residual_;

    /*
     * Residual for x == k:
     *
     * log(C(k,k)/begin_tail) - log(target/begin_tail).
     */
    long double k_tail_residual_;

    SInt local_begin_;
    SInt local_end_;
};
struct CorrectedMinimumCandidate {
    SInt        candidate;
    long double residual;

    bool        has_predecessor;
    long double predecessor_residual;
};

inline CorrectedMinimumCandidate CorrectMinimumCandidateByRecurrence(
    const SInt initial_candidate, const SInt local_begin, const SInt local_end, const SInt n, const SInt k,
    const long double candidate_residual, const long double k_tail_residual, std::uint64_t* steps) {
    constexpr SInt kMaxLinearCorrections = 8;

    LogBinomialCorrectionState state(
        initial_candidate, n, k, candidate_residual, k_tail_residual, local_begin, local_end);

    SInt corrections = 0;

    while (state.Candidate() > local_begin && state.PredecessorIsAtOrRight()) {
        state.MoveLeft();
        ++corrections;

        if (steps != nullptr) {
            ++(*steps);
        }

        if (corrections >= kMaxLinearCorrections) {
            break;
        }
    }

    while (corrections < kMaxLinearCorrections && state.Candidate() < local_end - 1 && !state.CandidateIsAtOrRight()) {
        state.MoveRight();
        ++corrections;

        if (steps != nullptr) {
            ++(*steps);
        }
    }

    const bool has_predecessor = state.HasPredecessor();

    return {
        .candidate = state.Candidate(),
        .residual  = state.CandidateResidual(),

        .has_predecessor = has_predecessor,

        .predecessor_residual =
            has_predecessor ? state.CurrentPredecessorResidual() : -std::numeric_limits<long double>::infinity(),
    };
}

template <typename Predicate>
inline SInt CorrectMinimumCandidateDirect(
    SInt candidate, const SInt local_begin, const SInt local_end, Predicate&& at_or_right_of_answer,
    std::uint64_t* steps) {
    constexpr SInt kMaxLinearCorrections = 8;

    SInt corrections = 0;

    while (candidate > local_begin && at_or_right_of_answer(candidate - 1)) {
        --candidate;
        ++corrections;

        if (steps != nullptr) {
            ++(*steps);
        }

        if (corrections >= kMaxLinearCorrections) {
            break;
        }
    }

    while (corrections < kMaxLinearCorrections && candidate < local_end - 1 && !at_or_right_of_answer(candidate)) {
        ++candidate;
        ++corrections;

        if (steps != nullptr) {
            ++(*steps);
        }
    }

    return candidate;
}

template <typename Predicate>
inline SInt ResolveMinimumCandidate(
    const SInt candidate, const SInt local_begin, const SInt local_end, Predicate&& at_or_right_of_answer,
    std::uint64_t* steps) {
    /*
     * Verify that candidate is the first position satisfying the monotone
     * predicate.
     */
    const bool candidate_valid = at_or_right_of_answer(candidate);

    const bool predecessor_invalid = candidate == local_begin || !at_or_right_of_answer(candidate - 1);

    if (candidate_valid && predecessor_invalid) {
        return candidate;
    }

    std::cerr << "Expensive minOwner[ candidate:" << candidate << " , local_begin:" << local_begin
              << " , local_end:" << local_end;

    /*
     * Restrict the fallback search to the unresolved side of candidate.
     */
    SInt left;
    SInt right;

    if (!candidate_valid) {
        /*
         * Candidate lies strictly before the answer.
         */
        left  = candidate + 1;
        right = local_end - 1;
    } else {
        /*
         * Candidate satisfies the predicate, but its predecessor does too,
         * so the answer lies strictly to the left.
         */
        left  = local_begin;
        right = candidate - 1;
    }

    while (left < right) {
        if (steps != nullptr) {
            ++(*steps);
        }

        const SInt mid = left + ((right - left) / 2);

        if (at_or_right_of_answer(mid)) {
            right = mid;
        } else {
            left = mid + 1;
        }
    }

    std::cerr << " , resolved:" << left
              << " , abs(candidate-resolved):" << (candidate > left ? candidate - left : left - candidate) << "]\n";
    return left;
}

inline SInt EstimateMinimumCandidate(
    const SInt local_begin, const SInt local_end, const SInt n, const SInt k, const long double log_target_tail,
    const long double log_k_factorial, const long double inv_k) {
    SInt candidate = local_begin;

    const bool estimate_is_usable = std::isfinite(log_target_tail) && std::isfinite(log_k_factorial);

    if (!estimate_is_usable) {
        return candidate;
    }

    const long double inverse_exponent = (log_target_tail + log_k_factorial) * inv_k;

    const long double base_estimate = std::exp(inverse_exponent);

    const long double x_estimate = base_estimate + (0.5L * static_cast<long double>(k - 1));

    const long double s_estimate = static_cast<long double>(n) - 1.0L - x_estimate;

    /*
     * Check floating-point bounds before converting to SInt.
     */
    if (!std::isfinite(s_estimate) || s_estimate <= static_cast<long double>(local_begin)) {
        return local_begin;
    }

    if (s_estimate >= static_cast<long double>(local_end - 1)) {
        return local_end - 1;
    }

    return static_cast<SInt>(std::floor(s_estimate));
}

template <typename RNG>
SInt SampleMinimumImplicitK2(const SInt local_begin, const SInt local_end, const SInt n, RNG& rng, SInt& seed) {
    const SInt draw_seed = sampling::Spooky::hash(seed++);

    const long double u = std::min<long double>(
        static_cast<long double>(rng.GenerateUniform(draw_seed, 0.0L, 1.0L)), std::nextafter(1.0L, 0.0L));

    const auto choose2 = [](const long double x) {
        return (x * (x - 1.0L)) / 2.0L;
    };

    const long double begin_tail = choose2(static_cast<long double>(n - local_begin));
    const long double end_tail   = (n - local_end >= 2) ? choose2(static_cast<long double>(n - local_end)) : 0.0L;

    const long double target_tail = begin_tail - (u * (begin_tail - end_tail));

    long double x_real = (1.0L + std::sqrt(1.0L + (8.0L * target_tail))) / 2.0L;

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

    if (local_begin >= local_end) {
        throw ConfigurationError("SampleMinimumImplicit requires a nonempty local minimum range");
    }

    /*
     * Draw a target from the exact rank-local tail interval.
     *
     * For a minimum S,
     *
     *     tail(s) = C(n - s, k).
     *
     * The returned minimum is the smallest s for which
     *
     *     C(n - (s + 1), k) <= target_tail.
     */
    const SInt draw_seed = sampling::Spooky::hash(seed++);

    const long double u = std::min<long double>(
        static_cast<long double>(rng.GenerateUniform(draw_seed, 0.0L, 1.0L)), std::nextafter(1.0L, 0.0L));

    auto cached_log_binomial = [&](const SInt x) -> long double {
        if (cache_gets) {
            ++(*cache_gets);
        }

        return cache.Get(x, k);
    };

    assert(cache.range_initialized);
    assert(cache.fixed_k == k);
    assert(cache.range_n == n);
    assert(cache.range_begin == local_begin);
    assert(cache.range_end == local_end);

    const long double log_begin_tail = cache.log_begin_tail;
    const long double log_end_tail   = cache.log_end_tail;

    /*
     * Normalize all tails by begin_tail.
     *
     * end_ratio = end_tail / begin_tail.
     *
     * local_mass_fraction
     *     = (begin_tail - end_tail) / begin_tail
     *     = 1 - end_ratio.
     */
    const long double local_mass_fraction =
        std::isinf(log_end_tail) && log_end_tail < 0.0L ? 1.0L : -std::expm1l(log_end_tail - log_begin_tail);

    /*
     * target_tail
     *     = begin_tail
     *       - uniform01 * (begin_tail - end_tail)
     *
     * Therefore
     *
     * target_tail / begin_tail
     *     = 1 - uniform01 * local_mass_fraction.
     */
    const long double log_target_ratio = std::log1pl(-u * local_mass_fraction);

    const SInt begin_x = n - local_begin;

    auto normalized_log_binomial = [&](const SInt x) -> long double {
        if (cache_gets != nullptr) {
            ++(*cache_gets);
        }

        return cache.LogBinomialRatioSmallK(x, begin_x);
    };
    /*
     * Monotone predicate used by both local correction and binary search.
     *
     * false: s is strictly before the selected minimum.
     * true:  s is the selected minimum or lies to its right.
     */
    auto at_or_right_of_answer = [&](const SInt s) -> bool {
        const SInt remaining = n - (s + 1);

        if (remaining < k) {
            return true;
        }

        const long double log_tail_ratio = normalized_log_binomial(remaining);

        return log_tail_ratio <= log_target_ratio;
    };

    /*
     * Approximate inverse:
     *
     *     C(x, k)
     *       = x(x-1)...(x-k+1) / k!
     *       ~= (x - (k-1)/2)^k / k!.
     *
     * Therefore
     *
     *     x ~= exp((log(target) + log(k!)) / k) + (k-1)/2,
     *
     * and because x = n - (s + 1),
     *
     *     s ~= n - 1 - x.
     */
    /*
     * EstimateMinimumCandidate still expects an absolute log-tail value.
     * Reconstruct it only for the coarse approximation.
     */
    const long double log_target_tail_for_estimate = log_begin_tail + log_target_ratio;

    SInt candidate = EstimateMinimumCandidate(
        local_begin, local_end, n, k, log_target_tail_for_estimate, cache.log_k_factorial, cache.inv_k);

    const SInt candidate_x = n - (candidate + 1);

    /*
     * Evaluate the candidate in normalized coordinates:
     *
     *   log_candidate_ratio
     *     = log(C(candidate_x, k) / C(begin_x, k)).
     */
    const long double log_candidate_ratio =
        candidate_x >= k ? normalized_log_binomial(candidate_x) : -std::numeric_limits<long double>::infinity();

    /*
     * residual
     *   = log(candidate_tail / begin_tail)
     *     - log(target_tail / begin_tail)
     *
     * The candidate satisfies the predicate iff residual <= 0.
     */
    const long double candidate_residual = log_candidate_ratio - log_target_ratio;

    /*
     * Residual for x == k.
     *
     * Since C(k,k) == 1:
     *
     *   log(C(k,k) / begin_tail)
     *     = -log_begin_tail.
     */
    const long double k_tail_residual = -log_begin_tail - log_target_ratio;

    const CorrectedMinimumCandidate corrected = CorrectMinimumCandidateByRecurrence(
        candidate, local_begin, local_end, n, k, candidate_residual, k_tail_residual, steps);

    candidate = corrected.candidate;

    const bool candidate_valid = corrected.residual <= 0.0L;

    const bool predecessor_invalid = !corrected.has_predecessor || corrected.predecessor_residual > 0.0L;

    if (candidate_valid && predecessor_invalid) {
        return candidate;
    }
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

enum class HugeHypergeometricMode : std::uint8_t {
    ExactCorrectedBinomial,
    ExactSequential,
    BinomialApproximation,
};

inline void ValidateHypergeometricArguments(const CountInt& population, const CountInt& successes, const SInt draws) {
    if (population < 0) {
        throw ConfigurationError("Hypergeometric population must be nonnegative");
    }

    if (successes < 0 || successes > population) {
        throw ConfigurationError("Hypergeometric successes must be in [0, population]");
    }

    if (CountInt(draws) > population) {
        throw ConfigurationError("Hypergeometric draws exceed population");
    }
}

inline std::pair<SInt, SInt>
HypergeometricSupport(const CountInt& population, const CountInt& successes, const SInt draws) {
    const CountInt failures = population - successes;

    const CountInt lower = std::max(CountInt{0}, CountInt(draws) - failures);

    const CountInt upper = std::min(CountInt(draws), successes);

    return {
        lower.convert_to<SInt>(),
        upper.convert_to<SInt>(),
    };
}

/*
 * Deterministic pseudorandom 64-bit word.
 *
 * The rank-tree node seed is identical on every PE visiting the same node,
 * so each PE reconstructs the same split.
 */
inline std::uint64_t HypergeometricRandomWord(SInt& seed) {
    const auto value = sampling::Spooky::hash(static_cast<unsigned long long>(seed++));

    return static_cast<std::uint64_t>(value);
}

/*
 * Exact Bernoulli(numerator / denominator).
 *
 * No floating-point probability is formed.
 *
 * Write:
 *
 *   numerator * 2^64 = quotient * denominator + remainder.
 *
 * The first 64 random bits resolve the decision except for the one interval
 * represented by "remainder". In that case, refine recursively.
 */
inline bool BernoulliExactRationalWithFirstWord(
    CountInt numerator, const CountInt& denominator, SInt& seed, std::uint64_t uniform) {
    if (denominator <= 0) {
        throw ConfigurationError("BernoulliExactRational requires denominator > 0");
    }

    if (numerator <= 0) {
        return false;
    }

    if (numerator >= denominator) {
        return true;
    }

    while (true) {
        const CountInt scaled = numerator << 64;

        const CountInt quotient_big = scaled / denominator;
        const CountInt remainder    = scaled % denominator;

        /*
         * numerator < denominator implies quotient < 2^64.
         */
        const std::uint64_t quotient = quotient_big.convert_to<std::uint64_t>();

        if (uniform < quotient) {
            return true;
        }

        if (uniform > quotient) {
            return false;
        }

        /*
         * uniform == quotient: resolve the residual probability exactly.
         */
        if (remainder == 0) {
            return false;
        }

        numerator = remainder;
        uniform   = HypergeometricRandomWord(seed);
    }
}

inline bool BernoulliExactRational(CountInt numerator, const CountInt& denominator, SInt& seed) {
    const std::uint64_t first_word = HypergeometricRandomWord(seed);
    return BernoulliExactRationalWithFirstWord(std::move(numerator), denominator, seed, first_word);
}

/*
 * Path 2: arbitrary-population exact hypergeometric sampler.
 *
 * This preserves the exact without-replacement distribution. Its runtime is
 * O(draws), so it is the exact fallback, not the scalable path.
 */
inline SInt GenerateHugeHypergeometricExact(CountInt population, CountInt successes, const SInt draws, SInt seed) {
    ValidateHypergeometricArguments(population, successes, draws);

    if (draws == 0 || successes == 0) {
        return 0;
    }

    if (successes == population) {
        return draws;
    }

    SInt sampled_successes = 0;

    for (SInt draw = 0; draw < draws; ++draw) {
        const SInt remaining_draws = draws - draw;

        if (successes == 0) {
            break;
        }

        if (successes == population) {
            sampled_successes += remaining_draws;
            break;
        }

        const CountInt failures = population - successes;

        /*
         * All remaining failures must be selected. Therefore every remaining
         * success must also be selected.
         */
        if (failures == CountInt(remaining_draws)) {
            sampled_successes += successes.convert_to<SInt>();
            break;
        }

        if (BernoulliExactRational(successes, population, seed)) {
            ++sampled_successes;
            --successes;
        }

        --population;
    }

    return sampled_successes;
}

/*
 * Approximate numerator / denominator for the explicitly approximate path.
 *
 * Only leading bits are converted, avoiding overflow when CountInt contains
 * an enormous combinatorial population.
 */
inline long double ApproximateCountIntRatio(const CountInt& numerator, const CountInt& denominator) {
    if (denominator <= 0) {
        throw ConfigurationError("ApproximateCountIntRatio requires denominator > 0");
    }

    if (numerator <= 0) {
        return 0.0L;
    }

    if (numerator >= denominator) {
        return 1.0L;
    }

    constexpr unsigned kRetainedBits = 63;

    const unsigned numerator_bits = boost::multiprecision::msb(numerator) + 1;

    const unsigned denominator_bits = boost::multiprecision::msb(denominator) + 1;

    const unsigned numerator_shift = numerator_bits > kRetainedBits ? numerator_bits - kRetainedBits : 0;

    const unsigned denominator_shift = denominator_bits > kRetainedBits ? denominator_bits - kRetainedBits : 0;

    const std::uint64_t numerator_top = (numerator >> numerator_shift).convert_to<std::uint64_t>();

    const std::uint64_t denominator_top = (denominator >> denominator_shift).convert_to<std::uint64_t>();

    const long double mantissa = static_cast<long double>(numerator_top) / static_cast<long double>(denominator_top);

    const long long exponent = static_cast<long long>(numerator_shift) - static_cast<long long>(denominator_shift);

    if (exponent < std::numeric_limits<long double>::min_exponent) {
        return 0.0L;
    }

    /*
     * numerator < denominator means the real ratio is below one, regardless
     * of intermediate floating-point rounding.
     */
    const long double ratio = std::ldexp(mantissa, static_cast<int>(exponent));

    return std::clamp(ratio, 0.0L, 1.0L);
}

struct ExactDyadicProbability {
    CountInt numerator;
    CountInt denominator;
};

/*
 * Recover the exact binary rational represented by a double in [0, 1].
 *
 * For a binary64 value q, this returns integers Q and D such that:
 *
 *     q = Q / D
 *
 * exactly.
 */
inline ExactDyadicProbability ExactProbabilityFromDouble(const double probability) {
    if (!(probability >= 0.0 && probability <= 1.0)) {
        throw ConfigurationError("ExactProbabilityFromDouble requires probability in [0, 1]");
    }

    if (probability == 0.0) {
        return {
            .numerator   = 0,
            .denominator = 1,
        };
    }

    if (probability == 1.0) {
        return {
            .numerator   = 1,
            .denominator = 1,
        };
    }

    int exponent = 0;

    /*
     * probability = fraction * 2^exponent,
     * with fraction in [0.5, 1).
     *
     * Binary64 has 53 significand bits, including the implicit bit.
     */
    const double fraction = std::frexp(probability, &exponent);

    constexpr int kSignificandBits = 53;

    const std::uint64_t significand = static_cast<std::uint64_t>(std::ldexp(fraction, kSignificandBits));

    CountInt numerator   = significand;
    CountInt denominator = 1;

    const int binary_exponent = exponent - kSignificandBits;

    if (binary_exponent >= 0) {
        numerator <<= binary_exponent;
    } else {
        denominator <<= -binary_exponent;
    }

    /*
     * Reduce powers of two. Both values are nonnegative and the denominator
     * is a power of two.
     */
    while ((numerator & 1) == 0 && denominator > 1) {
        numerator >>= 1;
        denominator >>= 1;
    }

    return {
        .numerator   = std::move(numerator),
        .denominator = std::move(denominator),
    };
}

inline CountInt CountIntGCD(CountInt lhs, CountInt rhs) {
    if (lhs < 0) {
        lhs = -lhs;
    }

    if (rhs < 0) {
        rhs = -rhs;
    }

    while (rhs != 0) {
        CountInt remainder = lhs % rhs;
        lhs                = std::move(rhs);
        rhs                = std::move(remainder);
    }

    return lhs;
}

/*
 * Multiply fraction numerator/denominator by factor_numerator/factor_denominator
 * while cancelling cross factors first.
 */
inline void MultiplyReducedFraction(
    CountInt& numerator, CountInt& denominator, CountInt factor_numerator, CountInt factor_denominator) {
    if (factor_denominator <= 0 || factor_numerator < 0) {
        throw ConfigurationError("Invalid factor in corrected hypergeometric ratio");
    }

    if (factor_numerator == 0) {
        numerator   = 0;
        denominator = 1;
        return;
    }

    CountInt divisor = CountIntGCD(factor_numerator, factor_denominator);

    factor_numerator /= divisor;
    factor_denominator /= divisor;

    /*
     * Cancel the new numerator against the accumulated denominator.
     */
    divisor = CountIntGCD(factor_numerator, denominator);

    factor_numerator /= divisor;
    denominator /= divisor;

    /*
     * Cancel the new denominator against the accumulated numerator.
     */
    divisor = CountIntGCD(factor_denominator, numerator);

    factor_denominator /= divisor;
    numerator /= divisor;

    numerator *= factor_numerator;
    denominator *= factor_denominator;
}

/*
 * Path 3: huge-population binomial approximation.
 *
 * Hypergeometric(N, K, draws) is approximated by
 *
 *   Binomial(draws, K / N).
 *
 * The recursive tree still conserves the exact total because the right child
 * receives draws - left_draws.
 */
template <typename RNG>
SInt GenerateHugeHypergeometricApprox(
    const CountInt& population, const CountInt& successes, const SInt draws, RNG& rng, const SInt seed) {
    ValidateHypergeometricArguments(population, successes, draws);

    if (draws == 0 || successes == 0) {
        return 0;
    }

    if (successes == population) {
        return draws;
    }

    const CountInt failures = population - successes;

    SInt result;

    /*
     * Sample the smaller category to avoid a probability close to one.
     */
    if (successes <= failures) {
        const double probability = static_cast<double>(ApproximateCountIntRatio(successes, population));

        result = rng.GenerateBinomial(seed, draws, probability);
    } else {
        const double failure_probability = static_cast<double>(ApproximateCountIntRatio(failures, population));

        const SInt sampled_failures = rng.GenerateBinomial(seed, draws, failure_probability);

        result = draws - sampled_failures;
    }

    /*
     * A binomial variate can theoretically lie outside the finite-population
     * hypergeometric support. This matters only outside the intended sparse
     * huge-population regime.
     */
    const auto [minimum, maximum] = HypergeometricSupport(population, successes, draws);

    return std::clamp(result, minimum, maximum);
}

inline bool CorrectedRatioIncreasesAt(
    const CountInt& population, const CountInt& successes, const SInt draws, const SInt value,
    const ExactDyadicProbability& proposal_probability) {
    const CountInt failures = population - successes;

    const CountInt proposal_success = proposal_probability.numerator;

    const CountInt proposal_failure = proposal_probability.denominator - proposal_probability.numerator;

    /*
     * R(value + 1) / R(value) >= 1 iff:
     *
     *   (K - value)(1 - q)
     *       >=
     *   (N - K - m + value + 1)q.
     */
    const CountInt lhs = (successes - value) * proposal_failure;

    const CountInt rhs = (failures - draws + value + 1) * proposal_success;

    return lhs >= rhs;
}

inline SInt CorrectedHypergeometricRatioMode(
    const CountInt& population, const CountInt& successes, const SInt draws,
    const ExactDyadicProbability& proposal_probability) {
    const auto [support_begin, support_end] = HypergeometricSupport(population, successes, draws);

    if (support_begin >= support_end) {
        return support_begin;
    }

    /*
     * Find the last value x for which R(x + 1) >= R(x).
     */
    SInt left  = support_begin;
    SInt right = support_end - 1;

    if (!CorrectedRatioIncreasesAt(population, successes, draws, left, proposal_probability)) {
        return support_begin;
    }

    while (left < right) {
        const SInt middle = left + ((right - left + 1) / 2);

        if (CorrectedRatioIncreasesAt(population, successes, draws, middle, proposal_probability)) {
            left = middle;
        } else {
            right = middle - 1;
        }
    }

    return std::min<SInt>(left + 1, support_end);
}

/*
 * Three-path dispatcher.
 *
 * Path 1:
 *   population fits SInt -> existing exact KaGen hypergeometric sampler.
 *
 * Path 2:
 *   huge population + Exact -> arbitrary-precision exact sequential sampler.
 *
 * Path 3:
 *   huge population + BinomialApproximation -> fast binomial split.
 */
template <typename RNG>
SInt GenerateHypergeometricHybrid(
    const CountInt& population, const CountInt& successes, const SInt draws, RNG& rng, const SInt seed,
    const HugeHypergeometricMode huge_mode) {
    ValidateHypergeometricArguments(population, successes, draws);

    if (draws == 0 || successes == 0) {
        return 0;
    }

    if (successes == population) {
        return draws;
    }

    const CountInt sint_max = std::numeric_limits<SInt>::max();

    /*
     * Path 1: machine-sized exact sampler.
     *
     * Since successes <= population, checking population is sufficient.
     */
    if (population <= sint_max) {
        const SInt result =
            rng.GenerateHypergeometric(seed, successes.convert_to<SInt>(), draws, population.convert_to<SInt>());

        const auto [minimum, maximum] = HypergeometricSupport(population, successes, draws);

        if (result < minimum || result > maximum) {
            throw ConfigurationError(
                "Native hypergeometric sampler returned "
                "a value outside its support");
        }

        return result;
    }

    if (huge_mode == HugeHypergeometricMode::BinomialApproximation) {
        return GenerateHugeHypergeometricApprox(population, successes, draws, rng, seed);
    }

    if (huge_mode == HugeHypergeometricMode::ExactCorrectedBinomial) {
        return GenerateHugeHypergeometricCorrectedBinomial(population, successes, draws, rng, seed);
    }

    /*
     * Explicit slow reference/fallback path.
     */
    const SInt result = GenerateHugeHypergeometricExact(population, successes, draws, seed);

    const auto [minimum, maximum] = HypergeometricSupport(population, successes, draws);

    if (result < minimum || result > maximum) {
        throw ConfigurationError(
            "Huge exact hypergeometric sampler returned "
            "a value outside its support");
    }

    return result;
}

/*
 * Fast prefilter for the corrected-binomial rejection step.
 *
 * It evaluates the acceptance ratio in log space and compares it with the
 * interval represented by the first 64 random bits. Decisions too close to
 * the boundary are deliberately sent to the original exact rational path.
 */
enum class FastAcceptanceDecision : std::uint8_t {
    Accept,
    Reject,
    NeedExactFallback,
};

inline long double LogPositiveCountIntApprox(const CountInt& value) {
    if (value <= 0) {
        return -std::numeric_limits<long double>::infinity();
    }

    constexpr unsigned kRetainedBits = 64;

    const unsigned bit_count = boost::multiprecision::msb(value) + 1;
    const unsigned shift     = bit_count > kRetainedBits ? bit_count - kRetainedBits : 0;

    const std::uint64_t top = (value >> shift).convert_to<std::uint64_t>();

    constexpr long double kLogTwo = 0.693147180559945309417232121458176568L;
    return std::log(static_cast<long double>(top)) + static_cast<long double>(shift) * kLogTwo;
}

inline long double CorrectedHypergeometricLogAcceptanceRatio(
    const CountInt& population, const CountInt& successes, const SInt draws, const SInt value, const SInt mode,
    const double proposal_probability) {
    const CountInt failures = population - successes;

    const long double q               = static_cast<long double>(proposal_probability);
    const long double log_q           = std::log(q);
    const long double log_one_minus_q = std::log1pl(-q);
    long double       log_ratio       = 0.0L;

    if (value > mode) {
        for (SInt j = mode; j < value; ++j) {
            log_ratio += LogPositiveCountIntApprox(successes - j);
            log_ratio += log_one_minus_q;
            log_ratio -= LogPositiveCountIntApprox(failures - draws + j + 1);
            log_ratio -= log_q;
        }
    } else if (value < mode) {
        for (SInt j = value; j < mode; ++j) {
            log_ratio += LogPositiveCountIntApprox(failures - draws + j + 1);
            log_ratio += log_q;
            log_ratio -= LogPositiveCountIntApprox(successes - j);
            log_ratio -= log_one_minus_q;
        }
    }

    // The exact value is at most zero. Suppress tiny positive roundoff.
    return std::min(0.0L, log_ratio);
}

inline FastAcceptanceDecision FilterCorrectedHypergeometricAcceptance(
    const long double log_acceptance, const std::uint64_t random_word, const SInt ratio_distance) {
    constexpr long double kLogTwo     = 0.693147180559945309417232121458176568L;
    constexpr long double kLogTwoTo64 = 64.0L * kLogTwo;

    const long double lower = random_word == 0 ? -std::numeric_limits<long double>::infinity()
                                               : std::log(static_cast<long double>(random_word)) - kLogTwoTo64;

    const long double upper = random_word == std::numeric_limits<std::uint64_t>::max()
                                  ? 0.0L
                                  : std::log(static_cast<long double>(random_word) + 1.0L) - kLogTwoTo64;

    /*
     * Conservative engineering guard. Any close comparison falls back to
     * exact cpp_int arithmetic. Increase this multiplier if validation on a
     * platform with a different long-double implementation finds ambiguity.
     */
    const long double finite_lower = std::isfinite(lower) ? std::abs(lower) : 0.0L;
    const long double scale        = 1.0L + std::abs(log_acceptance) + finite_lower + std::abs(upper);
    const long double operations   = 64.0L + 32.0L * static_cast<long double>(ratio_distance);
    const long double guard        = operations * std::numeric_limits<long double>::epsilon() * scale;

    if (upper + guard < log_acceptance) {
        return FastAcceptanceDecision::Accept;
    }

    if (lower - guard > log_acceptance) {
        return FastAcceptanceDecision::Reject;
    }

    return FastAcceptanceDecision::NeedExactFallback;
}

template <typename RNG>
SInt GenerateBinomialHybrid(
    const CountInt& population, const long double log_population, const long double log_probability, RNG& rng,
    const SInt seed, const char* error_context) {
    if (population < 0) {
        throw ConfigurationError(std::string(error_context) + " binomial population must be nonnegative");
    }

    if (population == 0) {
        return 0;
    }

    if (log_probability == -std::numeric_limits<long double>::infinity()) {
        return 0;
    }

    if (!(log_probability <= 0.0L)) {
        throw ConfigurationError(std::string(error_context) + " binomial probability exceeds one");
    }

    const CountInt native_limit = std::numeric_limits<SInt>::max();

    if (population <= native_limit) {
        const SInt trials = population.convert_to<SInt>();

        const double probability = static_cast<double>(std::exp(log_probability));

        if (probability <= 0.0) {
            return 0;
        }

        return rng.GenerateBinomial(seed, trials, probability);
    }

    return GenerateHugeBinomialCorrectedPoisson(population, log_population, log_probability, rng, seed, error_context);
}

template <typename RNG>
SInt GenerateBinomialHybrid(
    const CountInt& population, const long double log_probability, RNG& rng, const SInt seed,
    const char* error_context) {
    if (population <= CountInt(std::numeric_limits<SInt>::max())) {
        return GenerateBinomialHybrid(
            population,
            0.0L, // unused by the native branch
            log_probability, rng, seed, error_context);
    }

    return GenerateBinomialHybrid(
        population, LogPositiveCountIntApprox(population), log_probability, rng, seed, error_context);
}

inline long double TransformCountIntToLongDoubleFast(const CountInt& value) {
    if (value <= 0) {
        return 0.0L;
    }

    constexpr unsigned kRetainedBits = std::numeric_limits<long double>::digits;

    const unsigned bit_count = boost::multiprecision::msb(value) + 1;

    const unsigned shift = bit_count > kRetainedBits ? bit_count - kRetainedBits : 0;

    const std::uint64_t top = (value >> shift).convert_to<std::uint64_t>();

    return std::ldexp(static_cast<long double>(top), static_cast<int>(shift));
}

template <typename RNG>
SInt GenerateHugeBinomialCorrectedPoisson(
    const CountInt& population, const long double log_population, const long double log_probability, RNG& rng,
    const SInt seed, const char* error_context) {
    if (population < 0) {
        throw ConfigurationError(std::string(error_context) + " binomial population must be nonnegative");
    }

    if (population == 0) {
        return 0;
    }

    if (!std::isfinite(log_probability)) {
        if (log_probability == -std::numeric_limits<long double>::infinity()) {
            return 0;
        }

        throw ConfigurationError(std::string(error_context) + " invalid log probability");
    }

    if (log_probability >= 0.0L) {
        // p >= 1. CIGAM validation should normally prevent p > 1.
        if (log_probability > 0.0L) {
            throw ConfigurationError(std::string(error_context) + " binomial probability exceeds one");
        }

        if (population > CountInt(std::numeric_limits<SInt>::max())) {
            throw ConfigurationError(std::string(error_context) + " binomial result exceeds SInt");
        }

        return population.convert_to<SInt>();
    }

    /*
     * Compute
     *
     *     log(mean) = log(N) + log(p)
     *
     * without converting the complete population directly to a native
     * integer.
     */
    const long double log_mean = log_population + log_probability;

    const long double log_denorm_min = std::log(static_cast<long double>(std::numeric_limits<double>::denorm_min()));

    if (log_mean <= log_denorm_min) {
        return 0;
    }

    const long double log_sint_max = std::log(static_cast<long double>(std::numeric_limits<SInt>::max()));

    if (log_mean > log_sint_max) {
        throw ConfigurationError(std::string(error_context) + " expected block count exceeds SInt");
    }

    const double proposal_mean = std::exp(static_cast<double>(log_mean));

    if (!std::isfinite(proposal_mean) || proposal_mean <= 0.0) {
        return 0;
    }

    /*
     * The correction targets Binomial(N, p_eff), where
     *
     *     p_eff = proposal_mean / N.
     *
     * This ensures that the Poisson proposal mean and the corrected
     * binomial mean agree despite the proposal_mean being rounded to
     * binary64 for GeneratePoisson().
     */
    const long double population_ld = TransformCountIntToLongDoubleFast(population);

    if (!std::isfinite(population_ld) || population_ld <= 0.0L) {
        throw ConfigurationError(std::string(error_context) + " population cannot be represented as long double");
    }

    const long double effective_probability = static_cast<long double>(proposal_mean) / population_ld;

    if (effective_probability <= 0.0L || effective_probability >= 1.0L) {
        /*
         * This path is intended for huge N and small p. A native-sized or
         * dense population should use the native binomial sampler.
         */
        throw ConfigurationError(std::string(error_context) + " corrected-Poisson path requires 0 < p_eff < 1");
    }

    /*
     * For lambda = N p_eff,
     *
     *     R(k + 1) / R(k)
     *       = (N - k) / (N (1 - p_eff)),
     *
     * where R(k) is BinomialPMF(k) / PoissonPMF(k).
     *
     * A maximizing value is ceil(lambda). If lambda is integral, the two
     * adjacent values are both maxima and either is suitable.
     */
    const SInt        ratio_mode                = static_cast<SInt>(std::ceil(static_cast<long double>(proposal_mean)));
    const long double log_one_minus_probability = std::log1pl(-effective_probability);
    auto              log_acceptance_ratio      = [&](const SInt value) -> long double {
        long double log_acceptance = 0.0L;

        if (value > ratio_mode) {
            for (SInt j = ratio_mode; j < value; ++j) {
                /*
                 * log((N-j) / N) - log(1-p_eff)
                 *
                 * Use log1p so that tiny j/N and p_eff are retained.
                 */
                const long double j_over_population = static_cast<long double>(j) / population_ld;

                log_acceptance += std::log1pl(-j_over_population) - log_one_minus_probability;
            }
        } else if (value < ratio_mode) {
            for (SInt j = value; j < ratio_mode; ++j) {
                const long double j_over_population = static_cast<long double>(j) / population_ld;

                /*
                 * Reciprocal of the rightward ratio.
                 */
                log_acceptance += log_one_minus_probability - std::log1pl(-j_over_population);
            }
        }

        /*
         * The exact value is at most zero. Tiny positive values can occur
         * through floating-point roundoff.
         */
        return std::min(0.0L, log_acceptance);
    };

    constexpr SInt kMaximumAttempts = 1'000'000;

    for (SInt attempt = 0; attempt < kMaximumAttempts; ++attempt) {
        const SInt proposal_seed = sampling::Spooky::hash(
            static_cast<unsigned long long>(seed) + (static_cast<unsigned long long>(attempt) * 0x9e3779b97f4a7c15ULL));

        const SInt proposal = rng.GeneratePoisson(proposal_seed, proposal_mean);

        const long double log_acceptance = log_acceptance_ratio(proposal);

        const SInt acceptance_seed = sampling::Spooky::hash(
            static_cast<unsigned long long>(seed) ^ 0xd1b54a32d192ed03ULL
            ^ (static_cast<unsigned long long>(attempt) * 0x94d049bb133111ebULL));

        const long double uniform = std::min<long double>(
            static_cast<long double>(rng.GenerateUniform(acceptance_seed, 0.0, 1.0)), std::nextafter(1.0L, 0.0L));

        const long double acceptance_probability = std::exp(log_acceptance);

        if (uniform <= acceptance_probability) {
            return proposal;
        }
    }

    throw ConfigurationError(std::string(error_context) + " corrected-Poisson binomial sampler exceeded attempt limit");
}

inline std::pair<CountInt, CountInt> CorrectedHypergeometricAcceptanceRatio(
    const CountInt& population, const CountInt& successes, const SInt draws, const SInt value, const SInt mode,
    const ExactDyadicProbability& proposal_probability) {
    const CountInt failures = population - successes;

    const CountInt proposal_success = proposal_probability.numerator;

    const CountInt proposal_failure = proposal_probability.denominator - proposal_probability.numerator;

    CountInt numerator   = 1;
    CountInt denominator = 1;

    if (value > mode) {
        /*
         * Walk right:
         *
         * R(j + 1) / R(j)
         *   = (K - j)(1-q)
         *     / ((N-K-m+j+1)q).
         */
        for (SInt j = mode; j < value; ++j) {
            const CountInt factor_numerator = (successes - j) * proposal_failure;

            const CountInt factor_denominator = (failures - draws + j + 1) * proposal_success;

            MultiplyReducedFraction(numerator, denominator, factor_numerator, factor_denominator);
        }
    } else if (value < mode) {
        /*
         * Walk left using the reciprocal:
         *
         * R(j) / R(j + 1)
         *   = (N-K-m+j+1)q
         *     / ((K-j)(1-q)).
         */
        for (SInt j = value; j < mode; ++j) {
            const CountInt factor_numerator = (failures - draws + j + 1) * proposal_success;

            const CountInt factor_denominator = (successes - j) * proposal_failure;

            MultiplyReducedFraction(numerator, denominator, factor_numerator, factor_denominator);
        }
    }

    /*
     * The mode maximizes R, so this ratio must be in [0, 1].
     */
    if (numerator < 0 || denominator <= 0 || numerator > denominator) {
        throw ConfigurationError("Invalid corrected-binomial acceptance ratio");
    }

    return {
        std::move(numerator),
        std::move(denominator),
    };
}

template <typename RNG>
SInt GenerateHugeHypergeometricCorrectedBinomial(
    const CountInt& population, const CountInt& successes, const SInt draws, RNG& rng, const SInt seed) {
    ValidateHypergeometricArguments(population, successes, draws);

    if (draws == 0 || successes == 0) {
        return 0;
    }

    if (successes == population) {
        return draws;
    }

    const auto [support_begin, support_end] = HypergeometricSupport(population, successes, draws);

    /*
     * GenerateBinomial receives a double, so this double defines the actual
     * proposal distribution. The rejection correction below uses its exact
     * binary value.
     */
    const double proposal_probability = static_cast<double>(ApproximateCountIntRatio(successes, population));

    if (proposal_probability <= 0.0) {
        /*
         * The rounded proposal cannot represent the success probability.
         * Retain the sequential exact fallback.
         */
        return GenerateHugeHypergeometricExact(population, successes, draws, seed);
    }

    if (proposal_probability >= 1.0) {
        return GenerateHugeHypergeometricExact(population, successes, draws, seed);
    }

    const ExactDyadicProbability exact_proposal_probability = ExactProbabilityFromDouble(proposal_probability);

    const SInt ratio_mode = CorrectedHypergeometricRatioMode(population, successes, draws, exact_proposal_probability);

    constexpr SInt kMaximumAttempts = 1'000'000;

    for (SInt attempt = 0; attempt < kMaximumAttempts; ++attempt) {
        const SInt proposal_seed = sampling::Spooky::hash(
            static_cast<unsigned long long>(seed) + static_cast<unsigned long long>(attempt) * 0x9e3779b97f4a7c15ULL);

        const SInt proposal = rng.GenerateBinomial(proposal_seed, draws, proposal_probability);

        /*
         * Rejection, not clamping. Clamping would change the distribution.
         */
        if (proposal < support_begin || proposal > support_end) {
            continue;
        }

        /*
         * Use a separate deterministic random stream for acceptance. The
         * fast filter and exact fallback inspect the same first random word.
         */
        SInt acceptance_seed = sampling::Spooky::hash(
            static_cast<unsigned long long>(seed) ^ 0xd1b54a32d192ed03ULL
            ^ (static_cast<unsigned long long>(attempt) * 0x94d049bb133111ebULL));

        const std::uint64_t first_acceptance_word = HypergeometricRandomWord(acceptance_seed);

        const long double log_acceptance = CorrectedHypergeometricLogAcceptanceRatio(
            population, successes, draws, proposal, ratio_mode, proposal_probability);

        const SInt ratio_distance = proposal >= ratio_mode ? proposal - ratio_mode : ratio_mode - proposal;

        const FastAcceptanceDecision fast_decision =
            FilterCorrectedHypergeometricAcceptance(log_acceptance, first_acceptance_word, ratio_distance);

        if (fast_decision == FastAcceptanceDecision::Accept) {
            return proposal;
        }

        if (fast_decision == FastAcceptanceDecision::Reject) {
            continue;
        }

        /*
         * Ambiguous boundary case: retain the original exact rational
         * computation and reuse the already-consumed first random word.
         */
        const auto [acceptance_numerator, acceptance_denominator] = CorrectedHypergeometricAcceptanceRatio(
            population, successes, draws, proposal, ratio_mode, exact_proposal_probability);

        if (BernoulliExactRationalWithFirstWord(
                acceptance_numerator, acceptance_denominator, acceptance_seed, first_acceptance_word)) {
            return proposal;
        }
    }

    /*
     * This should be practically unreachable for the sparse regime.
     */
    return GenerateHugeHypergeometricExact(population, successes, draws, seed);
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

struct HypergraphMemoryStats {
    SInt max_pins_per_hyperedge = 0;
};
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
struct HGNPInstrumentation {
    // Whole-generator phases.
    double input_seconds      = 0.0;
    double budget_seconds     = 0.0;
    double planning_seconds   = 0.0;
    double reserve_seconds    = 0.0;
    double generation_seconds = 0.0;

    // Planning subphases.
    double boundary_seconds       = 0.0;
    double population_seconds     = 0.0;
    double count_sampling_seconds = 0.0;

    // Generation subphases.
    double cache_init_seconds      = 0.0;
    double minimum_sample_seconds  = 0.0;
    double pin_sample_seconds      = 0.0;
    double duplicate_check_seconds = 0.0;
    double csr_write_seconds       = 0.0;

    // Work counters.
    std::uint64_t active_sizes     = 0;
    std::uint64_t zero_count_sizes = 0;
    std::uint64_t planned_edges    = 0;
    std::uint64_t planned_pins     = 0;

    std::uint64_t generated_edges   = 0;
    std::uint64_t generated_pins    = 0;
    std::uint64_t attempts          = 0;
    std::uint64_t duplicate_rejects = 0;

    // Minimum sampler and cache behavior.
    std::uint64_t minimum_samples      = 0;
    std::uint64_t minimum_search_steps = 0;
    std::uint64_t minimum_cache_gets   = 0;

    // Map-backed LogBinomCache::Get().
    std::uint64_t cache_map_calls   = 0;
    std::uint64_t cache_map_hits    = 0;
    std::uint64_t cache_map_misses  = 0;
    std::uint64_t cache_map_inserts = 0;

    // Cursor-backed LogBinomCache::GetCandidate().
    std::uint64_t cache_candidate_calls                 = 0;
    std::uint64_t cache_candidate_exact_hits            = 0;
    std::uint64_t cache_candidate_recurrence_hits       = 0;
    std::uint64_t cache_candidate_direct_evals          = 0;
    std::uint64_t cache_candidate_below_range           = 0;
    std::uint64_t cache_candidate_distance_sum          = 0;
    std::uint64_t cache_candidate_max_distance          = 0;
    std::uint64_t cache_candidate_distance_observations = 0;

    std::uint64_t cache_backward_steps = 0;
    std::uint64_t cache_forward_steps  = 0;

    std::uint64_t cache_max_size = 0;

    SInt cache_min_key = std::numeric_limits<SInt>::max();
    SInt cache_max_key = 0;

    // Important rank-local maxima.
    std::uint64_t max_local_edges_for_size  = 0;
    std::uint64_t max_cache_initial_reserve = 0;
};
class ScopedMPITimer {
public:
    explicit ScopedMPITimer(double& accumulator) : accumulator_(accumulator), start_(MPI_Wtime()) {}

    ScopedMPITimer(const ScopedMPITimer&)            = delete;
    ScopedMPITimer& operator=(const ScopedMPITimer&) = delete;

    ~ScopedMPITimer() {
        accumulator_ += MPI_Wtime() - start_;
    }

private:
    double& accumulator_;
    double  start_;
};
#endif

struct FixedCountLocalRange {
    SInt min_begin = 0;
    SInt min_end   = 0;
    SInt local_m   = 0;
};

// Result of the zero-communication minimum-partition count allocation.
// Each rank follows only its own root-to-leaf path in the deterministic
// binary split tree. No global rank-count vector is materialized.
struct MinimumPartitionLocalCount {
    SInt     local_m          = 0;
    CountInt local_population = 0;
};

template <typename BigInt>
class ExactFixedCountHyperedgeGenerator {
public:
    ExactFixedCountHyperedgeGenerator(
        const PGeneratorConfig& config, SInt partition_id, SInt num_partitions, RNGWrapper<>& rng, Mersenne& mersenne,
        Graph& graph, HypergraphMemoryStats& memory_stats, FloydScratchSet& floyd_scratch,
        ErdosHypergraphDebugLogger* debug_logger, SInt* next_debug_hyperedge_id
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        ,
        HGNPInstrumentation* instrumentation
#endif
    );

    FixedCountLocalRange PrepareLocalRange(SInt hyperedge_size, SInt global_m);

    void Generate(SInt hyperedge_size, const FixedCountLocalRange& range);

private:
    MinimumPartitionLocalCount ComputeLocalCount(SInt hyperedge_size, SInt global_m);

    MinimumPartitionLocalCount ComputeLocalCountRecursive(
        SInt hyperedge_size, SInt rank_begin, SInt rank_end, SInt target_rank, SInt min_begin, SInt min_end,
        const CountInt& population, SInt draws, SInt level);

    double SplitProbability(const CountInt& child_population, const CountInt& parent_population) const;

    void ValidateDuplicateCheckingIsFeasible(SInt local_m) const;

    void ValidateExactSparseDensity(SInt local_m, const CountInt& local_population) const;

    SInt EdgeSeed(SInt hyperedge_size) const;

    SInt RankSplitSeed(SInt hyperedge_size, SInt rank_begin, SInt rank_end, SInt level) const;

    void GenerateSampledHyperedges(
        SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, SInt local_m, SInt& edge_seed,
        LogBinomCache& log_binom_cache);

    void SampleHyperedgeInto(SInt minimum_vertex, SInt hyperedge_size, std::vector<SInt>& pins);

    bool TryPushHyperedge(const std::vector<SInt>& pins, HyperedgeSeenSet& local_seen);

    const PGeneratorConfig& config_;
    SInt                    partition_id_;
    SInt                    num_partitions_;

    RNGWrapper<>&               rng_;
    Mersenne&                   mersenne_;
    Graph&                      graph_;
    HypergraphMemoryStats&      memory_stats_;
    FloydScratchSet&            floyd_scratch_;
    ErdosHypergraphDebugLogger* debug_logger_            = nullptr;
    SInt*                       next_debug_hyperedge_id_ = nullptr;
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    HGNPInstrumentation* instrumentation_;
#endif
    std::size_t log_binom_cache_reserve_hint_ = 4096;
};

inline void ObservePinsAndMaybeReserve(Graph& graph, HypergraphMemoryStats& stats, const std::size_t pins) {
    stats.max_pins_per_hyperedge = std::max(stats.max_pins_per_hyperedge, static_cast<SInt>(pins));

    const std::size_t needed_pins = graph.hyperedge_pins.size() + pins;
    if (needed_pins > graph.hyperedge_pins.capacity()) {
        const std::size_t slack = std::min<std::size_t>(needed_pins / 8, 1'000'000);
        graph.hyperedge_pins.reserve(needed_pins + slack);
    }
}

inline void PushUncompressedHyperedge(Graph& graph, HypergraphMemoryStats& stats, const std::vector<SInt>& pins) {
    ObservePinsAndMaybeReserve(graph, stats, pins.size());

    if (graph.hyperedge_offsets.empty()) {
        graph.hyperedge_offsets.push_back(0);
    }

    graph.hyperedge_pins.insert(graph.hyperedge_pins.end(), pins.begin(), pins.end());
    graph.hyperedge_offsets.push_back(static_cast<SInt>(graph.hyperedge_pins.size()));
}

inline std::pair<SInt, SInt> BalancedVertexRange(const SInt n, const PEID rank, const PEID size) {
    const SInt vertices_per_pe = n / size;
    const SInt remainder       = n % size;

    const SInt start = (rank * vertices_per_pe) + std::min<SInt>(rank, remainder);
    const SInt end   = start + vertices_per_pe + static_cast<SInt>(static_cast<SInt>(rank) < remainder);

    return {start, end};
}

inline HyperedgeSeenSet MakeLocalSeenSet(const bool allow_duplicates, const SInt local_edge_count) {
    HyperedgeSeenSet local_seen;

    if (!allow_duplicates) {
        local_seen.max_load_factor(0.5f);
        local_seen.reserve(static_cast<std::size_t>(local_edge_count));
    }

    return local_seen;
}

class UncompressedHyperedgeBuilder {
public:
    void Begin();
    void AddPin(SInt pin);
    void Commit(Graph& graph, HypergraphMemoryStats& stats);
    void Abort();

private:
    std::vector<SInt> pins_;
};

} // namespace kagen