#include "kagen/generators/hyper/h_erdos/hyper_gnp.h"

#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"

#include <sched.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace kagen {

std::unique_ptr<Generator>
HyperGNPFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    if (config.n <= (1ull << 31)) {
        return std::make_unique<HyperGNPSmall>(config, rank, size);
    }
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

    if (config_.debug) {
        debug_logger_.emplace(
            config_.output_graph.filename + "_hgnp_debug_pe_" + std::to_string(rank_) + ".csv", false);
    }
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateCSR() {
    const SInt lower_bound = config_.size_dist_lower_bound;
    const SInt hard_upper_bound =
        config_.size_dist_upper_bound != 0 ? std::min(config_.size_dist_upper_bound, config_.n) : config_.n;

    if (lower_bound > hard_upper_bound) {
        throw ConfigurationError("Invalid HGNP hyperedge size range");
    }

    {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        ScopedMPITimer timer(instrumentation_.input_seconds);
#endif
        LoadSizeDistributionInputs();
        SelectProbabilityMode();
    }

    const SInt upper_bound = EffectiveUpperBound(hard_upper_bound, lower_bound);

    if (probs_type_ == ProbabilityMode::EdgeAndPinBudget) {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        ScopedMPITimer timer(instrumentation_.budget_seconds);
#endif
        GenerateBudgetSizeCounts(lower_bound, hard_upper_bound);
    }

    std::vector<HGNPSizePlan> plan;
    {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        ScopedMPITimer timer(instrumentation_.planning_seconds);
#endif
        plan = BuildGenerationPlan(lower_bound, upper_bound);
    }

    {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        ScopedMPITimer timer(instrumentation_.reserve_seconds);
#endif
        ReserveCSRForPlan(plan);
    }

    graph_.hyperedge_offsets.push_back(0);

    {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        ScopedMPITimer timer(instrumentation_.generation_seconds);
#endif
        GenerateHyperedgesFromPlan(plan);
    }

    SetLocalVertexRange();
}

template <typename BigInt>
void HyperGNP<BigInt>::ReserveCSRForPlan(const std::vector<HGNPSizePlan>& plan) {
    std::size_t expected_edges = 0;
    std::size_t expected_pins  = 0;

    for (const HGNPSizePlan& entry: plan) {
        const std::size_t local_m = static_cast<std::size_t>(entry.range.local_m);

        const std::size_t reserve_m = local_m + (local_m / 8) + 16;

        const std::size_t k = static_cast<std::size_t>(entry.hyperedge_size);

        if (reserve_m > std::numeric_limits<std::size_t>::max() - expected_edges) {
            throw ConfigurationError("HGNP local edge reservation overflows size_t");
        }

        if (k != 0 && reserve_m > (std::numeric_limits<std::size_t>::max() - expected_pins) / k) {
            throw ConfigurationError("HGNP local pin reservation overflows size_t");
        }

        expected_edges += reserve_m;
        expected_pins += reserve_m * k;
    }

    graph_.hyperedge_offsets.reserve(expected_edges + 1);
    graph_.hyperedge_pins.reserve(expected_pins);
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateHyperedgesFromPlan(const std::vector<HGNPSizePlan>& plan) {
    for (const HGNPSizePlan& entry: plan) {
        GenerateHyperedgesFromSizePlan(entry);
    }
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateHyperedgesFromSizePlan(const HGNPSizePlan& entry) {
    if (entry.range.local_m == 0) {
        return;
    }

    if (config_.approx) {
        GenerateApproxHyperedgesFromPlan(entry);
        return;
    }

    ExactFixedCountHyperedgeGenerator<BigInt> fixed_count_generator(
        config_, entry.partition_id, config_.k, rng_, mersenne_, graph_, memory_stats_, floyd_scratch_,
        debug_logger_ ? &*debug_logger_ : nullptr, &next_debug_hyperedge_id_
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        ,
        &instrumentation_
#endif
    );

    fixed_count_generator.Generate(
        entry.hyperedge_size, FixedCountLocalRange{
                                  .min_begin = entry.range.begin,
                                  .min_end   = entry.range.end,
                                  .local_m   = entry.range.local_m,
                              });
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateApproxHyperedgesFromPlan(const HGNPSizePlan& entry) {
    const std::size_t cache_size = ComputeLogBinomCacheSize(entry.range.local_m, entry.range.begin, entry.range.end);

    const double cache_start = MPI_Wtime();

    LogBinomCache cache(entry.hyperedge_size, cache_size);

    cache.InitializeRange(config_.n, entry.range.begin, entry.range.end);
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    instrumentation_.cache_init_seconds += MPI_Wtime() - cache_start;

    instrumentation_.max_cache_initial_reserve =
        std::max<std::uint64_t>(instrumentation_.max_cache_initial_reserve, cache_size);
#endif
    std::vector<SInt> pins;
    pins.reserve(static_cast<std::size_t>(entry.hyperedge_size));

    auto seen = MakeLocalSeenSet(config_.allow_duplicates, entry.range.local_m);

    SInt edge_seed = LocalEdgeSeed(entry.hyperedge_size, entry.partition_id);

    mersenne_.RandomInit(edge_seed);

    GenerateLocalHyperedges(entry, edge_seed, cache, pins, seen);

    const LogBinomCacheStats& stats = cache.stats;
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    instrumentation_.cache_map_calls += stats.map_calls;

    instrumentation_.cache_map_hits += stats.map_hits;

    instrumentation_.cache_map_misses += stats.map_misses;

    instrumentation_.cache_map_inserts += stats.map_inserts;

    instrumentation_.cache_candidate_calls += stats.candidate_calls;

    instrumentation_.cache_candidate_exact_hits += stats.candidate_exact_hits;

    instrumentation_.cache_candidate_recurrence_hits += stats.candidate_recurrence_hits;

    instrumentation_.cache_candidate_direct_evals += stats.candidate_direct_evals;

    instrumentation_.cache_candidate_below_range += stats.candidate_below_range;

    instrumentation_.cache_candidate_distance_sum += stats.candidate_distance_sum;

    instrumentation_.cache_candidate_max_distance =
        std::max(instrumentation_.cache_candidate_max_distance, stats.candidate_max_distance);

    instrumentation_.cache_candidate_distance_observations += stats.candidate_distance_observations;

    instrumentation_.cache_backward_steps += stats.backward_steps;

    instrumentation_.cache_forward_steps += stats.forward_steps;

    instrumentation_.cache_max_size = std::max(instrumentation_.cache_max_size, stats.max_size);

    const std::uint64_t total_calls = stats.map_calls + stats.candidate_calls;

    if (total_calls > 0) {
        instrumentation_.cache_min_key = std::min(instrumentation_.cache_min_key, stats.min_key);

        instrumentation_.cache_max_key = std::max(instrumentation_.cache_max_key, stats.max_key);
    }
#endif
}

template <typename BigInt>
void HyperGNP<BigInt>::FinalizeCSR(MPI_Comm comm) {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    std::cerr << "rank=" << rank_ << " pid=" << getpid() << " cpu=" << sched_getcpu() << '\n';
    struct TimingRow {
        const char* name;
        double HGNPInstrumentation::* member;
    };

    constexpr std::array timing_rows{
        TimingRow{"Input", &HGNPInstrumentation::input_seconds},
        TimingRow{"Budget", &HGNPInstrumentation::budget_seconds},
        TimingRow{"Planning", &HGNPInstrumentation::planning_seconds},
        TimingRow{"  Boundaries", &HGNPInstrumentation::boundary_seconds},
        TimingRow{"  Exact population", &HGNPInstrumentation::population_seconds},
        TimingRow{"  Count sampling", &HGNPInstrumentation::count_sampling_seconds},
        TimingRow{"Reserve", &HGNPInstrumentation::reserve_seconds},
        TimingRow{"Generation", &HGNPInstrumentation::generation_seconds},
        TimingRow{"  Cache init", &HGNPInstrumentation::cache_init_seconds},
        TimingRow{"  Minimum sampling", &HGNPInstrumentation::minimum_sample_seconds},
        TimingRow{"  Pin sampling", &HGNPInstrumentation::pin_sample_seconds},
        TimingRow{"  Duplicate checking", &HGNPInstrumentation::duplicate_check_seconds},
        TimingRow{"  CSR writes", &HGNPInstrumentation::csr_write_seconds},
    };

    if (rank_ == 0) {
        std::cout << "\nHGNP instrumentation:\n"
                  << std::left << std::setw(24) << "Phase" << std::right << std::setw(12) << "Min" << std::setw(12)
                  << "Mean" << std::setw(12) << "Max" << std::setw(12) << "Max/Mean" << '\n';
    }

    for (const auto& row: timing_rows) {
        const double local = instrumentation_.*(row.member);

        double minimum = 0.0;
        double maximum = 0.0;
        double sum     = 0.0;

        MPI_Reduce(&local, &minimum, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
        MPI_Reduce(&local, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&local, &maximum, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

        if (rank_ == 0) {
            const double mean = sum / static_cast<double>(size_);

            const double imbalance = mean > 0.0 ? maximum / mean : 0.0;

            std::cout << std::left << std::setw(24) << row.name << std::right << std::fixed << std::setprecision(6)
                      << std::setw(12) << minimum << std::setw(12) << mean << std::setw(12) << maximum << std::setw(12)
                      << imbalance << '\n';
        }
    }

    auto reduce_sum_u64 = [&](const std::uint64_t local) {
        std::uint64_t global = 0;

        MPI_Reduce(&local, &global, 1, MPI_UINT64_T, MPI_SUM, 0, comm);

        return global;
    };

    auto reduce_max_u64 = [&](const std::uint64_t local) {
        std::uint64_t global = 0;

        MPI_Reduce(&local, &global, 1, MPI_UINT64_T, MPI_MAX, 0, comm);

        return global;
    };

    auto reduce_min_sint = [&](const SInt local) {
        SInt global = 0;

        MPI_Reduce(&local, &global, 1, MPI_UINT64_T, MPI_MIN, 0, comm);

        return global;
    };

    auto reduce_max_sint = [&](const SInt local) {
        SInt global = 0;

        MPI_Reduce(&local, &global, 1, MPI_UINT64_T, MPI_MAX, 0, comm);

        return global;
    };

    const std::uint64_t active_sizes = reduce_sum_u64(instrumentation_.active_sizes);

    const std::uint64_t zero_count_sizes = reduce_sum_u64(instrumentation_.zero_count_sizes);

    const std::uint64_t planned_edges = reduce_sum_u64(instrumentation_.planned_edges);

    const std::uint64_t planned_pins = reduce_sum_u64(instrumentation_.planned_pins);

    const std::uint64_t generated_edges = reduce_sum_u64(instrumentation_.generated_edges);

    const std::uint64_t generated_pins = reduce_sum_u64(instrumentation_.generated_pins);

    const std::uint64_t attempts = reduce_sum_u64(instrumentation_.attempts);

    const std::uint64_t duplicate_rejects = reduce_sum_u64(instrumentation_.duplicate_rejects);

    const std::uint64_t minimum_samples = reduce_sum_u64(instrumentation_.minimum_samples);

    const std::uint64_t minimum_search_steps = reduce_sum_u64(instrumentation_.minimum_search_steps);

    const std::uint64_t minimum_cache_gets = reduce_sum_u64(instrumentation_.minimum_cache_gets);

    const std::uint64_t cache_map_calls = reduce_sum_u64(instrumentation_.cache_map_calls);

    const std::uint64_t cache_map_hits = reduce_sum_u64(instrumentation_.cache_map_hits);

    const std::uint64_t cache_map_misses = reduce_sum_u64(instrumentation_.cache_map_misses);

    const std::uint64_t cache_map_inserts = reduce_sum_u64(instrumentation_.cache_map_inserts);

    const std::uint64_t cache_candidate_calls = reduce_sum_u64(instrumentation_.cache_candidate_calls);

    const std::uint64_t cache_candidate_exact_hits = reduce_sum_u64(instrumentation_.cache_candidate_exact_hits);

    const std::uint64_t cache_candidate_recurrence_hits =
        reduce_sum_u64(instrumentation_.cache_candidate_recurrence_hits);

    const std::uint64_t cache_candidate_direct_evals = reduce_sum_u64(instrumentation_.cache_candidate_direct_evals);

    const std::uint64_t cache_candidate_below_range = reduce_sum_u64(instrumentation_.cache_candidate_below_range);

    const std::uint64_t cache_candidate_distance_sum = reduce_sum_u64(instrumentation_.cache_candidate_distance_sum);

    const std::uint64_t cache_candidate_max_distance = reduce_max_u64(instrumentation_.cache_candidate_max_distance);

    const std::uint64_t cache_candidate_distance_observations =
        reduce_sum_u64(instrumentation_.cache_candidate_distance_observations);

    const std::uint64_t cache_backward_steps = reduce_sum_u64(instrumentation_.cache_backward_steps);

    const std::uint64_t cache_forward_steps = reduce_sum_u64(instrumentation_.cache_forward_steps);
    const std::uint64_t cache_max_size      = reduce_max_u64(instrumentation_.cache_max_size);

    const std::uint64_t max_cache_initial_reserve = reduce_max_u64(instrumentation_.max_cache_initial_reserve);

    const std::uint64_t max_local_edges_for_size = reduce_max_u64(instrumentation_.max_local_edges_for_size);

    const SInt cache_min_key = reduce_min_sint(instrumentation_.cache_min_key);

    const SInt cache_max_key = reduce_max_sint(instrumentation_.cache_max_key);

    if (rank_ != 0) {
        return;
    }

    const auto safe_ratio = [](const std::uint64_t numerator, const std::uint64_t denominator) -> double {
        return denominator > 0 ? static_cast<double>(numerator) / static_cast<double>(denominator) : 0.0;
    };

    const std::uint64_t cache_total_calls = cache_map_calls + cache_candidate_calls;

    const std::uint64_t cache_total_direct_evals = cache_map_misses + cache_candidate_direct_evals;

    const std::uint64_t cache_candidate_hits = cache_candidate_exact_hits + cache_candidate_recurrence_hits;

    const std::uint64_t recurrence_steps = cache_backward_steps + cache_forward_steps;

    const bool has_cache_keys = cache_total_calls > 0 && cache_min_key != std::numeric_limits<SInt>::max();

    const std::uint64_t candidate_distance_observations =
        cache_candidate_calls >= cache_candidate_direct_evals
            ? cache_candidate_calls - cache_candidate_direct_evals
                  + (cache_candidate_direct_evals > 0 ? cache_candidate_direct_evals - 1 : 0)
            : 0;

    std::cout << "\nHGNP work statistics:\n";

    std::cout << "  Active size/rank entries:       " << active_sizes << '\n';
    std::cout << "  Zero-count size/rank entries:   " << zero_count_sizes << '\n';
    std::cout << "  Planned edges:                  " << planned_edges << '\n';
    std::cout << "  Generated edges:                " << generated_edges << '\n';
    std::cout << "  Planned pins:                   " << planned_pins << '\n';
    std::cout << "  Generated pins:                 " << generated_pins << '\n';
    std::cout << "  Maximum local edges for size:   " << max_local_edges_for_size << '\n';
    std::cout << "  Attempts:                       " << attempts << '\n';
    std::cout << "  Attempts per generated edge:    " << std::fixed << std::setprecision(6)
              << safe_ratio(attempts, generated_edges) << '\n';
    std::cout << "  Duplicate rejects:              " << duplicate_rejects << '\n';
    std::cout << "  Duplicate rejection rate:       " << safe_ratio(duplicate_rejects, attempts) << '\n';

    std::cout << "\nMinimum sampler:\n";
    std::cout << "  Samples:                        " << minimum_samples << '\n';
    std::cout << "  Search steps:                   " << minimum_search_steps << '\n';
    std::cout << "  Search steps per sample:        " << safe_ratio(minimum_search_steps, minimum_samples) << '\n';
    std::cout << "  Predicate cache gets:           " << minimum_cache_gets << '\n';
    std::cout << "  Predicate cache gets/sample:    " << safe_ratio(minimum_cache_gets, minimum_samples) << '\n';

    std::cout << "\nLogBinomCache map:\n";

    std::cout << "  Calls:                          " << cache_map_calls << '\n';

    std::cout << "  Hits:                           " << cache_map_hits << '\n';

    std::cout << "  Misses:                         " << cache_map_misses << '\n';

    std::cout << "  Hit rate:                       " << safe_ratio(cache_map_hits, cache_map_calls) << '\n';

    std::cout << "  Miss rate:                      " << safe_ratio(cache_map_misses, cache_map_calls) << '\n';

    std::cout << "  Inserts:                        " << cache_map_inserts << '\n';

    std::cout << "  Calls per generated edge:       " << safe_ratio(cache_map_calls, generated_edges) << '\n';

    std::cout << "  Misses per generated edge:      " << safe_ratio(cache_map_misses, generated_edges) << '\n';

    std::cout << "  Inserts per generated edge:     " << safe_ratio(cache_map_inserts, generated_edges) << '\n';

    std::cout << "\nLogBinomCache candidate cursor:\n";

    std::cout << "  Calls:                          " << cache_candidate_calls << '\n';

    std::cout << "  Exact cursor hits:              " << cache_candidate_exact_hits << '\n';

    std::cout << "  Recurrence hits:                " << cache_candidate_recurrence_hits << '\n';

    std::cout << "  Total cursor hit rate:          " << safe_ratio(cache_candidate_hits, cache_candidate_calls)
              << '\n';

    std::cout << "  Direct evaluations:             " << cache_candidate_direct_evals << '\n';

    std::cout << "  Direct evaluation rate:         " << safe_ratio(cache_candidate_direct_evals, cache_candidate_calls)
              << '\n';

    std::cout << "  Below-range requests:           " << cache_candidate_below_range << '\n';

    std::cout << "  Average cursor distance:        "
              << safe_ratio(cache_candidate_distance_sum, cache_candidate_distance_observations) << '\n';

    std::cout << "  Maximum cursor distance:        " << cache_candidate_max_distance << '\n';

    std::cout << "  Forward recurrence steps:       " << cache_forward_steps << '\n';

    std::cout << "  Backward recurrence steps:      " << cache_backward_steps << '\n';

    std::cout << "  Recurrence steps per hit:       " << safe_ratio(recurrence_steps, cache_candidate_recurrence_hits)
              << '\n';

    std::cout << "\nLogBinomCache storage:\n";

    std::cout << "  Total calls:                    " << cache_total_calls << '\n';

    std::cout << "  Total direct evaluations:       " << cache_total_direct_evals << '\n';

    std::cout << "  Maximum live entries:           " << cache_max_size << '\n';

    std::cout << "  Maximum reserved capacity:      " << max_cache_initial_reserve << '\n';

    std::cout << "  Maximum size/capacity ratio:    " << safe_ratio(cache_max_size, max_cache_initial_reserve) << '\n';

    if (has_cache_keys) {
        std::cout << "  Minimum accessed key:           " << cache_min_key << '\n';

        std::cout << "  Maximum accessed key:           " << cache_max_key << '\n';

        std::cout << "  Accessed key span:              " << cache_max_key - cache_min_key + 1 << '\n';
    } else {
        std::cout << "  Accessed key range:             n/a\n";
    }

    if (planned_edges != generated_edges) {
        std::cout << "  WARNING: planned/generated edge mismatch\n";
    }

    if (planned_pins != generated_pins) {
        std::cout << "  WARNING: planned/generated pin mismatch\n";
    }

#endif
}

// #### Size /Probability Setup ####
template <typename BigInt>
SInt HyperGNP<BigInt>::SampleExactEdgeCountNative(const SInt population, const double probability, const SInt seed) {
    if (population == 0 || probability <= 0.0) {
        return 0;
    }

    if (probability >= 1.0) {
        return population;
    }

    return rng_.GenerateBinomial(seed, population, probability);
}

template <typename BigInt>
SInt HyperGNP<BigInt>::SampleExactEdgeCountHuge(
    const CountInt& exact_population, const double probability, const SInt seed) {
    if (exact_population <= 0 || probability <= 0.0) {
        return 0;
    }

    if (exact_population <= std::numeric_limits<SInt>::max()) {
        return SampleExactEdgeCountNative(exact_population.convert_to<SInt>(), probability, seed);
    }

    if (probability >= 1.0) {
        throw ConfigurationError("HGNP local population with p=1 does not fit in SInt");
    }

    const long double population = exact_population.convert_to<long double>();

    if (!std::isfinite(population) || population <= 0.0L) {
        throw ConfigurationError("HGNP population cannot be represented as long double");
    }

    const long double log_population = std::log(population);
    const long double log_lambda     = log_population + std::log(static_cast<long double>(probability));

    // Keep the rest of the existing rejection sampler.

    if (!std::isfinite(population) || population <= 0.0L) {
        throw ConfigurationError("HGNP population cannot be represented as long double");
    }

    if (log_lambda <= std::log(static_cast<long double>(std::numeric_limits<double>::denorm_min()))) {
        return 0;
    }

    const long double log_sint_max = std::log(static_cast<long double>(std::numeric_limits<SInt>::max()));

    if (log_lambda > log_sint_max) {
        throw ConfigurationError("HGNP sampled edge count may exceed SInt");
    }

    const long double lambda_ld = std::exp(log_lambda);

    if (!std::isfinite(lambda_ld) || lambda_ld < 0.0L) {
        throw ConfigurationError("Invalid HGNP binomial expected edge count");
    }

    const double lambda = static_cast<double>(lambda_ld);

    /*
     * R(y + 1) / R(y) >= 1 exactly while y <= Np.
     * floor(Np) is therefore a maximizing mode of the
     * binomial-to-Poisson ratio.
     */
    const SInt mode = lambda_ld >= static_cast<long double>(std::numeric_limits<SInt>::max())
                          ? std::numeric_limits<SInt>::max()
                          : static_cast<SInt>(std::floor(lambda_ld));

    if (CountInt(mode) > exact_population) {
        throw ConfigurationError("HGNP binomial correction mode exceeds population");
    }

    Mersenne acceptance_rng;
    acceptance_rng.RandomInit(sampling::Spooky::hash(static_cast<unsigned long long>(seed) + 0x9e3779b97f4a7c15ULL));

    SInt proposal_seed = seed;

    constexpr SInt kMaximumAttempts = 1'000'000;

    for (SInt attempt = 0; attempt < kMaximumAttempts; ++attempt) {
        const SInt value = rng_.GeneratePoisson(sampling::Spooky::hash(proposal_seed++), lambda);

        /*
         * A sampled count cannot exceed the population.
         * This is practically impossible in the sparse regime,
         * but the condition is required for correctness.
         */
        if (CountInt(value) > exact_population) {
            continue;
        }

        long double log_acceptance =
            LogBinomialPoissonRatioRelativeToMode(value, mode, population, static_cast<long double>(probability));

        /*
         * Numerical roundoff can produce a tiny positive value
         * even though the theoretical ratio is at most one.
         */
        log_acceptance = std::min<long double>(log_acceptance, 0.0L);

        const long double uniform = MersenneUniform01(acceptance_rng);

        /*
         * MersenneUniform01 can return zero. In that case,
         * log(U) is -infinity and acceptance is automatic.
         */
        const long double log_uniform =
            uniform > 0.0L ? std::log(uniform) : -std::numeric_limits<long double>::infinity();

        if (log_uniform <= log_acceptance) {
            return value;
        }
    }

    throw ConfigurationError(
        "HGNP exact huge-population binomial rejection sampler "
        "exceeded its attempt limit");
}

template <typename BigInt>
void HyperGNP<BigInt>::PrepareSampledExactPlan(HGNPSizePlan& entry, const double probability) {
    const SInt k = entry.hyperedge_size;

    const auto [min_begin, min_end] = LocalMinOwnerRange(k, entry.partition_id);

    entry.range.begin   = min_begin;
    entry.range.end     = min_end;
    entry.range.local_m = 0;

    if (min_begin >= min_end || probability <= 0.0) {
        return;
    }

    CountInt local_population;
    {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        ScopedMPITimer timer(instrumentation_.population_seconds);
#endif

        local_population = MinRangeMassExact(min_begin, min_end, config_.n, k);
    }

    if (local_population == 0) {
        return;
    }

    SInt local_m;
    {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        ScopedMPITimer timer(instrumentation_.count_sampling_seconds);
#endif

        local_m = SampleExactEdgeCount(local_population, probability, LocalCountSeed(k, entry.partition_id));
    }

    entry.range.local_m = local_m;
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    instrumentation_.max_local_edges_for_size =
        std::max<std::uint64_t>(instrumentation_.max_local_edges_for_size, static_cast<std::uint64_t>(local_m));
#endif
}

template <typename BigInt>
SInt HyperGNP<BigInt>::SampleExactEdgeCount(const CountInt& population, const double probability, const SInt seed) {
    if (population <= std::numeric_limits<SInt>::max()) {
        return SampleExactEdgeCountNative(population.convert_to<SInt>(), probability, seed);
    }

    return SampleExactEdgeCountHuge(population, probability, seed);
}

template <typename BigInt>
long double HyperGNP<BigInt>::LogBinomialPoissonRatioRelativeToMode(
    const SInt value, const SInt mode, const long double population, const long double probability) const {
    if (!(population > 0.0L) || !(probability > 0.0L) || !(probability < 1.0L)) {
        throw ConfigurationError("Invalid parameters for HGNP binomial count correction");
    }

    const long double log_q = std::log1pl(-probability);

    long double log_ratio = 0.0L;

    if (value > mode) {
        for (SInt j = mode; j < value; ++j) {
            const long double fraction = static_cast<long double>(j) / population;

            /*
             * R(j + 1) / R(j)
             *     = (1 - j / N) / (1 - p).
             */
            log_ratio += std::log1pl(-fraction) - log_q;
        }

        return log_ratio;
    }

    for (SInt j = value; j < mode; ++j) {
        const long double fraction = static_cast<long double>(j) / population;

        /*
         * Moving from mode down to value requires the inverse
         * of each forward ratio.
         */
        log_ratio -= std::log1pl(-fraction) - log_q;
    }

    return log_ratio;
}

template <typename BigInt>
std::vector<HGNPSizePlan> HyperGNP<BigInt>::BuildGenerationPlan(const SInt lower_bound, const SInt upper_bound) {
    std::vector<HGNPSizePlan> plan;

    /*
     * EdgeAndPinBudget produces a sparse set of nonzero size counts.
     * Iterating through every k in [lower_bound, upper_bound] would
     * take Theta(n) time and reserve Theta(n) plan entries.
     */
    if (probs_type_ == ProbabilityMode::EdgeAndPinBudget) {
        std::vector<SInt> active_sizes;
        active_sizes.reserve(budget_size_counts_.size());

        for (const auto& [k, count]: budget_size_counts_) {
            if (count > 0 && k >= lower_bound && k <= upper_bound) {
                active_sizes.push_back(k);
            }
        }

        /*
         * budget_size_counts_ is apparently unordered. Sorting preserves
         * the previous increasing-k generation and CSR order.
         */
        std::sort(active_sizes.begin(), active_sizes.end());

        plan.reserve(active_sizes.size());

        for (const SInt k: active_sizes) {
            /*
             * This reuses the existing probability resolution, exact-count
             * sampling, instrumentation and plan insertion logic.
             */
            const bool continue_planning = AppendSizePlanIfNeeded(k, lower_bound, plan);

            if (!continue_planning) {
                /*
                 * This should not currently happen for EdgeAndPinBudget,
                 * but retaining the check makes the contract explicit.
                 */
                break;
            }
        }

        return plan;
    }

    /*
     * Existing dense/terminating behavior for all other probability modes.
     */
    if (upper_bound >= lower_bound) {
        const SInt number_of_sizes = upper_bound - lower_bound + 1;

        plan.reserve(static_cast<std::size_t>(number_of_sizes));
    }

    // Avoid overflowing k when upper_bound == SInt::max().
    for (SInt k = lower_bound;; ++k) {
        if (!AppendSizePlanIfNeeded(k, lower_bound, plan)) {
            break;
        }

        if (k == upper_bound) {
            break;
        }
    }

    return plan;
}
template <typename BigInt>
bool HyperGNP<BigInt>::AppendSizePlanIfNeeded(
    const SInt hyperedge_size, const SInt lower_bound, std::vector<HGNPSizePlan>& plan) {
    const SizeGenerationParameters params = GetSizeGenerationParameters(hyperedge_size, lower_bound);

    if (ExpectedCountSkipsRemainingSizes(params)) {
        return false;
    }

    if (ExpectedCountSkipsCurrentSize(params)) {
        return true;
    }

    double probability = params.probability;

    if (params.expected_count > 0.0) {
        const auto resolved = ResolveProbabilityForSize(hyperedge_size, params);

        if (!resolved) {
            return probs_type_ != ProbabilityMode::EdgeBudget;
        }

        probability = *resolved;
    }

    ValidateProbability(probability);

    if (ShouldSkipSizeGeneration(hyperedge_size, probability)) {
        return true;
    }

    const auto partitions = AssignPartitionsToPE(config_.k, rank_, size_);

    for (SInt partition_id = partitions.begin; partition_id < partitions.end; ++partition_id) {
        HGNPSizePlan entry = PrepareSizePlan(hyperedge_size, probability, partition_id);

        if (entry.range.local_m > 0) {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
            ++instrumentation_.active_sizes;

            instrumentation_.planned_edges += static_cast<std::uint64_t>(entry.range.local_m);

            instrumentation_.planned_pins +=
                static_cast<std::uint64_t>(entry.range.local_m) * static_cast<std::uint64_t>(hyperedge_size);
#endif

            plan.push_back(std::move(entry));
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        } else {
            ++instrumentation_.zero_count_sizes;
#endif
        }
    }
    return true;
}

template <typename BigInt>
HGNPSizePlan HyperGNP<BigInt>::PrepareSizePlan(const SInt hyperedge_size, const double probability, SInt partition_id) {
    HGNPSizePlan entry;
    entry.hyperedge_size = hyperedge_size;
    entry.partition_id   = partition_id;

    if (config_.approx) {
        entry.range = PrepareApproxLocalRange(hyperedge_size, probability, partition_id);
    } else {
        PrepareSampledExactPlan(entry, probability);
    }

    return entry;
}

template <typename BigInt>
SInt HyperGNP<BigInt>::LocalCountSeed(const SInt hyperedge_size, const SInt partition_id) const {
    return sampling::Spooky::hash(
        static_cast<unsigned long long>(config_.seed)
        + (static_cast<unsigned long long>(hyperedge_size) * kCountSeedMultiplier)
        + (static_cast<unsigned long long>(partition_id) * kRankSeedMultiplier));
}

template <typename BigInt>
HGNPLocalGenerationRange
HyperGNP<BigInt>::PrepareApproxLocalRange(const SInt hyperedge_size, const double probability, SInt partition_id) {
    const auto [begin, end] = LocalMinOwnerRange(hyperedge_size, partition_id);

    HGNPLocalGenerationRange range;
    range.begin = begin;
    range.end   = end;

    if (begin >= end || probability <= 0.0) {
        return range;
    }

    LogBinomCache cache(hyperedge_size);

    const long double log_total = cache.Get(config_.n, hyperedge_size);

    const long double mass = MinRangeMassApproxCached(begin, end, config_.n, hyperedge_size, cache);

    if (mass <= 0.0L) {
        return range;
    }

    const long double log_lambda = std::log(static_cast<long double>(probability)) + log_total + std::log(mass);

    if (log_lambda > std::log(static_cast<long double>(std::numeric_limits<SInt>::max()))) {
        throw ConfigurationError("HGNP approximate local edge count exceeds SInt");
    }

    const double lambda = static_cast<double>(std::exp(log_lambda));

    range.local_m = lambda > 0.0 ? rng_.GeneratePoisson(LocalCountSeed(hyperedge_size, partition_id), lambda) : 0;

    return range;
}

template <typename BigInt>
void HyperGNP<BigInt>::LoadSizeDistributionInputs() {
    LoadSizeProbabilitiesFromFile();
    LoadExpectedSizeCountsFromFile();
}

template <typename BigInt>
void HyperGNP<BigInt>::LoadSizeProbabilitiesFromFile() {
    if (config_.size_probabilities_file.empty()) {
        return;
    }

    const auto values = ReadNumericLines(config_.size_probabilities_file, "Could not open HGNP size probability file");

    config_.size_probabilities.insert(config_.size_probabilities.end(), values.begin(), values.end());
}

template <typename BigInt>
void HyperGNP<BigInt>::LoadExpectedSizeCountsFromFile() {
    if (config_.size_expected_counts_file.empty()) {
        return;
    }

    const auto values =
        ReadNumericLines(config_.size_expected_counts_file, "Could not open HGNP size expected counts file");

    config_.size_expected_counts.insert(config_.size_expected_counts.end(), values.begin(), values.end());
}

template <typename BigInt>
void HyperGNP<BigInt>::SelectProbabilityMode() {
    const bool explicit_probs       = !config_.size_probabilities.empty();
    const bool explicit_expected    = !config_.size_expected_counts.empty();
    const bool edge_budget_mode     = config_.edge_budget > 0.0;
    const bool edge_pin_budget_mode = config_.edge_budget > 0.0 && config_.size_dist_pin_budget > 0;
    if (explicit_probs) {
        probs_type_ = ProbabilityMode::ExplicitProbabilities;
        return;
    }
    if (explicit_expected) {
        probs_type_ = ProbabilityMode::ExplicitExpectedCounts;
        return;
    }
    if (edge_pin_budget_mode) {
        probs_type_ = ProbabilityMode::EdgeAndPinBudget;
        return;
    }
    if (edge_budget_mode) {
        probs_type_ = ProbabilityMode::EdgeBudget;
        return;
    }

    probs_type_ = ProbabilityMode::GlobalProbability;
}

template <typename BigInt>
SizeGenerationParameters
HyperGNP<BigInt>::GetSizeGenerationParameters(const SInt hyperedge_size, const SInt lower_bound) {
    SizeGenerationParameters params;

    switch (probs_type_) {
        case ProbabilityMode::ExplicitProbabilities:
            params.probability = config_.size_probabilities[hyperedge_size - lower_bound];
            break;

        case ProbabilityMode::ExplicitExpectedCounts:
            params.expected_count = config_.size_expected_counts[hyperedge_size - lower_bound];
            break;

        case ProbabilityMode::EdgeBudget:
            params.expected_count =
                config_.edge_budget * config_.size_decay
                * std::pow(1.0 - config_.size_decay, static_cast<double>(hyperedge_size - lower_bound));
            break;

        case ProbabilityMode::EdgeAndPinBudget: {
            auto it = budget_size_counts_.find(hyperedge_size);
            if (it != budget_size_counts_.end()) {
                params.expected_count = static_cast<double>(it->second);
            } else {
                params.expected_count = 0.0;
            }
            break;
        }

        case ProbabilityMode::GlobalProbability:
            params.probability = config_.p;
            break;
    }

    return params;
}

// Convert an expected number of hyperedges of this size into the
// equivalent independent edge probability p = expected / C(n, k).
// Returns nullopt if p underflows to zero.
template <typename BigInt>
std::optional<double>
HyperGNP<BigInt>::ResolveProbabilityForSize(SInt hyperedge_size, const SizeGenerationParameters& params) const {
    if (params.expected_count > 0.0) {
        const long double log_p =
            std::log(static_cast<long double>(params.expected_count)) - LogBinomialApprox(config_.n, hyperedge_size);

        if (log_p < std::log(static_cast<long double>(std::numeric_limits<double>::denorm_min()))) {
            return std::nullopt;
        }

        return std::clamp(static_cast<double>(expl(log_p)), 0.0, 1.0);
    }

    return params.probability;
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateBudgetSizeCounts(const SInt lower_bound, const SInt upper_bound) {
    budget_size_counts_.clear();

    if (config_.edge_budget <= 0.0 || config_.size_dist_pin_budget <= 0) {
        return;
    }

    const SInt edge_budget = static_cast<SInt>(std::llround(config_.edge_budget));

    const auto min_pins = static_cast<unsigned __int128>(lower_bound) * static_cast<unsigned __int128>(edge_budget);
    const auto max_pins = static_cast<unsigned __int128>(upper_bound) * static_cast<unsigned __int128>(edge_budget);
    const auto pins     = static_cast<unsigned __int128>(config_.size_dist_pin_budget);

    if (pins < min_pins || pins > max_pins) {
        throw ConfigurationError("Invalid HGNP pin budget");
    }

    BoltzmannPinBudgetSizeCountSampler sampler(edge_budget, lower_bound, upper_bound, config_.size_dist_pin_budget);

    budget_size_counts_ = sampler.Sample();
}

template <typename BigInt>
SInt HyperGNP<BigInt>::EffectiveUpperBound(SInt hard_upper_bound, SInt lower_bound) {
    switch (probs_type_) {
        case ProbabilityMode::ExplicitProbabilities:
            return std::min<SInt>(lower_bound + static_cast<SInt>(config_.size_probabilities.size()) - 1, config_.n);
        case ProbabilityMode::ExplicitExpectedCounts:
            return std::min<SInt>(lower_bound + static_cast<SInt>(config_.size_expected_counts.size()) - 1, config_.n);
        default:
            return hard_upper_bound;
    }
}

template <typename BigInt>
bool HyperGNP<BigInt>::UsesExpectedCount() const {
    return probs_type_ == ProbabilityMode::ExplicitExpectedCounts || probs_type_ == ProbabilityMode::EdgeBudget
           || probs_type_ == ProbabilityMode::EdgeAndPinBudget;
}

template <typename BigInt>
bool HyperGNP<BigInt>::ExpectedCountSkipsRemainingSizes(const SizeGenerationParameters& params) const {
    switch (probs_type_) {
        case ProbabilityMode::EdgeBudget:
            // Geometric decay produces monotonically decreasing counts, so once
            // we hit zero we can stop.
            return ExpectedCountSkipsCurrentSize(params);

        case ProbabilityMode::EdgeAndPinBudget:
            // Boltzmann-generated counts are not guaranteed to be contiguous.
            // A later size may still receive a nonzero count.
            return false;

        default:
            return false;
    }
}

template <typename BigInt>
bool HyperGNP<BigInt>::ExpectedCountSkipsCurrentSize(const SizeGenerationParameters& params) const {
    return UsesExpectedCount() && (params.expected_count == 0.0 || ExpectedCountIsNegligible(params.expected_count));
}

template <typename BigInt>
void HyperGNP<BigInt>::ValidateProbability(const double probability) const {
    if (probability < 0.0 || probability > 1.0) {
        throw ConfigurationError("HGNP probability must be in [0, 1]");
    }
}

template <typename BigInt>
bool HyperGNP<BigInt>::ShouldSkipSizeGeneration(const SInt hyperedge_size, const double probability) const {
    return hyperedge_size == 0 || hyperedge_size > config_.n || probability <= 0.0;
}

// #### Partition helpers ####
template <typename BigInt>
void HyperGNP<BigInt>::SetLocalVertexRange() {
    const auto [start, end] = BalancedVertexRange(config_.n, rank_, size_);
    SetVertexRange(start, end);
}

template <typename BigInt>
std::pair<SInt, SInt> HyperGNP<BigInt>::LocalMinOwnerRange(const SInt hyperedge_size, const SInt partition_id) const {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    const double start = MPI_Wtime();
#endif
    const SInt begin = FindMinBoundaryByMass(config_.n, hyperedge_size, partition_id, config_.k);

    const SInt end = FindMinBoundaryByMass(config_.n, hyperedge_size, partition_id + 1, config_.k);
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    instrumentation_.boundary_seconds += MPI_Wtime() - start;
#endif
    return {begin, end};
}

// #### Approx generation ####

template <typename BigInt>
void HyperGNP<BigInt>::GenerateLocalHyperedges(
    const HGNPSizePlan& entry, SInt& edge_seed, LogBinomCache& cache, std::vector<SInt>& pins, HyperedgeSeenSet& seen) {
    const SInt local_m = entry.range.local_m;

    const std::uint64_t local_m_u64 = static_cast<std::uint64_t>(local_m);

    const std::uint64_t max_attempts = local_m_u64 > std::numeric_limits<std::uint64_t>::max() / 10
                                           ? std::numeric_limits<std::uint64_t>::max()
                                           : std::max<std::uint64_t>(local_m_u64 * 10, 1000);

    SInt generated = 0;

    while (generated < local_m) {
        const auto event_start = std::chrono::steady_clock::now();

        std::uint64_t sampling_attempts    = 0;
        std::uint64_t duplicate_rejections = 0;
        std::uint64_t minimum_search_steps = 0;
        std::uint64_t minimum_cache_gets   = 0;

        while (true) {
            std::uint64_t attempt_search_steps = 0;
            std::uint64_t attempt_cache_gets   = 0;

            SampleLocalHyperedgeInto(
                entry.hyperedge_size, entry.range.begin, entry.range.end, edge_seed, cache, pins, attempt_search_steps,
                attempt_cache_gets);

            ++sampling_attempts;
            minimum_search_steps += attempt_search_steps;
            minimum_cache_gets += attempt_cache_gets;

            bool accepted = true;

            if (!config_.allow_duplicates) {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
                const double duplicate_start = MPI_Wtime();
#endif

                accepted = seen.insert(FingerprintPins(pins)).second;

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
                instrumentation_.duplicate_check_seconds += MPI_Wtime() - duplicate_start;
#endif
            }

            if (!accepted) {
                ++duplicate_rejections;

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
                ++instrumentation_.duplicate_rejects;
#endif

                if (sampling_attempts > max_attempts) {
                    throw ConfigurationError("HGNP sampling exceeded attempt limit");
                }

                continue;
            }

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
            const double write_start = MPI_Wtime();
#endif

            PushUncompressedHyperedge(graph_, memory_stats_, pins);

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
            instrumentation_.csr_write_seconds += MPI_Wtime() - write_start;
            ++instrumentation_.generated_edges;
            instrumentation_.generated_pins += static_cast<std::uint64_t>(pins.size());
            instrumentation_.attempts += sampling_attempts;
#endif

            ++generated;

            if (debug_logger_) {
                const auto duration_ns =
                    std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - event_start)
                        .count();

                debug_logger_->LogHyperedge({
                    .hyperedge_id         = next_debug_hyperedge_id_++,
                    .hyperedge_size       = entry.hyperedge_size,
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
}

template <typename BigInt>
void HyperGNP<BigInt>::SampleLocalHyperedgeInto(
    const SInt hyperedge_size, const SInt local_min_begin, const SInt local_min_end, SInt& edge_seed,
    LogBinomCache& log_binom_cache, std::vector<SInt>& pins, std::uint64_t& minimum_search_steps,
    std::uint64_t& minimum_cache_gets) {
    pins.clear();

    minimum_search_steps = 0;
    minimum_cache_gets   = 0;

    SInt minimum_vertex;

    {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        ScopedMPITimer timer(instrumentation_.minimum_sample_seconds);
#endif
        minimum_vertex = SampleMinimumImplicit(
            local_min_begin, local_min_end, config_.n, hyperedge_size, rng_, edge_seed, log_binom_cache,
            &minimum_search_steps, &minimum_cache_gets);
    }
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    ++instrumentation_.minimum_samples;
    instrumentation_.minimum_search_steps += minimum_search_steps;
    instrumentation_.minimum_cache_gets += minimum_cache_gets;
#endif
    pins.push_back(minimum_vertex);

    {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        ScopedMPITimer timer(instrumentation_.pin_sample_seconds);
#endif
        FloydSampleInto(
            minimum_vertex + 1, config_.n - minimum_vertex - 1, hyperedge_size - 1, mersenne_, pins, floyd_scratch_, 1);
    }
}
template <typename BigInt>
SInt HyperGNP<BigInt>::LocalEdgeSeed(const SInt hyperedge_size, const SInt partition_id) const {
    return sampling::Spooky::hash(
        static_cast<unsigned long long>(config_.seed)
        + (static_cast<unsigned long long>(partition_id) * kEdgeRankSeedMultiplier)
        + (static_cast<unsigned long long>(hyperedge_size) * kEdgeSeedMultiplier));
}

// #### Logging/graph mutation ####
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
template class HyperGNP<__int128>;
template class HyperGNP<SInt>;
#pragma GCC diagnostic pop
} // namespace kagen
