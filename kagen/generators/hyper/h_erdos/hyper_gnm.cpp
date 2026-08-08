#include "kagen/generators/hyper/h_erdos/hyper_gnm.h"

#include "kagen/generators/generator.h"
#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"
#include "kagen/tools/rng_wrapper.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace kagen {

std::unique_ptr<Generator>
HyperGNMFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    if (config.n <= (1ull << 31)) {
        return std::make_unique<HyperGNMSmall>(config, rank, size);
    }
    return std::make_unique<HyperGNMBig>(config, rank, size);
}

void ValidateHGNMConfig(const PGeneratorConfig& config) {
    if (config.n == 0) {
        throw ConfigurationError("HGNM requires n > 0");
    }

    if (config.m == 0) {
        throw ConfigurationError("HGNM requires m > 0");
    }

    if (config.k == 0) {
        throw ConfigurationError("HGNM requires at least one chunk");
    }
}

template <typename BigInt>
HyperGNM<BigInt>::HyperGNM(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size),
      rng_(config) {
    ValidateHGNMConfig(config);

    if (config_.debug) {
        debug_logger_.emplace(
            config_.output_graph.filename + "_hgnm_debug_pe_" + std::to_string(rank_) + ".csv", false);
    }
}

template <typename BigInt>
std::pair<SInt, SInt> HyperGNM<BigInt>::ResolveSizeRange() const {
    const SInt lower_bound = config_.size_dist_lower_bound;
    const SInt hard_upper_bound =
        config_.size_dist_upper_bound != 0 ? std::min(config_.size_dist_upper_bound, config_.n) : config_.n;

    if (lower_bound > hard_upper_bound) {
        throw ConfigurationError("Invalid hyperedge size range");
    }

    return {lower_bound, hard_upper_bound};
}

template <typename BigInt>
double HyperGNM<BigInt>::ValidateAndGetSizeAlpha() const {
    const double alpha = config_.size_dist_alpha;
    if (alpha <= 0.0 || alpha > 1.0) {
        throw ConfigurationError("Invalid size distribution alpha");
    }

    return alpha;
}

template <typename BigInt>
void HyperGNM<BigInt>::GenerateSizeCountsForConfig(
    SInt lower_bound, SInt upper_bound, double alpha, std::unordered_map<SInt, SInt>& size_counts) {
    if (config_.size_dist_pin_budget > 0) {
        GenerateBoltzmannPinBudgetSizeCounts(lower_bound, upper_bound, config_.size_dist_pin_budget, size_counts);
    } else if (config_.size_dist_deterministic) {
        GenerateDeterministicGeometricSizeCounts(lower_bound, upper_bound, config_.size_decay, size_counts);
    } else {
        GenerateRandomGeometricSizeCounts(lower_bound, upper_bound, alpha, size_counts);
    }
}
template <typename BigInt>
std::pair<SInt, SInt> HyperGNM<BigInt>::ComputeLocalVertexRange() const {
    return BalancedVertexRange(config_.n, rank_, size_);
}

template <typename BigInt>
void HyperGNM<BigInt>::LogSizeCounts(const std::unordered_map<SInt, SInt>& size_counts) const {
    if (!config_.debug) {
        return;
    }

    std::vector<std::pair<SInt, SInt>> entries(size_counts.begin(), size_counts.end());
    std::sort(entries.begin(), entries.end());

    std::cout << "size_counts:\n";
    for (const auto& [size, count]: entries) {
        std::cout << "  size=" << size << " count=" << count << '\n';
    }
}

template <typename BigInt>
void HyperGNM<BigInt>::GenerateHyperedgesFromPlan(const std::vector<HGNMSizePlan>& plan) {
    for (const auto& entry: plan) {
        GenerateHyperedgesOfSize(entry.hyperedge_size, entry.m_k, entry.partition_id, entry.range);
    }
}

template <typename BigInt>
void HyperGNM<BigInt>::GenerateCSR() {
    auto [lower_bound, hard_upper_bound] = ResolveSizeRange();
    double alpha                         = ValidateAndGetSizeAlpha();

    std::unordered_map<SInt, SInt> size_counts;

    GenerateSizeCountsForConfig(lower_bound, hard_upper_bound, alpha, size_counts);

    LogSizeCounts(size_counts);

    std::vector<HGNMSizePlan> plan;
    plan.reserve(size_counts.size());

    std::size_t expected_local_edges = 0;
    std::size_t expected_local_pins  = 0;

    const auto partitions = AssignPartitionsToPE(config_.k, rank_, size_);

    for (const auto& [k, m_k]: size_counts) {
        for (SInt partition = partitions.begin; partition < partitions.end; ++partition) {
            auto range = PrepareLocalGenerationRange(k, m_k, partition);

            expected_local_edges += static_cast<std::size_t>(range.local_m);

            expected_local_pins += static_cast<std::size_t>(range.local_m) * static_cast<std::size_t>(k);

            plan.push_back({
                .hyperedge_size = k,
                .m_k            = m_k,
                .partition_id   = partition,
                .range          = range,
            });
        }
    }

    graph_.hyperedge_offsets.reserve(expected_local_edges + 1);
    graph_.hyperedge_pins.reserve(expected_local_pins);
    graph_.hyperedge_offsets.push_back(0);

    GenerateHyperedgesFromPlan(plan);

    auto [start, end] = ComputeLocalVertexRange();
    SetVertexRange(start, end);
}

template <typename BigInt>
void HyperGNM<BigInt>::GenerateBoltzmannPinBudgetSizeCounts(
    const SInt lower_bound, const SInt upper_bound, const SInt pin_budget,
    std::unordered_map<SInt, SInt>& size_counts) {
    size_counts.clear();

    const SInt m = config_.m;

    if (lower_bound > upper_bound) {
        throw ConfigurationError("Invalid hyperedge size range");
    }

    const auto pins     = static_cast<unsigned __int128>(pin_budget);
    const auto min_pins = static_cast<unsigned __int128>(lower_bound) * static_cast<unsigned __int128>(m);
    const auto max_pins = static_cast<unsigned __int128>(upper_bound) * static_cast<unsigned __int128>(m);

    if (pins < min_pins || pins > max_pins) {
        throw ConfigurationError("Invalid HGNM pin budget");
    }

    BoltzmannPinBudgetSizeCountSampler sampler(m, lower_bound, upper_bound, pin_budget);

    size_counts = sampler.Sample();
}

template <typename BigInt>
void HyperGNM<BigInt>::GenerateRandomGeometricSizeCounts(
    const SInt lower, const SInt upper, const double alpha, std::unordered_map<SInt, SInt>& size_counts) {
    SInt remaining_m = config_.m;

    long double remaining_mass = 1.0L;

    SInt seed = sampling::Spooky::hash(static_cast<SInt>(config_.seed) * kGNMCountSeedMultiplier);

    const long double q     = 1.0L - static_cast<long double>(alpha);
    const SInt        max_t = upper - lower;

    const long double normalizer = q == 0.0L ? 1.0L : 1.0L - std::pow(q, static_cast<long double>(max_t + 1));

    for (SInt k = lower; k <= upper && remaining_m > 0; ++k) {
        const SInt t = k - lower;

        const long double pk =
            q == 0.0L ? (k == lower ? 1.0L : 0.0L)
                      : (static_cast<long double>(alpha) * std::pow(q, static_cast<long double>(t))) / normalizer;

        const long double conditional_p = std::clamp(pk / remaining_mass, 0.0L, 1.0L);

        SInt count = 0;

        if (k == upper) {
            count = remaining_m;
        } else {
            const SInt draw_seed = sampling::Spooky::hash(seed++);
            count                = rng_.GenerateBinomial(draw_seed, remaining_m, conditional_p);
        }

        if (count > 0) {
            size_counts[k] = count;
        }

        remaining_m -= count;
        remaining_mass -= pk;
    }
}

template <typename BigInt>
void HyperGNM<BigInt>::GenerateDeterministicGeometricSizeCounts(
    const SInt lower, const SInt upper, const double decay, std::unordered_map<SInt, SInt>& size_counts) {
    if (decay <= 0.0 || decay > 1.0) {
        throw ConfigurationError("Invalid deterministic size decay");
    }

    SInt remaining_m = config_.m;

    long double remaining_factor = 1.0L;

    for (SInt k = lower; k <= upper && remaining_m > 0; ++k) {
        long double expected = static_cast<long double>(config_.m) * static_cast<long double>(decay) * remaining_factor;

        SInt count = static_cast<SInt>(std::floor(expected));

        if (count <= 0 || k == upper) {
            count = remaining_m;
        }

        count = std::min(count, remaining_m);

        if (count > 0) {
            size_counts[k] = count;
        }

        remaining_m -= count;
        remaining_factor *= (1.0L - static_cast<long double>(decay));
    }
}

template <typename BigInt>
void HyperGNM<BigInt>::FinalizeCSR(MPI_Comm comm) {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
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

    const std::uint64_t map_calls = reduce_sum_u64(instrumentation_.cache_map_calls);

    const std::uint64_t map_hits = reduce_sum_u64(instrumentation_.cache_map_hits);

    const std::uint64_t map_misses = reduce_sum_u64(instrumentation_.cache_map_misses);

    const std::uint64_t map_inserts = reduce_sum_u64(instrumentation_.cache_map_inserts);

    const std::uint64_t candidate_calls = reduce_sum_u64(instrumentation_.cache_candidate_calls);

    const std::uint64_t candidate_exact_hits = reduce_sum_u64(instrumentation_.cache_candidate_exact_hits);

    const std::uint64_t candidate_recurrence_hits = reduce_sum_u64(instrumentation_.cache_candidate_recurrence_hits);

    const std::uint64_t candidate_direct_evals = reduce_sum_u64(instrumentation_.cache_candidate_direct_evals);

    const std::uint64_t candidate_distance_sum = reduce_sum_u64(instrumentation_.cache_candidate_distance_sum);

    const std::uint64_t candidate_distance_observations =
        reduce_sum_u64(instrumentation_.cache_candidate_distance_observations);

    const std::uint64_t candidate_max_distance = reduce_max_u64(instrumentation_.cache_candidate_max_distance);

    const std::uint64_t max_size = reduce_max_u64(instrumentation_.cache_max_size);

    const std::uint64_t max_capacity = reduce_max_u64(instrumentation_.max_cache_initial_reserve);

    const std::uint64_t candidate_below_range = reduce_sum_u64(instrumentation_.cache_candidate_below_range);

    const std::uint64_t backward_steps = reduce_sum_u64(instrumentation_.cache_backward_steps);

    const std::uint64_t forward_steps = reduce_sum_u64(instrumentation_.cache_forward_steps);

    const std::uint64_t recurrence_steps = forward_steps + backward_steps;

    if (rank_ != 0) {
        return;
    }

    const auto ratio = [](const std::uint64_t numerator, const std::uint64_t denominator) {
        return denominator > 0 ? static_cast<double>(numerator) / static_cast<double>(denominator) : 0.0;
    };

    std::cout << "\nHGNM LogBinomCache statistics:\n"
              << "  Map calls:                 " << map_calls << '\n'
              << "  Map hits:                  " << map_hits << '\n'
              << "  Map misses:                " << map_misses << '\n'
              << "  Map hit rate:              " << ratio(map_hits, map_calls) << '\n'
              << "  Map inserts:               " << map_inserts << '\n'
              << "  Candidate calls:           " << candidate_calls << '\n'
              << "  Candidate exact hits:      " << candidate_exact_hits << '\n'
              << "  Candidate recurrence hits: " << candidate_recurrence_hits << '\n'
              << "  Candidate direct evals:    " << candidate_direct_evals << '\n'
              << "  Candidate direct rate:     " << ratio(candidate_direct_evals, candidate_calls) << '\n'
              << "  Average candidate distance:" << ratio(candidate_distance_sum, candidate_distance_observations)
              << '\n'
              << "  Maximum candidate distance:" << candidate_max_distance << '\n'
              << "  Maximum live entries:      " << max_size << '\n'
              << "  Maximum requested capacity:" << max_capacity << '\n'
              << "  Candidate below range:     " << candidate_below_range << '\n'

              << "  Forward recurrence steps:  " << forward_steps << '\n'

              << "  Backward recurrence steps: " << backward_steps << '\n'
              << "  Steps per recurrence hit:  " << ratio(recurrence_steps, candidate_recurrence_hits) << '\n';

#endif
}

template <typename BigInt>
void HyperGNM<BigInt>::ValidateExactSparseDensity(
    SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, SInt local_m) const {
    const CountInt local_population = MinRangeMassExact(local_min_begin, local_min_end, config_.n, hyperedge_size);

    const long double density = static_cast<long double>(local_m) / local_population.convert_to<long double>();

    if (density > 0.25L) {
        throw ConfigurationError("Dense exact HGNM not implemented yet");
    }

    if (CountInt(local_m) > local_population) {
        throw ConfigurationError("Exact HGNM local_m exceeds local population");
    }
}

template <typename BigInt>
SInt HyperGNM<BigInt>::EdgeSeed(SInt hyperedge_size) const {
    return sampling::Spooky::hash(
        static_cast<unsigned long long>(config_.seed)
        + (static_cast<unsigned long long>(rank_) * kEdgeRankSeedMultiplier)
        + (static_cast<unsigned long long>(hyperedge_size) * kEdgeSeedMultiplier));
}
template <typename BigInt>
void HyperGNM<BigInt>::SampleHyperedgeInto(
    const SInt minimum_vertex, const SInt hyperedge_size, std::vector<SInt>& pins) {
    pins.clear();
    pins.push_back(minimum_vertex);

    FloydSampleInto(
        minimum_vertex + 1, config_.n - minimum_vertex - 1, hyperedge_size - 1, mersenne_, pins, floyd_scratch_, 1);
}

template <typename BigInt>
bool HyperGNM<BigInt>::TryPushHyperedge(const std::vector<SInt>& pins, HyperedgeSeenSet& local_seen) {
    if (config_.allow_duplicates || local_seen.insert(FingerprintPins(pins)).second) {
        PushUncompressedHyperedge(graph_, memory_stats_, pins);
        return true;
    }

    return false;
}

template <typename BigInt>
void HyperGNM<BigInt>::GenerateSampledHyperedges(
    const SInt hyperedge_size, const SInt local_min_begin, const SInt local_min_end, const SInt local_m,
    SInt& edge_seed, LogBinomCache& log_binom_cache) {
    auto local_seen = MakeLocalSeenSet(config_.allow_duplicates, local_m);

    std::vector<SInt> pins;
    pins.reserve(static_cast<std::size_t>(hyperedge_size));

    const std::uint64_t max_attempts = std::max<std::uint64_t>(static_cast<std::uint64_t>(local_m) * 10, 1000);

    std::uint64_t total_attempts = 0;
    SInt          generated      = 0;

    while (generated < local_m) {
        const auto event_start = std::chrono::steady_clock::now();

        std::uint64_t sampling_attempts    = 0;
        std::uint64_t duplicate_rejections = 0;
        std::uint64_t minimum_search_steps = 0;
        std::uint64_t minimum_cache_gets   = 0;

        while (true) {
            std::uint64_t attempt_search_steps = 0;
            std::uint64_t attempt_cache_gets   = 0;

            const SInt minimum_vertex = SampleMinimumImplicit(
                local_min_begin, local_min_end, config_.n, hyperedge_size, rng_, edge_seed, log_binom_cache,
                &attempt_search_steps, &attempt_cache_gets);

            SampleHyperedgeInto(minimum_vertex, hyperedge_size, pins);

            ++sampling_attempts;
            ++total_attempts;

            minimum_search_steps += attempt_search_steps;
            minimum_cache_gets += attempt_cache_gets;

            if (!TryPushHyperedge(pins, local_seen)) {
                ++duplicate_rejections;

                if (total_attempts > max_attempts) {
                    throw ConfigurationError("HGNM rejection sampling exceeded attempt limit");
                }

                continue;
            }

            ++generated;

            if (debug_logger_) {
                const auto duration_ns =
                    std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - event_start)
                        .count();

                debug_logger_->LogHyperedge({
                    .hyperedge_id         = next_debug_hyperedge_id_++,
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
}

template <typename BigInt>
HGNMLocalGenerationRange
HyperGNM<BigInt>::PrepareLocalGenerationRange(const SInt hyperedge_size, const SInt m_k, const SInt partition_id) {
    if (!config_.approx) {
        ExactFixedCountHyperedgeGenerator<BigInt> exact_generator(
            config_, partition_id, config_.k, rng_, mersenne_, graph_, memory_stats_, floyd_scratch_,
            debug_logger_ ? &*debug_logger_ : nullptr, &next_debug_hyperedge_id_
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
            ,
            &instrumentation_
#endif
        );

        const FixedCountLocalRange exact_range = exact_generator.PrepareLocalRange(hyperedge_size, m_k);

        return {
            .min_begin  = exact_range.min_begin,
            .min_end    = exact_range.min_end,
            .local_m    = exact_range.local_m,
            .use_approx = false,
        };
    }

    const SInt min_begin = FindMinBoundaryByMass(config_.n, hyperedge_size, rank_, size_);

    const SInt min_end = FindMinBoundaryByMass(config_.n, hyperedge_size, rank_ + 1, size_);

    const SInt local_m = ApproximateLocalHyperedgeCount(hyperedge_size, m_k);

    return {
        .min_begin  = min_begin,
        .min_end    = min_end,
        .local_m    = local_m,
        .use_approx = true,
    };
}

template <typename BigInt>
void HyperGNM<BigInt>::GenerateHyperedgesOfSize(
    const SInt hyperedge_size, const SInt m_k, const SInt partition_id, const HGNMLocalGenerationRange& data) {
    if (hyperedge_size == 0 || hyperedge_size > config_.n || m_k == 0) {
        return;
    }

    if (!data.use_approx) {
        ExactFixedCountHyperedgeGenerator<BigInt> exact_generator(
            config_, partition_id, config_.k, rng_, mersenne_, graph_, memory_stats_, floyd_scratch_,
            debug_logger_ ? &*debug_logger_ : nullptr, &next_debug_hyperedge_id_
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
            ,
            &instrumentation_
#endif
        );

        exact_generator.Generate(
            hyperedge_size, FixedCountLocalRange{
                                .min_begin = data.min_begin,
                                .min_end   = data.min_end,
                                .local_m   = data.local_m,
                            });

        return;
    }

    const std::size_t cache_size = ComputeLogBinomCacheSize(data.local_m, data.min_begin, data.min_end);

    LogBinomCache log_binom_cache(hyperedge_size, cache_size);
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    instrumentation_.max_cache_initial_reserve =
        std::max<std::uint64_t>(instrumentation_.max_cache_initial_reserve, cache_size);
#endif
    log_binom_cache.InitializeRange(config_.n, data.min_begin, data.min_end);

    SInt edge_seed = EdgeSeed(hyperedge_size);

    mersenne_.RandomInit(edge_seed);

    GenerateSampledHyperedges(hyperedge_size, data.min_begin, data.min_end, data.local_m, edge_seed, log_binom_cache);

    AccumulateCacheStats(log_binom_cache, cache_size);
}

template <typename BigInt>
SInt HyperGNM<BigInt>::ApproximateLocalHyperedgeCount(const SInt hyperedge_size, const SInt m_k) {
    if (size_ == 1) {
        return m_k;
    }

    constexpr std::size_t cache_size = 4096;

    LogBinomCache log_binom_cache(hyperedge_size, cache_size);

    const SInt result =
        ApproximateLocalHyperedgeCountRecursive(hyperedge_size, 0, size_, rank_, 1.0L, m_k, 0, log_binom_cache);

    AccumulateCacheStats(log_binom_cache, cache_size);

    return result;
    return result;
}

template <typename BigInt>
SInt HyperGNM<BigInt>::RankSplitSeed(SInt hyperedge_size, SInt rank_begin, SInt rank_end, SInt level) const {
    return sampling::Spooky::hash(
        (static_cast<unsigned long long>(config_.seed) * kGNMCountSeedMultiplier)
        + (static_cast<unsigned long long>(hyperedge_size) * kCountSeedMultiplier)
        + (static_cast<unsigned long long>(rank_begin) * kRankSeedMultiplier)
        + (static_cast<unsigned long long>(rank_end) * kEdgeRankSeedMultiplier)
        + static_cast<unsigned long long>(level));
}

template <typename BigInt>
void HyperGNM<BigInt>::AccumulateCacheStats(const LogBinomCache& cache, const std::size_t requested_capacity) {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    const LogBinomCacheStats& stats = cache.stats;

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

    instrumentation_.cache_candidate_distance_observations += stats.candidate_distance_observations;

    instrumentation_.cache_candidate_max_distance =
        std::max(instrumentation_.cache_candidate_max_distance, stats.candidate_max_distance);

    instrumentation_.cache_backward_steps += stats.backward_steps;

    instrumentation_.cache_forward_steps += stats.forward_steps;

    instrumentation_.cache_max_size = std::max(instrumentation_.cache_max_size, stats.max_size);

    instrumentation_.max_cache_initial_reserve =
        std::max<std::uint64_t>(instrumentation_.max_cache_initial_reserve, requested_capacity);

    if (stats.map_calls + stats.candidate_calls > 0) {
        instrumentation_.cache_min_key = std::min(instrumentation_.cache_min_key, stats.min_key);

        instrumentation_.cache_max_key = std::max(instrumentation_.cache_max_key, stats.max_key);
    }
#endif
}

template <typename BigInt>
SInt HyperGNM<BigInt>::ApproximateLocalHyperedgeCountRecursive(
    const SInt hyperedge_size, const SInt rank_begin, const SInt rank_end, const SInt target_rank,
    const long double population_mass, const SInt draws, const SInt level, LogBinomCache& log_binom_cache) {
    if (draws == 0 || population_mass <= 0.0L) {
        return 0;
    }

    if (rank_end - rank_begin == 1) {
        return draws;
    }

    const SInt rank_mid = rank_begin + ((rank_end - rank_begin) / 2);

    const SInt left_min_begin = FindMinBoundaryByMass(config_.n, hyperedge_size, rank_begin, size_);
    const SInt left_min_end   = FindMinBoundaryByMass(config_.n, hyperedge_size, rank_mid, size_);

    const long double left_mass =
        MinRangeMassApproxCached(left_min_begin, left_min_end, config_.n, hyperedge_size, log_binom_cache);

    if (left_mass <= 0.0L) {
        if (target_rank < rank_mid) {
            return 0;
        }

        return ApproximateLocalHyperedgeCountRecursive(
            hyperedge_size, rank_mid, rank_end, target_rank, population_mass, draws, level + 1, log_binom_cache);
    }

    if (left_mass >= population_mass) {
        if (target_rank >= rank_mid) {
            return 0;
        }

        return ApproximateLocalHyperedgeCountRecursive(
            hyperedge_size, rank_begin, rank_mid, target_rank, population_mass, draws, level + 1, log_binom_cache);
    }

    const long double p = std::clamp(left_mass / population_mass, 0.0L, 1.0L);

    const SInt seed = RankSplitSeed(hyperedge_size, rank_begin, rank_end, level);

    const SInt left_draws = rng_.GenerateBinomial(seed, draws, static_cast<double>(p));

    if (target_rank < rank_mid) {
        return ApproximateLocalHyperedgeCountRecursive(
            hyperedge_size, rank_begin, rank_mid, target_rank, left_mass, left_draws, level + 1, log_binom_cache);
    }

    return ApproximateLocalHyperedgeCountRecursive(
        hyperedge_size, rank_mid, rank_end, target_rank, population_mass - left_mass, draws - left_draws, level + 1,
        log_binom_cache);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
template class HyperGNM<__int128>;
template class HyperGNM<SInt>;
#pragma GCC diagnostic pop

} // namespace kagen
