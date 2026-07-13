#include "kagen/generators/hyper/h_erdos/hyper_gnm.h"

#include "kagen/generators/generator.h"
#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"
#include "kagen/tools/rng_wrapper.h"

#include <algorithm>
#include <cmath>
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
            config_.output_graph.filename + "_hgnm_debug_rank_" + std::to_string(rank_) + ".csv", false);
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
        GenerateHyperedgesOfSize(entry.hyperedge_size, entry.m_k, entry.range);
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

    for (const auto& [k, m_k]: size_counts) {
        auto range = PrepareLocalGenerationRange(k, m_k);

        expected_local_edges += static_cast<std::size_t>(range.local_m);
        expected_local_pins += static_cast<std::size_t>(range.local_m) * static_cast<std::size_t>(k);

        plan.push_back({
            .hyperedge_size = k,
            .m_k            = m_k,
            .range          = range,
        });
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
void HyperGNM<BigInt>::FinalizeCSR(MPI_Comm /*comm*/) {}

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
std::size_t HyperGNM<BigInt>::ComputeCacheSize(SInt local_m, SInt local_min_begin, SInt local_min_end) const {
    const SInt search_width = local_min_end - local_min_begin;

    return static_cast<std::size_t>(
        std::min<SInt>(config_.n, std::max<SInt>(4096, std::min<SInt>(search_width, local_m * 8))));
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
    pins.resize(hyperedge_size);
    pins[0] = minimum_vertex;

    FloydSampleInto(
        minimum_vertex + 1, config_.n - minimum_vertex - 1, hyperedge_size - 1, mersenne_, pins, floyd_scratch_, 1);
}

template <typename BigInt>
bool HyperGNM<BigInt>::TryPushHyperedge(const std::vector<SInt>& pins, std::unordered_set<std::uint64_t>& local_seen) {
    if (config_.allow_duplicates || local_seen.insert(FingerprintPins(pins)).second) {
        PushUncompressedHyperedge(graph_, memory_stats_, pins);
        return true;
    }

    return false;
}

template <typename BigInt>
void HyperGNM<BigInt>::GenerateSampledHyperedges(
    SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, SInt local_m, SInt& edge_seed,
    LogBinomCache& log_binom_cache) {
    auto              local_seen = MakeLocalSeenSet(config_.allow_duplicates, local_m);
    std::vector<SInt> pins;
    pins.reserve(hyperedge_size);

    const CountInt max_attempts = std::max(CountInt(local_m) * 10, CountInt(1000));

    CountInt attempts  = 0;
    SInt     generated = 0;

    while (generated < local_m) {
        const SInt minimum_vertex = SampleMinimumImplicit(
            local_min_begin, local_min_end, config_.n, hyperedge_size, rng_, edge_seed, log_binom_cache);

        SampleHyperedgeInto(minimum_vertex, hyperedge_size, pins);

        if (TryPushHyperedge(pins, local_seen)) {
            ++generated;
        }

        if (++attempts > max_attempts) {
            throw ConfigurationError("HGNM rejection sampling exceeded attempt limit");
        }
    }
}

template <typename BigInt>
HGNMLocalGenerationRange HyperGNM<BigInt>::PrepareLocalGenerationRange(const SInt hyperedge_size, const SInt m_k) {
    if (!config_.approx) {
        ExactFixedCountHyperedgeGenerator<BigInt> exact_generator(
            config_, rank_, size_, rng_, mersenne_, graph_, memory_stats_, floyd_scratch_);

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
    const SInt hyperedge_size, const SInt m_k, const HGNMLocalGenerationRange& data) {
    if (hyperedge_size == 0 || hyperedge_size > config_.n || m_k == 0) {
        return;
    }

    if (!data.use_approx) {
        ExactFixedCountHyperedgeGenerator<BigInt> exact_generator(
            config_, rank_, size_, rng_, mersenne_, graph_, memory_stats_, floyd_scratch_);

        exact_generator.Generate(
            hyperedge_size, FixedCountLocalRange{
                                .min_begin = data.min_begin,
                                .min_end   = data.min_end,
                                .local_m   = data.local_m,
                            });

        return;
    }

    LogBinomCache log_binom_cache(hyperedge_size, ComputeCacheSize(data.local_m, data.min_begin, data.min_end));

    log_binom_cache.InitializeRange(config_.n, data.min_begin, data.min_end);

    SInt edge_seed = EdgeSeed(hyperedge_size);

    mersenne_.RandomInit(edge_seed);

    GenerateSampledHyperedges(hyperedge_size, data.min_begin, data.min_end, data.local_m, edge_seed, log_binom_cache);
}

template <typename BigInt>
SInt HyperGNM<BigInt>::ApproximateLocalHyperedgeCount(const SInt hyperedge_size, const SInt m_k) {
    if (size_ == 1) {
        return m_k;
    }

    LogBinomCache log_binom_cache(hyperedge_size);

    return ApproximateLocalHyperedgeCountRecursive(hyperedge_size, 0, size_, rank_, 1.0L, m_k, 0, log_binom_cache);
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
