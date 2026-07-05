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
    const SInt vertices_per_pe = config_.n / size_;
    const SInt remainder       = config_.n % size_;

    const SInt start = (rank_ * vertices_per_pe) + std::min<SInt>(rank_, remainder);
    const SInt end   = start + vertices_per_pe + static_cast<SInt>(static_cast<SInt>(rank_) < remainder);

    return {start, end};
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
void HyperGNM<BigInt>::GenerateHyperedgesFromSizeCounts(const std::unordered_map<SInt, SInt>& size_counts) {
    for (const auto& [k, m_k]: size_counts) {
        if (debug_logger_) {
            const char* size_mode = config_.size_dist_pin_budget > 0 ? "boltzmann_pin_budget" : "implicit_min_owner";
            debug_logger_->LogSize(rank_, k, m_k, config_.m, size_mode);
        }

        GenerateHyperedgesOfSize(k, m_k);
    }
}

template <typename BigInt>
void HyperGNM<BigInt>::GenerateCSR() {
    graph_.hyperedge_offsets.push_back(0);

    auto [lower_bound, hard_upper_bound] = ResolveSizeRange();
    double alpha                         = ValidateAndGetSizeAlpha();

    std::unordered_map<SInt, SInt> size_counts;

    GenerateSizeCountsForConfig(lower_bound, hard_upper_bound, alpha, size_counts);

    LogSizeCounts(size_counts);

    GenerateHyperedgesFromSizeCounts(size_counts);

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
std::unordered_set<std::vector<SInt>, VectorHash> HyperGNM<BigInt>::MakeLocalSeenSet(const SInt local_m) const {
    std::unordered_set<std::vector<SInt>, VectorHash> local_seen;

    if (!config_.allow_duplicates) {
        local_seen.max_load_factor(0.5);
        local_seen.reserve(static_cast<std::size_t>(local_m));
    }

    return local_seen;
}

template <typename BigInt>
void HyperGNM<BigInt>::ValidateDuplicateCheckingIsFeasible(const SInt local_m) const {
    if (!config_.allow_duplicates && local_m > (SInt{1} << 26)) {
        throw ConfigurationError("Duplicate checking for huge hypergraph generation is infeasible; use --fast");
    }
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
SInt HyperGNM<BigInt>::ComputeLocalHyperedgeCount(const SInt hyperedge_size, const SInt m_k, const bool use_approx) {
    return use_approx ? ApproximateLocalHyperedgeCount(hyperedge_size, m_k)
                      : ExactLocalHyperedgeCount(hyperedge_size, m_k);
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
std::vector<SInt> HyperGNM<BigInt>::SampleHyperedge(SInt minimum_vertex, SInt hyperedge_size) {
    auto pins = FloydSample(minimum_vertex + 1, config_.n - minimum_vertex - 1, hyperedge_size - 1, mersenne_);

    pins.push_back(minimum_vertex);
    std::sort(pins.begin(), pins.end());

    return pins;
}

template <typename BigInt>
bool HyperGNM<BigInt>::TryPushHyperedge(
    const std::vector<SInt>& pins, std::unordered_set<std::vector<SInt>, VectorHash>& local_seen) {
    if (config_.allow_duplicates || local_seen.insert(pins).second) {
        PushHyperedge(pins);
        return true;
    }

    return false;
}

template <typename BigInt>
void HyperGNM<BigInt>::GenerateSampledHyperedges(
    SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, SInt local_m, SInt& edge_seed,
    LogBinomCache& log_binom_cache) {
    auto local_seen = MakeLocalSeenSet(local_m);

    const CountInt max_attempts = std::max(CountInt(local_m) * 10, CountInt(1000));

    CountInt attempts  = 0;
    SInt     generated = 0;

    while (generated < local_m) {
        const SInt minimum_vertex = SampleMinimumImplicit(
            local_min_begin, local_min_end, config_.n, hyperedge_size, rng_, edge_seed, log_binom_cache);

        auto pins = SampleHyperedge(minimum_vertex, hyperedge_size);

        if (TryPushHyperedge(pins, local_seen)) {
            ++generated;
        }

        if (++attempts > max_attempts) {
            throw ConfigurationError("HGNM rejection sampling exceeded attempt limit");
        }
    }
}

template <typename BigInt>
HGNMLocalGenerationRange HyperGNM<BigInt>::PrepareLocalGenerationRange(SInt hyperedge_size, SInt m_k) {
    const bool use_approx = config_.approx || config_.size_dist_pin_budget > 0;

    SInt min_begin = FindMinBoundaryByMass(config_.n, hyperedge_size, rank_, size_);
    SInt min_end   = FindMinBoundaryByMass(config_.n, hyperedge_size, rank_ + 1, size_);
    SInt local_m   = ComputeLocalHyperedgeCount(hyperedge_size, m_k, use_approx);

    ValidateDuplicateCheckingIsFeasible(local_m);

    if (!use_approx) {
        ValidateExactSparseDensity(hyperedge_size, min_begin, min_end, local_m);
    }

    return {.min_begin = min_begin, .min_end = min_end, .local_m = local_m, .use_approx = use_approx};
}

template <typename BigInt>
void HyperGNM<BigInt>::GenerateHyperedgesOfSize(const SInt hyperedge_size, const SInt m_k) {
    if (hyperedge_size > config_.n || m_k == 0) {
        return;
    }

    auto data = PrepareLocalGenerationRange(hyperedge_size, m_k);

    LogBinomCache log_binom_cache(hyperedge_size, ComputeCacheSize(data.local_m, data.min_begin, data.min_end));

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

template <typename BigInt>
SInt HyperGNM<BigInt>::ExactLocalHyperedgeCount(const SInt hyperedge_size, const SInt m_k) {
    const CountInt total_population = BinomialExact(config_.n, hyperedge_size);

    if (CountInt(m_k) > total_population) {
        throw ConfigurationError("Exact HGNM m_k exceeds hyperedge population");
    }

    return ExactLocalHyperedgeCountRecursive(hyperedge_size, 0, size_, rank_, total_population, m_k, 0);
}

template <typename BigInt>
SInt HyperGNM<BigInt>::ExactLocalHyperedgeCountRecursive(
    const SInt hyperedge_size, const SInt rank_begin, const SInt rank_end, const SInt target_rank,
    const CountInt population, const SInt draws, const SInt level) {
    if (draws == 0 || population == 0) {
        return 0;
    }

    if (rank_end - rank_begin == 1) {
        return draws;
    }

    const SInt rank_mid = rank_begin + ((rank_end - rank_begin) / 2);

    const SInt left_min_begin = FindMinBoundaryByMass(config_.n, hyperedge_size, rank_begin, size_);
    const SInt left_min_end   = FindMinBoundaryByMass(config_.n, hyperedge_size, rank_mid, size_);

    const CountInt left_population = MinRangeMassExact(left_min_begin, left_min_end, config_.n, hyperedge_size);

    if (left_population == 0) {
        if (target_rank < rank_mid) {
            return 0;
        }

        return ExactLocalHyperedgeCountRecursive(
            hyperedge_size, rank_mid, rank_end, target_rank, population, draws, level + 1);
    }

    if (left_population == population) {
        if (target_rank >= rank_mid) {
            return 0;
        }

        return ExactLocalHyperedgeCountRecursive(
            hyperedge_size, rank_begin, rank_mid, target_rank, population, draws, level + 1);
    }

    SInt seed = RankSplitSeed(hyperedge_size, rank_begin, rank_end, level);

    SInt left_draws = 0;
    if (population <= std::numeric_limits<SInt>::max() && left_population <= std::numeric_limits<SInt>::max()) {
        left_draws =
            rng_.GenerateHypergeometric(seed, left_population.convert_to<SInt>(), draws, population.convert_to<SInt>());
    } else {
        left_draws = HypergeometricCountSequential(population, left_population, draws, rng_, seed);
    }

    if (target_rank < rank_mid) {
        return ExactLocalHyperedgeCountRecursive(
            hyperedge_size, rank_begin, rank_mid, target_rank, left_population, left_draws, level + 1);
    }

    return ExactLocalHyperedgeCountRecursive(
        hyperedge_size, rank_mid, rank_end, target_rank, population - left_population, draws - left_draws, level + 1);
}

template <typename BigInt>
void HyperGNM<BigInt>::PushHyperedge(const std::vector<SInt>& pins) {
    graph_.hyperedge_pins.insert(graph_.hyperedge_pins.end(), pins.begin(), pins.end());
    graph_.hyperedge_offsets.push_back(static_cast<SInt>(graph_.hyperedge_pins.size()));
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
template class HyperGNM<__int128>;
template class HyperGNM<SInt>;
#pragma GCC diagnostic pop

} // namespace kagen
