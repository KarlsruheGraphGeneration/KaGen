#include "kagen/generators/hyper/h_erdos/hyper_gnp.h"

#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"

#include <algorithm>
#include <chrono>
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
            config_.output_graph.filename + "_hgnp_debug_rank_" + std::to_string(rank_) + ".csv", false);
    }
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateCSR() {
    graph_.hyperedge_offsets.push_back(0);

    const SInt lower_bound = config_.size_dist_lower_bound;
    const SInt hard_upper_bound =
        config_.size_dist_upper_bound != 0 ? std::min(config_.size_dist_upper_bound, config_.n) : config_.n;

    LoadSizeDistributionInputs();
    SelectProbabilityMode();

    SInt upper_bound = EffectiveUpperBound(hard_upper_bound, lower_bound);

    GenerateHyperedgesForAllSizes(lower_bound, upper_bound);

    SetLocalVertexRange();
}

template <typename BigInt>
void HyperGNP<BigInt>::FinalizeCSR(MPI_Comm /*comm*/) {}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateHyperedgesForAllSizes(const SInt lower_bound, const SInt upper_bound) {
    for (SInt hyperedge_size = lower_bound; hyperedge_size <= upper_bound; ++hyperedge_size) {
        if (!GenerateHyperedgesForSizeIfNeeded(hyperedge_size, lower_bound)) {
            break;
        }
    }
}

template <typename BigInt>
bool HyperGNP<BigInt>::GenerateHyperedgesForSizeIfNeeded(const SInt hyperedge_size, const SInt lower_bound) {
    const auto params = GetSizeGenerationParameters(hyperedge_size, lower_bound);

    if (ExpectedCountSkipsRemainingSizes(params)) {
        return false;
    }

    if (ExpectedCountSkipsCurrentSize(params)) {
        return true;
    }

    double probability = params.probability;

    if (params.expected_count > 0.0) {
        const auto resolved_probability = ResolveProbabilityForSize(hyperedge_size, params);

        if (!resolved_probability) {
            return probs_type_ != ProbabilityMode::EdgeBudget;
        }

        probability = *resolved_probability;
    }

    ValidateProbability(probability);

    if (probability <= 0.0) {
        return true;
    }

    LogSizeProbability(hyperedge_size, probability);
    GenerateHyperedgesOfSize(hyperedge_size, probability);

    return true;
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateHyperedgesOfSize(const SInt hyperedge_size, const double probability) {
    if (ShouldUseApproxGeneration(hyperedge_size)) {
        GenerateHyperedgesOfSizeApprox(hyperedge_size, probability);
    } else {
        GenerateHyperedgesOfSizeExact(hyperedge_size, probability);
    }
}

// #### Size /Probability Setup ####

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
    const bool explicit_probs    = !config_.size_probabilities.empty();
    const bool explicit_expected = !config_.size_expected_counts.empty();
    const bool budget_mode       = config_.edge_budget > 0.0;
    if (explicit_probs) {
        probs_type_ = ProbabilityMode::ExplicitProbabilities;
        return;
    }
    if (explicit_expected) {
        probs_type_ = ProbabilityMode::ExplicitExpectedCounts;
        return;
    }
    if (budget_mode) {
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
    return probs_type_ == ProbabilityMode::ExplicitExpectedCounts || probs_type_ == ProbabilityMode::EdgeBudget;
}

template <typename BigInt>
bool HyperGNP<BigInt>::ExpectedCountSkipsRemainingSizes(const SizeGenerationParameters& params) const {
    return probs_type_ == ProbabilityMode::EdgeBudget && ExpectedCountSkipsCurrentSize(params);
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

template <typename BigInt>
bool HyperGNP<BigInt>::ShouldUseApproxGeneration(const SInt hyperedge_size) const {
    return config_.approx || config_.n > (SInt{1} << 31) || hyperedge_size > 64;
}

// #### Partition helpers ####
template <typename BigInt>
void HyperGNP<BigInt>::SetLocalVertexRange() {
    const SInt vertices_per_pe = config_.n / size_;
    const SInt remainder       = config_.n % size_;

    const SInt start = (rank_ * vertices_per_pe) + std::min<SInt>(rank_, remainder);
    const SInt end   = start + vertices_per_pe + static_cast<SInt>(rank_ < remainder);

    SetVertexRange(start, end);
}

template <typename BigInt>
std::pair<SInt, SInt> HyperGNP<BigInt>::LocalMinOwnerRange(const SInt hyperedge_size) const {
    return {
        FindMinBoundaryByMass(config_.n, hyperedge_size, rank_, size_),
        FindMinBoundaryByMass(config_.n, hyperedge_size, rank_ + 1, size_)};
}

template <typename BigInt>
bool HyperGNP<BigInt>::ShouldSkipLocalRange(
    const SInt local_min_begin, const SInt local_min_end, const double probability) const {
    return local_min_begin >= local_min_end || probability <= 0.0;
}

// #### Approx generation ####
template <typename BigInt>
void HyperGNP<BigInt>::GenerateHyperedgesOfSizeApprox(const SInt hyperedge_size, const double probability) {
    if (ShouldSkipSizeGeneration(hyperedge_size, probability)) {
        return;
    }

    const auto [local_min_begin, local_min_end] = LocalMinOwnerRange(hyperedge_size);

    if (local_min_begin == local_min_end) {
        return;
    }

    LogBinomCache count_cache(hyperedge_size);
    GenerateApproxLocalRange(hyperedge_size, local_min_begin, local_min_end, probability, 0, count_cache);
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateApproxLocalRange(
    const SInt hyperedge_size, const SInt local_min_begin, const SInt local_min_end, const double probability,
    const SInt level, LogBinomCache& count_cache) {
    if (ShouldSkipLocalRange(local_min_begin, local_min_end, probability)) {
        return;
    }

    const auto stats =
        ComputeApproxRangeStats(hyperedge_size, local_min_begin, local_min_end, probability, count_cache);

    if (!stats) {
        return;
    }

    if (ShouldSplitApproxRange(local_min_begin, local_min_end, stats->expected)) {
        const SInt split =
            ChooseApproxRangeSplit(hyperedge_size, local_min_begin, local_min_end, stats->local_mass, count_cache);

        GenerateApproxLocalRange(hyperedge_size, local_min_begin, split, probability, level + 1, count_cache);
        GenerateApproxLocalRange(hyperedge_size, split, local_min_end, probability, level + 1, count_cache);

        return;
    }
    GenerateApproxLeafRange(hyperedge_size, local_min_begin, local_min_end, probability, level, *stats);
}

template <typename BigInt>
std::optional<ApproxRangeStats> HyperGNP<BigInt>::ComputeApproxRangeStats(
    const SInt hyperedge_size, const SInt local_min_begin, const SInt local_min_end, const double probability,
    LogBinomCache& count_cache) const {
    ApproxRangeStats stats;

    stats.local_mass = MinRangeMassApproxCached(local_min_begin, local_min_end, config_.n, hyperedge_size, count_cache);

    if (stats.local_mass <= 0.0L) {
        return std::nullopt;
    }

    stats.log_total = count_cache.Get(config_.n, hyperedge_size);

    const long double log_expected =
        std::log(static_cast<long double>(probability)) + std::log(stats.local_mass) + stats.log_total;

    stats.expected = log_expected > std::log(static_cast<long double>(std::numeric_limits<SInt>::max()))
                         ? static_cast<long double>(std::numeric_limits<SInt>::max())
                         : expl(log_expected);

    return stats;
}

template <typename BigInt>
bool HyperGNP<BigInt>::ShouldSplitApproxRange(
    const SInt local_min_begin, const SInt local_min_end, const long double expected) const {
    constexpr SInt kTargetEdgesPerBlock = 16384;
    constexpr SInt kTargetMinRangeWidth = 1ULL << 20;

    const SInt min_range_width = local_min_end - local_min_begin;

    return min_range_width > 1
           && (expected > static_cast<long double>(kTargetEdgesPerBlock) || min_range_width > kTargetMinRangeWidth);
}

template <typename BigInt>
SInt HyperGNP<BigInt>::ChooseApproxRangeSplit(
    const SInt hyperedge_size, const SInt local_min_begin, const SInt local_min_end, const long double local_mass,
    LogBinomCache& count_cache) {
    constexpr SInt kTargetMinRangeWidth = 1ULL << 20;

    const SInt min_range_width = local_min_end - local_min_begin;

    if (min_range_width > kTargetMinRangeWidth) {
        return local_min_begin + (min_range_width / 2);
    }

    return FindApproxMassSplit(local_min_begin, local_min_end, config_.n, hyperedge_size, local_mass, count_cache);
}

template <typename BigInt>
SInt HyperGNP<BigInt>::FindApproxMassSplit(
    const SInt begin, const SInt end, const SInt n, const SInt k, const long double total_mass, LogBinomCache& cache) {
    const long double half = total_mass * 0.5L;

    SInt left  = begin + 1;
    SInt right = end - 1;

    SInt best = begin + ((end - begin) / 2);

    while (left <= right) {
        const SInt mid = left + ((right - left) / 2);

        const long double left_mass = MinRangeMassApproxCached(begin, mid, n, k, cache);

        if (left_mass <= half) {
            best = mid;
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }

    return std::clamp<SInt>(best, begin + 1, end - 1);
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateApproxLeafRange(
    const SInt hyperedge_size, const SInt local_min_begin, const SInt local_min_end, const double probability,
    const SInt level, const ApproxRangeStats& stats) {
    const SInt count_seed = CalculateCountSeed(hyperedge_size, local_min_begin, level);

    const SInt local_edge_count =
        PoissonLocalCountFromScaledMass(stats.local_mass, stats.log_total, probability, rng_, count_seed);

    if (local_edge_count == 0) {
        return;
    }

    const SInt        search_width        = local_min_end - local_min_begin;
    const std::size_t expected_cache_size = static_cast<std::size_t>(
        std::min<SInt>(config_.n, std::max<SInt>(4096, std::min<SInt>(search_width, local_edge_count * 8))));

    LogBinomCache log_binom_cache(hyperedge_size, expected_cache_size);
    GenerateLocalHyperedges(hyperedge_size, local_min_begin, local_min_end, local_edge_count, log_binom_cache);
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateLocalHyperedges(
    const SInt hyperedge_size, const SInt local_min_begin, const SInt local_min_end, const SInt local_edge_count,
    LogBinomCache& log_binom_cache) {
    if (local_edge_count == 0) {
        return;
    }

    const auto time_begin = std::chrono::steady_clock::now();

    ValidateDuplicateCheckingIsFeasible(local_edge_count);

    auto local_seen = MakeLocalSeenSet(local_edge_count);

    SInt generated = 0;

    SInt edge_seed = LocalEdgeSeed(hyperedge_size);

    const CountInt max_attempts = std::max(CountInt(local_edge_count) * 10, CountInt(1000));
    CountInt       attempts     = 0;
    mersenne_.RandomInit(edge_seed);

    while (generated < local_edge_count) {
        auto pins = SampleLocalHyperedge(hyperedge_size, local_min_begin, local_min_end, edge_seed, log_binom_cache);

        if (config_.allow_duplicates || local_seen.insert(pins).second) {
            PushHyperedge(pins);
            ++generated;
        }

        if (++attempts > max_attempts) {
            throw ConfigurationError("HGNP rejection sampling exceeded attempt limit");
        }
    }

    const auto time_end = std::chrono::steady_clock::now();

    const double seconds = std::chrono::duration<double>(time_end - time_begin).count();

    const long double num_pins = static_cast<long double>(local_edge_count) * static_cast<long double>(hyperedge_size);

    const double us_per_pin = num_pins > 0 ? (seconds * 1e6) / static_cast<double>(num_pins) : 0.0;
    if (config_.debug) {
        std::cout << "HGNP rank=" << rank_ << " k=" << hyperedge_size << " local_m=" << local_edge_count
                  << " pins=" << num_pins << " attempts=" << attempts << " us_per_pin=" << us_per_pin
                  << " min_begin=" << local_min_begin << " min_end=" << local_min_end << '\n';
    }
}

template <typename BigInt>
std::vector<SInt> HyperGNP<BigInt>::SampleLocalHyperedge(
    const SInt hyperedge_size, const SInt local_min_begin, const SInt local_min_end, SInt& edge_seed,
    LogBinomCache& log_binom_cache) {
    const SInt minimum_vertex = SampleMinimumImplicit(
        local_min_begin, local_min_end, config_.n, hyperedge_size, rng_, edge_seed, log_binom_cache);

    auto pins = FloydSample(minimum_vertex + 1, config_.n - minimum_vertex - 1, hyperedge_size - 1, mersenne_);

    pins.insert(pins.begin(), minimum_vertex);

    return pins;
}

template <typename BigInt>
SInt HyperGNP<BigInt>::LocalEdgeSeed(const SInt hyperedge_size) const {
    return sampling::Spooky::hash(
        static_cast<unsigned long long>(config_.seed)
        + (static_cast<unsigned long long>(rank_) * kEdgeRankSeedMultiplier)
        + (static_cast<unsigned long long>(hyperedge_size) * kEdgeSeedMultiplier));
}

// #### Exact generation
template <typename BigInt>
void HyperGNP<BigInt>::GenerateHyperedgesOfSizeExact(const SInt hyperedge_size, const double probability) {
    if (ShouldSkipSizeGeneration(hyperedge_size, probability)) {
        return;
    }

    const auto [local_min_begin, local_min_end] = LocalMinOwnerRange(hyperedge_size);

    if (local_min_begin == local_min_end) {
        return;
    }

    ValidateExactMinOwnerPartition(hyperedge_size);

    GenerateExactOwnedPrefixRange(hyperedge_size, local_min_begin, local_min_end, probability, 0);
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateExactOwnedPrefixRange(
    const SInt hyperedge_size, const SInt local_min_begin, const SInt local_min_end, const double probability,
    const SInt level) {
    GenerateExactPrefixRange(hyperedge_size, {}, local_min_begin, local_min_end, hyperedge_size, probability, level);
}

// Prefix recursion invariant:
//
// prefix contains vertices fixed so far.
// remaining_pins is the number of vertices still to choose.
// The next selected vertex must lie in [candidate_begin, candidate_end).
//
// This subtree represents exactly all hyperedges extending prefix whose
// next selected vertex lies in that interval.
template <typename BigInt>
void HyperGNP<BigInt>::GenerateExactPrefixRange(
    const SInt hyperedge_size, std::vector<SInt> prefix, const SInt candidate_begin, const SInt candidate_end,
    const SInt remaining_pins, const double probability, const SInt level) {
    if (remaining_pins == 0) {
        PushHyperedge(prefix);
        return;
    }

    const CountInt population = PrefixRangePopulation(candidate_begin, candidate_end, remaining_pins);

    if (population == 0) {
        return;
    }

    if (population <= std::numeric_limits<SInt>::max()) {
        GenerateExactPrefixLeafRange(
            hyperedge_size, std::move(prefix), candidate_begin, candidate_end, remaining_pins, probability, level,
            population);
        return;
    }

    const SInt split = FindExactPrefixPopulationSplit(candidate_begin, candidate_end, remaining_pins, population);

    GenerateExactPrefixRange(hyperedge_size, prefix, candidate_begin, split, remaining_pins, probability, level + 1);

    GenerateExactPrefixRange(
        hyperedge_size, std::move(prefix), split, candidate_end, remaining_pins, probability, level + 1);
}

// Counts all completions represented by this prefix range.
// If the next selected vertex is v, the remaining pins can be chosen from
// vertices greater than v in C(n - v - 1, remaining_pins - 1) ways.
template <typename BigInt>
CountInt HyperGNP<BigInt>::PrefixRangePopulation(
    const SInt candidate_begin, const SInt candidate_end, const SInt remaining_pins) const {
    CountInt population = 0;

    for (SInt v = candidate_begin; v < candidate_end; ++v) {
        population += BinomialExact(config_.n - v - 1, remaining_pins - 1);
    }

    return population;
}

template <typename BigInt>
SInt HyperGNP<BigInt>::FindExactPrefixPopulationSplit(
    const SInt begin, const SInt end, const SInt remaining_pins, const CountInt& total_population) const {
    const CountInt half = total_population / 2;

    SInt left  = begin + 1;
    SInt right = end - 1;
    SInt best  = begin + ((end - begin) / 2);

    while (left <= right) {
        const SInt mid = left + ((right - left) / 2);

        const CountInt left_population = PrefixRangePopulation(begin, mid, remaining_pins);

        if (left_population <= half) {
            best = mid;
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }

    return std::clamp<SInt>(best, begin + 1, end - 1);
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateExactPrefixLeafRange(
    const SInt hyperedge_size, std::vector<SInt> prefix, const SInt candidate_begin, const SInt candidate_end,
    const SInt remaining_pins, const double probability, const SInt level, const CountInt& population) {
    const SInt count_seed = CalculateCountSeed(hyperedge_size, candidate_begin, level);

    const SInt local_population = population.convert_to<SInt>();
    const SInt local_m          = rng_.GenerateBinomial(count_seed, local_population, probability);

    if (local_m == 0) {
        return;
    }

    ValidateDuplicateCheckingIsFeasible(local_m);
    ValidateExactSparseDensity(local_m, population);

    Mersenne prefix_mersenne;
    prefix_mersenne.RandomInit(count_seed);

    GenerateSampledPrefixHyperedges(prefix, candidate_begin, candidate_end, remaining_pins, local_m, prefix_mersenne);
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateSampledPrefixHyperedges(
    const std::vector<SInt>& prefix, SInt candidate_begin, SInt candidate_end, SInt remaining_pins, SInt local_m,
    Mersenne& mersenne) {
    auto local_seen = MakeLocalSeenSet(local_m);

    SInt           generated    = 0;
    CountInt       attempts     = 0;
    const CountInt max_attempts = std::max(CountInt(local_m) * 10, CountInt(1000));

    while (generated < local_m) {
        auto pins = SamplePrefixHyperedge(prefix, candidate_begin, candidate_end, remaining_pins, mersenne);

        if (config_.allow_duplicates || local_seen.insert(pins).second) {
            PushHyperedge(pins);
            ++generated;
        }

        if (++attempts > max_attempts) {
            throw ConfigurationError("HGNP exact prefix rejection sampling exceeded attempt limit");
        }
    }
}

template <typename BigInt>
std::vector<SInt> HyperGNP<BigInt>::SamplePrefixHyperedge(
    const std::vector<SInt>& prefix, SInt candidate_begin, SInt candidate_end, SInt remaining_pins,
    Mersenne& mersenne) {
    auto pins = prefix;

    const SInt next_vertex = SampleNextPrefixVertex(candidate_begin, candidate_end, remaining_pins, mersenne);

    pins.push_back(next_vertex);

    auto suffix = FloydSample(next_vertex + 1, config_.n - next_vertex - 1, remaining_pins - 1, mersenne);

    pins.insert(pins.end(), suffix.begin(), suffix.end());

    return pins;
}

// Sample the next prefix vertex proportional to the number of suffix
// completions it represents, not uniformly over vertices.
template <typename BigInt>
SInt HyperGNP<BigInt>::SampleNextPrefixVertex(
    const SInt candidate_begin, const SInt candidate_end, const SInt remaining_pins, Mersenne& mersenne) const {
    const CountInt total = PrefixRangePopulation(candidate_begin, candidate_end, remaining_pins);

    const long double x      = MersenneUniform01(mersenne);
    CountInt          target = static_cast<CountInt>(x * total.convert_to<long double>());

    CountInt prefix_sum = 0;

    for (SInt v = candidate_begin; v < candidate_end; ++v) {
        const CountInt weight = BinomialExact(config_.n - v - 1, remaining_pins - 1);

        if (target < prefix_sum + weight) {
            return v;
        }

        prefix_sum += weight;
    }

    return candidate_end - 1;
}

template <typename BigInt>
void HyperGNP<BigInt>::ValidateExactSparseDensity(const SInt local_m, const CountInt& population) const {
    const long double density = static_cast<long double>(local_m) / population.convert_to<long double>();

    if (density > 0.25L) {
        throw ConfigurationError("Dense exact HGNP not implemented yet");
    }
}

template <typename BigInt>
void HyperGNP<BigInt>::ValidateExactMinOwnerPartition(const SInt hyperedge_size) {
    if (!config_.debug || config_.n > 10000 || hyperedge_size > 8) {
        return;
    }

    CountInt check_sum = 0;

    for (SInt r = 0; r < size_; ++r) {
        const SInt begin = FindMinBoundaryByMass(config_.n, hyperedge_size, r, size_);
        const SInt end   = FindMinBoundaryByMass(config_.n, hyperedge_size, r + 1, size_);

        if (begin > end) {
            throw ConfigurationError("HGNP exact validation failed: invalid min range");
        }

        check_sum += MinRangeMassExact(begin, end, config_.n, hyperedge_size);
    }

    const CountInt total = BinomialExact(config_.n, hyperedge_size);

    if (check_sum != total) {
        throw ConfigurationError("HGNP exact validation failed: local populations do not sum to total population");
    }
}

// #### Sampling/duplication helpers ####
template <typename BigInt>
void HyperGNP<BigInt>::ValidateDuplicateCheckingIsFeasible(const SInt local_edge_count) const {
    if (!config_.allow_duplicates && local_edge_count > (SInt{1} << 26)) {
        throw ConfigurationError("Duplicate checking for huge hypergraph generation is infeasible; use --fast");
    }
}

template <typename BigInt>
std::unordered_set<std::vector<SInt>, VectorHash>
HyperGNP<BigInt>::MakeLocalSeenSet(const SInt local_edge_count) const {
    std::unordered_set<std::vector<SInt>, VectorHash> local_seen;

    if (!config_.allow_duplicates) {
        local_seen.max_load_factor(0.5);
        local_seen.reserve(static_cast<std::size_t>(local_edge_count));
    }

    return local_seen;
}

template <typename BigInt>
SInt HyperGNP<BigInt>::CalculateCountSeed(SInt hyperedge_size, SInt candidate_begin, SInt level) const {
    return sampling::Spooky::hash(
        static_cast<unsigned long long>(config_.seed)
        + (static_cast<unsigned long long>(hyperedge_size) * kCountSeedMultiplier)
        + (static_cast<unsigned long long>(rank_) * kRankSeedMultiplier)
        + (static_cast<unsigned long long>(candidate_begin) * kEdgeSeedMultiplier)
        + (static_cast<unsigned long long>(level) * kEdgeRankSeedMultiplier));
}

// #### Logging/graph mutation ####

template <typename BigInt>
void HyperGNP<BigInt>::LogSizeProbability(const SInt hyperedge_size, const double probability) {
    if (!debug_logger_) {
        return;
    }

    std::ostringstream p_info;
    p_info << std::scientific << probability;

    debug_logger_->LogBlock(
        rank_, hyperedge_size, "size", "min_owner", 0, 0, 0, "unknown", 0, 0, "0", 0, 0, "p=" + p_info.str());
}

template <typename BigInt>
void HyperGNP<BigInt>::PushHyperedge(const std::vector<SInt>& pins) {
    graph_.hyperedge_pins.insert(graph_.hyperedge_pins.end(), pins.begin(), pins.end());
    graph_.hyperedge_offsets.push_back(static_cast<SInt>(graph_.hyperedge_pins.size()));
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
template class HyperGNP<__int128>;
template class HyperGNP<SInt>;
#pragma GCC diagnostic pop
} // namespace kagen
