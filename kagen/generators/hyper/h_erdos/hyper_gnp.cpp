#include "kagen/generators/hyper/h_erdos/hyper_gnp.h"

#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"

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
            config_.output_graph.filename + "_hgnp_debug_rank_" + std::to_string(rank_) + ".csv", false);
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

    LoadSizeDistributionInputs();
    SelectProbabilityMode();

    const SInt upper_bound = EffectiveUpperBound(hard_upper_bound, lower_bound);

    if (probs_type_ == ProbabilityMode::EdgeAndPinBudget) {
        GenerateBudgetSizeCounts(lower_bound, hard_upper_bound);
    }

    const std::vector<HGNPSizePlan> plan = BuildGenerationPlan(lower_bound, upper_bound);

    ReserveCSRForPlan(plan);

    graph_.hyperedge_offsets.push_back(0);
    GenerateHyperedgesFromPlan(plan);

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
        config_, rank_, size_, rng_, mersenne_, graph_, memory_stats_, floyd_scratch_);

    fixed_count_generator.Generate(
        entry.hyperedge_size, FixedCountLocalRange{
                                  .min_begin = entry.range.begin,
                                  .min_end   = entry.range.end,
                                  .local_m   = entry.range.local_m,
                              });
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateApproxHyperedgesFromPlan(const HGNPSizePlan& entry) {
    const std::size_t cache_size = ComputeCacheSize(entry.range.local_m, entry.range.begin, entry.range.end);

    LogBinomCache cache(entry.hyperedge_size, cache_size);

    cache.InitializeRange(config_.n, entry.range.begin, entry.range.end);

    std::vector<SInt> pins;
    pins.reserve(static_cast<std::size_t>(entry.hyperedge_size));

    auto seen = MakeLocalSeenSet(config_.allow_duplicates, entry.range.local_m);

    SInt edge_seed = LocalEdgeSeed(entry.hyperedge_size);

    mersenne_.RandomInit(edge_seed);

    GenerateLocalHyperedges(entry, edge_seed, cache, pins, seen);
}

template <typename BigInt>
void HyperGNP<BigInt>::FinalizeCSR(MPI_Comm /*comm*/) {}

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

    const auto [min_begin, min_end] = LocalMinOwnerRange(k);

    entry.range.begin   = min_begin;
    entry.range.end     = min_end;
    entry.range.local_m = 0;

    if (min_begin >= min_end || probability <= 0.0) {
        return;
    }

    const CountInt local_population = MinRangeMassExact(min_begin, min_end, config_.n, k);

    if (local_population == 0) {
        return;
    }

    const SInt local_m = SampleExactEdgeCount(local_population, probability, LocalCountSeed(k));

    if (CountInt(local_m) > local_population) {
        throw ConfigurationError("HGNP local binomial count exceeds local population");
    }

    entry.range.local_m = local_m;
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

    if (upper_bound >= lower_bound) {
        plan.reserve(static_cast<std::size_t>(upper_bound - lower_bound + 1));
    }

    // To avoid overflow of k for upper_bound = SInt.max()
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

    LogSizeProbability(hyperedge_size, probability);

    HGNPSizePlan entry = PrepareSizePlan(hyperedge_size, probability);

    // Zero-count entries need not be retained.
    if (entry.range.local_m > 0) {
        plan.push_back(std::move(entry));
    }

    return true;
}

template <typename BigInt>
HGNPSizePlan HyperGNP<BigInt>::PrepareSizePlan(const SInt hyperedge_size, const double probability) {
    HGNPSizePlan entry;
    entry.hyperedge_size = hyperedge_size;

    if (config_.approx) {
        entry.range = PrepareApproxLocalRange(hyperedge_size, probability);
    } else {
        PrepareSampledExactPlan(entry, probability);
    }

    return entry;
}

template <typename BigInt>
SInt HyperGNP<BigInt>::LocalCountSeed(const SInt hyperedge_size) const {
    return sampling::Spooky::hash(
        static_cast<unsigned long long>(config_.seed)
        + (static_cast<unsigned long long>(hyperedge_size) * kCountSeedMultiplier)
        + (static_cast<unsigned long long>(rank_) * kRankSeedMultiplier));
}

template <typename BigInt>
HGNPLocalGenerationRange
HyperGNP<BigInt>::PrepareApproxLocalRange(const SInt hyperedge_size, const double probability) {
    const auto [begin, end] = LocalMinOwnerRange(hyperedge_size);

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

    range.local_m = lambda > 0.0 ? rng_.GeneratePoisson(LocalCountSeed(hyperedge_size), lambda) : 0;

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
std::size_t HyperGNP<BigInt>::ComputeCacheSize(const SInt local_m, const SInt begin, const SInt end) const {
    const SInt search_width = end - begin;

    const SInt scaled_m =
        local_m > std::numeric_limits<SInt>::max() / 8 ? std::numeric_limits<SInt>::max() : local_m * 8;

    return static_cast<std::size_t>(
        std::min<SInt>(config_.n, std::max<SInt>(4096, std::min<SInt>(search_width, scaled_m))));
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
std::pair<SInt, SInt> HyperGNP<BigInt>::LocalMinOwnerRange(const SInt hyperedge_size) const {
    return {
        FindMinBoundaryByMass(config_.n, hyperedge_size, rank_, size_),
        FindMinBoundaryByMass(config_.n, hyperedge_size, rank_ + 1, size_)};
}

// #### Approx generation ####

template <typename BigInt>
void HyperGNP<BigInt>::GenerateLocalHyperedges(
    const HGNPSizePlan& entry, SInt& edge_seed, LogBinomCache& cache, std::vector<SInt>& pins,
    std::unordered_set<std::uint64_t>& seen) {
    const SInt local_m = entry.range.local_m;

    const CountInt max_attempts = std::max(CountInt(local_m) * 10, CountInt(1000));

    CountInt attempts  = 0;
    SInt     generated = 0;

    while (generated < local_m) {
        SampleLocalHyperedgeInto(entry.hyperedge_size, entry.range.begin, entry.range.end, edge_seed, cache, pins);

        if (config_.allow_duplicates || seen.insert(FingerprintPins(pins)).second) {
            PushUncompressedHyperedge(graph_, memory_stats_, pins);

            ++generated;
        }

        if (++attempts > max_attempts) {
            throw ConfigurationError("HGNP sampling exceeded attempt limit");
        }
    }
}

template <typename BigInt>
void HyperGNP<BigInt>::SampleLocalHyperedgeInto(
    const SInt hyperedge_size, const SInt local_min_begin, const SInt local_min_end, SInt& edge_seed,
    LogBinomCache& log_binom_cache, std::vector<SInt>& pins) {
    pins.clear();

    const SInt minimum_vertex = SampleMinimumImplicit(
        local_min_begin, local_min_end, config_.n, hyperedge_size, rng_, edge_seed, log_binom_cache);
    pins.push_back(minimum_vertex);

    FloydSampleInto(
        minimum_vertex + 1, config_.n - minimum_vertex - 1, hyperedge_size - 1, mersenne_, pins, floyd_scratch_, 1);
}

template <typename BigInt>
SInt HyperGNP<BigInt>::LocalEdgeSeed(const SInt hyperedge_size) const {
    return sampling::Spooky::hash(
        static_cast<unsigned long long>(config_.seed)
        + (static_cast<unsigned long long>(rank_) * kEdgeRankSeedMultiplier)
        + (static_cast<unsigned long long>(hyperedge_size) * kEdgeSeedMultiplier));
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
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
template class HyperGNP<__int128>;
template class HyperGNP<SInt>;
#pragma GCC diagnostic pop
} // namespace kagen
