#include "kagen/generators/hyper/h_erdos/hyper_gnp.h"

#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"
#include "kagen/sampling/hash.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

namespace kagen {

std::unique_ptr<Generator>
HyperGNPFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
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
std::pair<double, double> HyperGNP<BigInt>::CalculateProbability(SInt hyperedge_size, SInt lower_bound) {
    double p            = 0.0;
    double expected_m_k = -1.0;
    switch (probs_type_) {
        case EXPLICIT_PROBS:
            p = config_.size_probabilities[hyperedge_size - lower_bound];
            break;
        case EXPLICIT_EXPECTED:
            expected_m_k = config_.size_expected_counts[hyperedge_size - lower_bound];
            break;
        case BUDGET_MODE:
            expected_m_k = config_.edge_budget * config_.size_decay
                           * std::pow(1.0 - config_.size_decay, static_cast<double>(hyperedge_size - lower_bound));
            break;
        default:
            p = config_.p;
            break;
    }
    return {p, expected_m_k};
}

template <typename BigInt>
void HyperGNP<BigInt>::SetProbability() {
    const bool explicit_probs    = !config_.size_probabilities.empty();
    const bool explicit_expected = !config_.size_expected_counts.empty();
    const bool budget_mode       = config_.edge_budget > 0.0;
    if (explicit_probs) {
        probs_type_ = ProbabilityType::EXPLICIT_PROBS;
        return;
    }
    if (explicit_expected) {
        probs_type_ = ProbabilityType::EXPLICIT_EXPECTED;
        return;
    }
    if (budget_mode) {
        probs_type_ = ProbabilityType::BUDGET_MODE;
        return;
    }

    probs_type_ = ProbabilityType::GLOBAL_PROBABILITY;
}

template <typename BigInt>
SInt HyperGNP<BigInt>::SetUpperBound(SInt hard_upper_bound, SInt lower_bound) {
    switch (probs_type_) {
        case EXPLICIT_PROBS:
            return std::min<SInt>(lower_bound + static_cast<SInt>(config_.size_probabilities.size()) - 1, config_.n);
        case EXPLICIT_EXPECTED:
            return std::min<SInt>(lower_bound + static_cast<SInt>(config_.size_expected_counts.size()) - 1, config_.n);
        default:
            return hard_upper_bound;
    }
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateCSR() {
    graph_.hyperedge_offsets.push_back(0);

    const SInt lower_bound = config_.size_dist_lower_bound;
    const SInt hard_upper_bound =
        config_.size_dist_upper_bound != 0 ? std::min(config_.size_dist_upper_bound, config_.n) : config_.n;

    if (!config_.size_probabilities_file.empty()) {
        std::ifstream in(config_.size_probabilities_file);
        if (!in) {
            throw ConfigurationError("Could not open HGNP size probability file");
        }

        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }

            config_.size_probabilities.push_back(std::stod(line));
        }
    }

    if (!config_.size_expected_counts_file.empty()) {
        std::ifstream in(config_.size_expected_counts_file);
        if (!in) {
            throw ConfigurationError("Could not open HGNP size expected counts file");
        }

        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }

            config_.size_expected_counts.push_back(std::stod(line));
        }
    }

    SetProbability();

    SInt upper_bound = SetUpperBound(hard_upper_bound, lower_bound);

    for (SInt hyperedge_size = lower_bound; hyperedge_size <= upper_bound; ++hyperedge_size) {
        std::pair<double, double> probability_factors = CalculateProbability(hyperedge_size, lower_bound);
        double                    p                   = probability_factors.first;
        double                    expected_m_k        = probability_factors.second;

        if (expected_m_k == 0.0 || ExpectedCountIsNegligible(expected_m_k)) {
            if (probs_type_ == BUDGET_MODE) {
                break;
            }
            continue;
        }

        if (expected_m_k > 0.0) {
            const long double log_p =
                std::log(static_cast<long double>(expected_m_k)) - LogBinomialApprox(config_.n, hyperedge_size);

            if (log_p < std::log(static_cast<long double>(std::numeric_limits<double>::denorm_min()))) {
                if (probs_type_ == BUDGET_MODE) {
                    break;
                }
                continue;
            }

            p = std::clamp(static_cast<double>(expl(log_p)), 0.0, 1.0);
        }

        if (p < 0.0 || p > 1.0) {
            throw ConfigurationError("HGNP probability must be in [0, 1]");
        }

        if (p <= 0.0) {
            continue;
        }

        if (debug_logger_) {
            std::ostringstream p_info;
            p_info << std::scientific << p;

            debug_logger_->LogBlock(
                rank_, hyperedge_size, "size", "min_owner", 0, 0, 0, "unknown", 0, 0, "0", 0, 0, "p=" + p_info.str());
        }

        GenerateHyperedgesOfSize(hyperedge_size, p);
    }

    const SInt vertices_per_pe = config_.n / size_;
    const SInt remainder       = config_.n % size_;

    const SInt start = (rank_ * vertices_per_pe) + std::min<SInt>(rank_, remainder);
    const SInt end   = start + vertices_per_pe + static_cast<SInt>(static_cast<SInt>(rank_) < remainder);

    SetVertexRange(start, end);
}
template <typename BigInt>
void HyperGNP<BigInt>::GenerateHyperedgesOfSize(const SInt hyperedge_size, const double p) {
    if (config_.approx) {
        GenerateHyperedgesOfSizeApprox(hyperedge_size, p);
    } else {
        GenerateHyperedgesOfSizeExact(hyperedge_size, p);
    }
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateHyperedgesOfSizeApprox(const SInt hyperedge_size, const double p) {
    if (hyperedge_size == 0 || hyperedge_size > config_.n || p <= 0.0) {
        return;
    }

    if (p > 1.0) {
        throw ConfigurationError("HGNP probabilities must be in [0, 1]");
    }

    const SInt local_min_begin = FindMinBoundaryByMass(config_.n, hyperedge_size, rank_, size_);

    const SInt local_min_end = FindMinBoundaryByMass(config_.n, hyperedge_size, rank_ + 1, size_);

    if (local_min_begin == local_min_end) {
        return;
    }

    LogBinomCache count_cache(hyperedge_size);

    GenerateApproxLocalRange(hyperedge_size, local_min_begin, local_min_end, p, 0, count_cache);
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
void HyperGNP<BigInt>::GenerateApproxLocalRange(
    const SInt hyperedge_size, const SInt local_min_begin, const SInt local_min_end, const double p, const SInt level,
    LogBinomCache& count_cache) {
    constexpr SInt kTargetEdgesPerBlock = 16384;
    constexpr SInt kTargetMinRangeWidth = 1ULL << 20; // 1,048,576

    if (local_min_begin >= local_min_end || p <= 0.0) {
        return;
    }

    const long double local_mass =
        MinRangeMassApproxCached(local_min_begin, local_min_end, config_.n, hyperedge_size, count_cache);

    if (local_mass <= 0.0L) {
        return;
    }

    const long double log_total = count_cache.Get(config_.n, hyperedge_size);

    const long double log_expected = std::log(static_cast<long double>(p)) + std::log(local_mass) + log_total;

    const long double expected = log_expected > std::log(static_cast<long double>(std::numeric_limits<SInt>::max()))
                                     ? static_cast<long double>(std::numeric_limits<SInt>::max())
                                     : expl(log_expected);

    const SInt min_range_width = local_min_end - local_min_begin;

    const bool split_by_edges = expected > static_cast<long double>(kTargetEdgesPerBlock);

    const bool split_by_width = min_range_width > kTargetMinRangeWidth;

    if ((split_by_edges || split_by_width) && min_range_width > 1) {
        const SInt split =
            split_by_width ? local_min_begin + (min_range_width / 2)
                           : FindApproxMassSplit(
                                 local_min_begin, local_min_end, config_.n, hyperedge_size, local_mass, count_cache);

        GenerateApproxLocalRange(hyperedge_size, local_min_begin, split, p, level + 1, count_cache);

        GenerateApproxLocalRange(hyperedge_size, split, local_min_end, p, level + 1, count_cache);

        return;
    }

    const SInt count_seed = sampling::Spooky::hash(
        static_cast<unsigned long long>(config_.seed)
        + (static_cast<unsigned long long>(hyperedge_size) * kCountSeedMultiplier)
        + (static_cast<unsigned long long>(rank_) * kRankSeedMultiplier)
        + (static_cast<unsigned long long>(local_min_begin) * kEdgeSeedMultiplier)
        + (static_cast<unsigned long long>(level) * kEdgeRankSeedMultiplier));

    const SInt local_m = PoissonLocalCountFromScaledMass(local_mass, log_total, p, rng_, count_seed);

    if (local_m == 0) {
        return;
    }

    const SInt search_width = local_min_end - local_min_begin;

    const std::size_t expected_cache_size = static_cast<std::size_t>(
        std::min<SInt>(config_.n, std::max<SInt>(4096, std::min<SInt>(search_width, local_m * 8))));

    LogBinomCache log_binom_cache(hyperedge_size, expected_cache_size);

    GenerateLocalHyperedges(hyperedge_size, local_min_begin, local_min_end, local_m, log_binom_cache);
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateHyperedgesOfSizeExact(const SInt hyperedge_size, const double p) {
    if (hyperedge_size == 0 || hyperedge_size > config_.n || p <= 0.0) {
        return;
    }

    if (p > 1.0) {
        throw ConfigurationError("HGNP probabilities must be in [0, 1]");
    }

    const SInt local_min_begin = FindMinBoundaryByMass(config_.n, hyperedge_size, rank_, size_);

    const SInt local_min_end = FindMinBoundaryByMass(config_.n, hyperedge_size, rank_ + 1, size_);

    if (local_min_begin == local_min_end) {
        return;
    }

    const CountInt local_population = MinRangeMassExact(local_min_begin, local_min_end, config_.n, hyperedge_size);

    ValidateExactMinOwnerPartition(hyperedge_size);

    GenerateExactLocalRange(hyperedge_size, local_min_begin, local_min_end, p, 0);
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateExactLocalRange(
    const SInt hyperedge_size, const SInt local_min_begin, const SInt local_min_end, const double p, const SInt level) {
    if (local_min_begin >= local_min_end || p <= 0.0) {
        return;
    }

    const CountInt population = MinRangeMassExact(local_min_begin, local_min_end, config_.n, hyperedge_size);

    if (population == 0) {
        return;
    }

    if (population <= std::numeric_limits<SInt>::max()) {
        const SInt count_seed = sampling::Spooky::hash(
            static_cast<unsigned long long>(config_.seed)
            + (static_cast<unsigned long long>(hyperedge_size) * kCountSeedMultiplier)
            + (static_cast<unsigned long long>(rank_) * kRankSeedMultiplier)
            + (static_cast<unsigned long long>(local_min_begin) * kEdgeSeedMultiplier)
            + (static_cast<unsigned long long>(level) * kEdgeRankSeedMultiplier));

        const SInt local_population = population.convert_to<SInt>();
        const SInt local_m          = rng_.GenerateBinomial(count_seed, local_population, p);

        if (local_m == 0) {
            return;
        }

        const long double density = static_cast<long double>(local_m) / population.convert_to<long double>();

        if (density > 0.25L) {
            throw ConfigurationError("Dense exact HGNP not implemented yet");
        }

        const SInt search_width = local_min_end - local_min_begin;

        const std::size_t expected_cache_size = static_cast<std::size_t>(
            std::min<SInt>(config_.n, std::max<SInt>(4096, std::min<SInt>(search_width, local_m * 8))));

        LogBinomCache log_binom_cache(hyperedge_size, expected_cache_size);

        GenerateLocalHyperedges(hyperedge_size, local_min_begin, local_min_end, local_m, log_binom_cache);

        return;
    }

    if (local_min_end - local_min_begin <= 1) {
        throw ConfigurationError(
            "Exact HGNP min-owner atom population exceeds SInt; exact sub-minimum partitioning is required");
    }

    const SInt split = FindExactPopulationSplit(local_min_begin, local_min_end, config_.n, hyperedge_size, population);

    GenerateExactLocalRange(hyperedge_size, local_min_begin, split, p, level + 1);

    GenerateExactLocalRange(hyperedge_size, split, local_min_end, p, level + 1);
}

template <typename BigInt>
SInt HyperGNP<BigInt>::FindExactPopulationSplit(
    const SInt begin, const SInt end, const SInt n, const SInt k, const CountInt& total_population) {
    const CountInt half = total_population / 2;

    SInt left  = begin + 1;
    SInt right = end - 1;

    SInt best = begin + ((end - begin) / 2);

    while (left <= right) {
        const SInt mid = left + ((right - left) / 2);

        const CountInt left_population = MinRangeMassExact(begin, mid, n, k);

        if (left_population <= half) {
            best = mid;
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }

    best = std::clamp<SInt>(best, begin + 1, end - 1);
    return best;
}

template <typename BigInt>
void HyperGNP<BigInt>::GenerateLocalHyperedges(
    const SInt hyperedge_size, const SInt local_min_begin, const SInt local_min_end, const SInt local_m,
    LogBinomCache& log_binom_cache) {
    if (local_m == 0) {
        return;
    }
    const auto time_begin = std::chrono::steady_clock::now();

    std::unordered_set<std::vector<SInt>, VectorHash> local_seen;

    if (!config_.allow_duplicates) {
        local_seen.max_load_factor(0.5);
        local_seen.reserve(static_cast<std::size_t>(local_m));
    }

    SInt generated = 0;

    SInt edge_seed = sampling::Spooky::hash(
        static_cast<unsigned long long>(config_.seed)
        + (static_cast<unsigned long long>(rank_) * kEdgeRankSeedMultiplier)
        + (static_cast<unsigned long long>(hyperedge_size) * kEdgeSeedMultiplier));

    const CountInt max_attempts = std::max(CountInt(local_m) * 10, CountInt(1000));
    CountInt       attempts     = 0;

    while (generated < local_m) {
        const SInt s = SampleMinimumImplicit(
            local_min_begin, local_min_end, config_.n, hyperedge_size, rng_, edge_seed, log_binom_cache);

        auto pins = FloydSample(s + 1, config_.n - s - 1, hyperedge_size - 1, rng_, edge_seed);

        pins.push_back(s);
        std::sort(pins.begin(), pins.end());

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

    const SInt num_pins = local_m * hyperedge_size;

    const double us_per_pin = num_pins > 0 ? (seconds * 1e6) / static_cast<double>(num_pins) : 0.0;
    if (config_.debug) {
        std::cout << "HGNP rank=" << rank_ << " k=" << hyperedge_size << " local_m=" << local_m << " pins=" << num_pins
                  << " attempts=" << attempts << " us_per_pin=" << us_per_pin << " min_begin=" << local_min_begin
                  << " min_end=" << local_min_end << '\n';
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

template <typename BigInt>
void HyperGNP<BigInt>::FinalizeCSR(MPI_Comm /*comm*/) {}

template <typename BigInt>
void HyperGNP<BigInt>::PushHyperedge(const std::vector<SInt>& pins) {
    graph_.hyperedge_pins.insert(graph_.hyperedge_pins.end(), pins.begin(), pins.end());
    graph_.hyperedge_offsets.push_back(static_cast<SInt>(graph_.hyperedge_pins.size()));
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
template class HyperGNP<__int128>;
#pragma GCC diagnostic pop
} // namespace kagen