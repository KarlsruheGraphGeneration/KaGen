#include "kagen/generators/hyper/h_erdos/cigam.h"

#include "kagen/generators/generator.h"
#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"
#include "kagen/tools/rng_wrapper.h"

#include <set>

namespace kagen {
std::unique_ptr<Generator>
HyperCIGAMFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    if (config.n <= (1ull << 31)) {
        return std::make_unique<HyperCIGAMSmall>(config, rank, size);
    }
    return std::make_unique<HyperCIGAMBig>(config, rank, size);
}

namespace {

constexpr std::uint64_t kCIGAMPermutationDomain = 0x434947414d504552ULL; // "CIGAMPER"

std::uint64_t CIGAMPermutationSeed(const std::uint64_t global_seed) {
    return static_cast<std::uint64_t>(
        sampling::Spooky::hash(static_cast<unsigned long long>(global_seed ^ kCIGAMPermutationDomain)));
}

} // namespace

template <typename BigInt>
HyperCIGAM<BigInt>::HyperCIGAM(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size),
      vertex_permutation_(
          random_permutation::FeistelPseudoRandomPermutation::buildPermutation(
              static_cast<std::uint64_t>(config.n - 1), CIGAMPermutationSeed(static_cast<std::uint64_t>(config.seed)))),
      rng_(config) {
    if (config_.debug) {
        debug_logger_.emplace(
            config_.output_graph.filename + "_cigam_debug_pe_" + std::to_string(rank_) + ".csv", false);
    }
}

void ValidateCIGAMConfig(const PGeneratorConfig& config) {
    const SInt n = config.n;

    if (n < 2) {
        throw ConfigurationError("CIGAM requires at least two vertices");
    }

    if (!config.cigam_sizes.empty()) {
        for (const SInt k: config.cigam_sizes) {
            if (k < 2 || k > n) {
                throw ConfigurationError("Invalid CIGAM hyperedge size");
            }
        }
    } else {
        if (config.size_dist_upper_bound == 0) {
            throw ConfigurationError("CIGAM requires either -l/-u or an explicit list of hyperedge sizes");
        }

        if (config.size_dist_lower_bound < 2 || config.size_dist_lower_bound > config.size_dist_upper_bound
            || config.size_dist_upper_bound > n) {
            throw ConfigurationError("Invalid CIGAM hyperedge size range");
        }
    }

    if (config.cigam_lambda <= 0.0) {
        throw ConfigurationError("CIGAM lambda must be positive");
    }

    if (config.cigam_c.empty()) {
        throw ConfigurationError("CIGAM requires at least one density parameter");
    }
    if (config.cigam_breakpoints.size() != config.cigam_c.size()) {
        throw ConfigurationError("CIGAM requires |breakpoints| = |c| and last breakpoint 1");
    }

    if (std::abs(config.cigam_breakpoints.back() - 1.0L) > 1e-12L) {
        throw ConfigurationError("CIGAM last breakpoint must be 1");
    }
    long double prev_bp = 0.0L;
    for (const long double bp: config.cigam_breakpoints) {
        if (bp <= prev_bp || bp > 1.0L) {
            throw ConfigurationError("CIGAM breakpoints must be strictly increasing in (0, 1]");
        }
        prev_bp = bp;
    }

    long double prev_c = 1.0L;
    for (const long double c: config.cigam_c) {
        if (c <= 1.0L) {
            throw ConfigurationError("CIGAM density parameters must be > 1");
        }
        if (c <= prev_c) {
            throw ConfigurationError("CIGAM density parameters must be strictly increasing");
        }
        prev_c = c;
    }
}

template <typename BigInt>
void HyperCIGAM<BigInt>::InitQuantileGenerationState() {
    InitLayerBounds();
    InitProbabilityConstants();
    InitSizeWeights();
    if (config_.edge_budget > 0.0) {
        InitMassCache();
    }
    InitEdgeBudgetScaling();
}

// #### Entrypoints ####
template <typename BigInt>
void HyperCIGAM<BigInt>::GenerateCSR() {
    ValidateCIGAMConfig(config_);

    graph_.hyperedge_offsets.push_back(0);
    cigam_sizes_ = HyperedgeSizes();

    const auto [local_begin, local_end] = ComputeLocalVertexRange();

    SetVertexRange(local_begin, local_end);

#ifndef NDEBUG
    debug_edges_per_layer_.assign(static_cast<std::size_t>(NumLayers()), 0);
#endif

    switch (config_.cigam_mode) {
        case CIGAMMode::PAPER:
            GeneratePythonDistributedCSR();
            break;

        case CIGAMMode::EXACT:
            InitQuantileGenerationState();
            GenerateExactCSR();
            break;

        case CIGAMMode::APPROX:
            InitQuantileGenerationState();
            GenerateApproxCSR();
            break;
    }

#ifndef NDEBUG
    if (config_.cigam_mode != CIGAMMode::PAPER) {
        std::cerr << "CIGAM layer statistics:\n";

        for (SInt layer = 0; layer < NumLayers(); ++layer) {
            std::cerr << "  layer " << layer << " [" << layer_begin_[layer] << ", " << layer_end_[layer] << ") -> "
                      << debug_edges_per_layer_[static_cast<std::size_t>(layer)] << " edges\n";
        }
    }
#endif
}

template <typename BigInt>
void HyperCIGAM<BigInt>::FinalizeCSR(MPI_Comm /*comm*/) {
    if (pins_are_final_vertex_ids_) {
        return;
    }

    for (SInt& pin: graph_.hyperedge_pins) {
        if (pin < 0 || pin >= config_.n) {
            throw ConfigurationError("CIGAM generated a rank-position pin outside [0, n)");
        }

        pin = static_cast<SInt>(vertex_permutation_.f(static_cast<std::uint64_t>(pin)));
    }

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION

    if (rank_ == 0) {
        const auto& s = log_binom_stats_;

        std::cerr << "logbinom.map_calls=" << s.map_calls << '\n'
                  << "logbinom.map_hits=" << s.map_hits << '\n'
                  << "logbinom.map_misses=" << s.map_misses << '\n'
                  << "logbinom.candidate_calls=" << s.candidate_calls << '\n'
                  << "logbinom.candidate_exact_hits=" << s.candidate_exact_hits << '\n'
                  << "logbinom.candidate_recurrence_hits=" << s.candidate_recurrence_hits << '\n'
                  << "logbinom.candidate_direct_evals=" << s.candidate_direct_evals << '\n'
                  << "logbinom.distance_sum=" << s.candidate_distance_sum << '\n'
                  << "logbinom.distance_observations=" << s.candidate_distance_observations << '\n'
                  << "logbinom.max_distance=" << s.candidate_max_distance << '\n'
                  << "logbinom.forward_steps=" << s.forward_steps << '\n'
                  << "logbinom.backward_steps=" << s.backward_steps << '\n';
    }
#endif
}

// #### Setup ####
template <typename BigInt>
std::pair<SInt, SInt> HyperCIGAM<BigInt>::ComputeLocalVertexRange() {
    return BalancedVertexRange(config_.n, rank_, size_);
}

template <typename BigInt>
void HyperCIGAM<BigInt>::InitEdgeBudgetScaling() {
    if (config_.edge_budget > 0.0) {
        for (const SInt k: cigam_sizes_) {
            const long double log_Z_k = log_mass_by_size_.at(k);

            if (log_Z_k == -std::numeric_limits<long double>::infinity()) {
                throw ConfigurationError("CIGAM edge budget requested, but expected mass for size is zero");
            }

            log_edge_scaling_by_size_[k] =
                std::log(static_cast<long double>(config_.edge_budget)) + LogSizeWeight(k) - log_Z_k;
        }
    } else {
        for (const SInt k: cigam_sizes_) {
            log_edge_scaling_by_size_[k] = 0.0L;
        }
    }
}

template <typename BigInt>
void HyperCIGAM<BigInt>::InitProbabilityConstants() {
    const SInt num_layers = NumLayers();

    log_c_.resize(static_cast<std::size_t>(num_layers));

    for (SInt layer = 0; layer < num_layers; ++layer) {
        const long double c = static_cast<long double>(config_.cigam_c[layer]);

        if (!(c > 1.0L) || !std::isfinite(c)) {
            throw ConfigurationError("CIGAM density parameters must be finite and greater than one");
        }

        log_c_[static_cast<std::size_t>(layer)] = std::log(c);
    }

    const SInt local_dominant_count = config_.n > rank_ ? 1 + (config_.n - 1 - rank_) / size_ : 0;

    dominant_exponent_.resize(static_cast<std::size_t>(local_dominant_count));

    SInt local_index = 0;

    for (SInt dominant = rank_; dominant < config_.n; dominant += size_, ++local_index) {
        dominant_exponent_[static_cast<std::size_t>(local_index)] = -2.0L + RankValue(dominant);
    }
}

template <typename BigInt>
void HyperCIGAM<BigInt>::InitSizeWeights() {
    log_size_weight_.clear();

    if (config_.size_decay <= 0.0 || config_.size_decay >= 1.0) {
        const long double log_weight = -std::log(static_cast<long double>(cigam_sizes_.size()));

        for (const SInt k: cigam_sizes_) {
            log_size_weight_[k] = log_weight;
        }

        return;
    }

    const SInt        k_min = cigam_sizes_.front();
    const long double x     = static_cast<long double>(config_.size_decay);

    long double log_z = -std::numeric_limits<long double>::infinity();

    for (const SInt k: cigam_sizes_) {
        const long double log_w = std::log(x) + (static_cast<long double>(k - k_min) * std::log1pl(-x));
        log_size_weight_[k]     = log_w;
        log_z                   = LogAdd(log_z, log_w);
    }

    for (auto& [k, log_w]: log_size_weight_) {
        log_w -= log_z;
    }
}

template <typename BigInt>
void HyperCIGAM<BigInt>::InitMassCache() {
    log_mass_by_size_.clear();

    for (const SInt k: cigam_sizes_) {
        long double local_mass = 0.0L;

        LogBinomCache high_cursor(k - 1, 0);
        LogBinomCache low_cursor(k - 1, 0);

        const SInt dominant_end = config_.n - k + 1;

        for (SInt i = rank_; i < dominant_end; i += size_) {
            for (SInt layer = 0; layer < NumLayers(); ++layer) {
                const auto [j_min, j_max] = LayerEndpointRange(i, layer);

                if (j_min > j_max || j_max - i < k - 1) {
                    continue;
                }

                const long double log_mass =
                    LogBlockSize(k, i, j_min, j_max, high_cursor) + LogProbabilityForDominant(i, layer);

                local_mass += std::exp(log_mass);
            }
        }

        long double global_mass = 0.0L;

        MPI_Allreduce(&local_mass, &global_mass, 1, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        log_mass_by_size_[k] =
            global_mass > 0.0L ? std::log(global_mass) : -std::numeric_limits<long double>::infinity();
    }
}

template <typename BigInt>
long double HyperCIGAM<BigInt>::LogSizeWeight(const SInt k) const {
    return log_size_weight_.at(k);
}

template <typename BigInt>
std::vector<SInt> HyperCIGAM<BigInt>::HyperedgeSizes() const {
    if (!config_.cigam_sizes.empty()) {
        auto sizes = config_.cigam_sizes;
        std::sort(sizes.begin(), sizes.end());
        sizes.erase(std::unique(sizes.begin(), sizes.end()), sizes.end());
        return sizes;
    }

    std::vector<SInt> sizes;
    for (SInt k = config_.size_dist_lower_bound; k <= config_.size_dist_upper_bound; ++k) {
        sizes.push_back(k);
    }
    return sizes;
}
template <typename BigInt>
std::vector<double> HyperCIGAM<BigInt>::BuildDominantMassPrefix(const SInt k, const SInt layer) const {
    std::vector<double> prefix(config_.n + 1, 0.0);

    const SInt j_max = layer_end_[layer] - 1;

    /*
     * high and low form two independent monotone sequences as i grows.
     * Keep a separate recurrence cursor for each sequence.
     *
     * expected_size = 0 prevents allocating the unused hash table.
     */
    LogBinomCache high_cursor(k - 1, 0);
    LogBinomCache low_cursor(k - 1, 0);

    for (SInt i = 0; i < config_.n; ++i) {
        double mass = 0.0;

        if (config_.n - i - 1 >= k - 1) {
            const SInt j_min = std::max<SInt>(i + 1, layer_begin_[layer]);

            if (j_min <= j_max && j_max - i >= k - 1) {
                const SInt high = j_max - i;
                const SInt low  = j_min - i - 1;

                double log_block_size;

                if (k == 2) {
                    log_block_size = std::log(static_cast<double>(j_max - j_min + 1));
                } else {
                    const long double log_high = high_cursor.GetCandidate(high);

                    if (low < k - 1) {
                        log_block_size = static_cast<double>(log_high);
                    } else {
                        const long double log_low = low_cursor.GetCandidate(low);

                        log_block_size = static_cast<double>(LogDifferenceOfExponentials(log_high, log_low));
                    }
                }

                mass = std::exp(log_block_size + static_cast<double>(LogProbabilityForDominant(i, layer)));
            }
        }

        prefix[i + 1] = prefix[i] + mass;
    }

    return prefix;
}

template <typename BigInt>
SInt HyperCIGAM<BigInt>::FindDominantBoundaryInPrefix(
    const std::vector<double>& prefix, const PEID rank, const PEID size) const {
    if (rank == 0) {
        return 0;
    }

    if (rank == size) {
        return static_cast<SInt>(prefix.size() - 1);
    }

    const long double total = prefix.back();

    if (total <= 0.0L) {
        return static_cast<SInt>(prefix.size() - 1);
    }

    const long double target = total * static_cast<long double>(rank) / static_cast<long double>(size);

    auto it = std::lower_bound(prefix.begin(), prefix.end(), target);
    return static_cast<SInt>(std::distance(prefix.begin(), it));
}

// #### Exact generation ####
template <typename BigInt>
void HyperCIGAM<BigInt>::GenerateExactCSR() {
    for (const SInt k: cigam_sizes_) {
        LogBinomCache binom_cache(k - 1);
        const SInt    dominant_end = config_.n - k + 1;

        for (SInt dominant = rank_; dominant < dominant_end; dominant += size_) {
            for (SInt layer = 0; layer < NumLayers(); ++layer) {
                GenerateBoundedBlock(k, dominant, layer, binom_cache);
            }
        }

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
        AccumulateLogBinomStats(binom_cache.stats);
#endif
    }
}
template <typename BigInt>
void HyperCIGAM<BigInt>::GenerateBoundedBlock(
    const SInt k, const SInt i, const SInt layer, LogBinomCache& binom_cache) {
    const auto [j_min, j_max] = LayerEndpointRange(i, layer);

    if (j_min > j_max || j_max - i < k - 1) {
        return;
    }

    const long double log_block_size = LogBlockSize(k, i, j_min, j_max, binom_cache);

    const CountInt block_population = ExactBlockPopulation(k, i, j_min, j_max);

    if (block_population == 0) {
        return;
    }

    const long double log_p = LogProbabilityForDominant(i, layer) + log_edge_scaling_by_size_.at(k);

    const SInt count_seed = BlockCountSeed(k, i, layer);

    const SInt local_m = GenerateBinomialHybrid(block_population, log_block_size, log_p, rng_, count_seed, "CIGAM");

    if (local_m == 0) {
        return;
    }

    ValidateExactBlockDensity(local_m, log_block_size);

    SInt edge_seed = EdgeSeed(count_seed);
    mersenne_.RandomInit(edge_seed);

    auto local_seen = MakeLocalSeenSet(config_.allow_duplicates, local_m);

    SInt           generated    = 0;
    CountInt       attempts     = 0;
    const CountInt max_attempts = std::max(CountInt(local_m) * 10, CountInt(1000));

    std::vector<SInt> pins;
    pins.reserve(k);

    FloydScratchSet scratch;
    scratch.reserve(k);

    while (generated < local_m) {
        const auto event_start = std::chrono::steady_clock::now();

        SInt sampling_attempts    = 0;
        SInt duplicate_rejections = 0;

        while (true) {
            SampleBoundedHyperedge(k, i, j_min, j_max, log_block_size, binom_cache, pins, scratch);

            ++sampling_attempts;
            ++attempts;

            const bool accepted = config_.allow_duplicates || local_seen.insert(FingerprintPins(pins)).second;

            if (!accepted) {
                ++duplicate_rejections;

                if (attempts > max_attempts) {
                    throw ConfigurationError("CIGAM rejection sampling exceeded attempt limit; use --fast");
                }

                continue;
            }

            PushCheckedHyperedge(pins, layer);
            ++generated;

            if (debug_logger_) {
                const auto duration_ns =
                    std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - event_start)
                        .count();

                debug_logger_->LogHyperedge({
                    .hyperedge_id         = next_debug_hyperedge_id_++,
                    .mode                 = static_cast<SInt>(config_.cigam_mode),
                    .hyperedge_size       = k,
                    .dominant_vertex      = i,
                    .layer                = layer,
                    .endpoint_vertex      = pins.back(),
                    .block_j_min          = j_min,
                    .block_j_max          = j_max,
                    .block_log_population = log_block_size,
                    .log_probability      = LogProbabilityForDominant(i, layer) + log_edge_scaling_by_size_.at(k),
                    .sampling_attempts    = sampling_attempts,
                    .duplicate_rejections = duplicate_rejections,
                    .duration_ns          = duration_ns,
                });
            }

            break;
        }
    }
}

template <typename BigInt>
void HyperCIGAM<BigInt>::ValidateExactBlockDensity(const SInt local_m, const long double log_block_size) {
    if (!config_.allow_duplicates) {
        const long double density = static_cast<long double>(local_m) / expl(log_block_size);

        if (density > 0.25L) {
            throw ConfigurationError("Dense exact CIGAM block not implemented yet; use --fast");
        }
    }
}

// #### Approx generation ####
template <typename BigInt>
void HyperCIGAM<BigInt>::GenerateApproxCSR() {
    for (const SInt k: cigam_sizes_) {
        LogBinomCache cache(k - 1, std::min<SInt>(config_.n, SInt{1} << 20));

        const SInt dominant_end = config_.n - k + 1;

        std::vector<BlockInfo> blocks;
        blocks.reserve(NumLayers());

        for (SInt dominant = rank_; dominant < dominant_end; dominant += size_) {
            blocks.clear();

            for (SInt layer = 0; layer < NumLayers(); ++layer) {
                const auto [j_min, j_max] = LayerEndpointRange(dominant, layer);

                if (j_min > j_max || j_max - dominant < k - 1) {
                    continue;
                }

                blocks.push_back({layer, j_min, j_max, LogBlockSize(k, dominant, j_min, j_max, cache)});
            }

            for (const auto& block: blocks) {
                GenerateApproxBlock(k, dominant, block.layer, block.j_min, block.j_max, block.log_block_size, cache);
            }
        }
    }
}

template <typename BigInt>
void HyperCIGAM<BigInt>::GenerateApproxBlock(
    const SInt k, const SInt dominant, const SInt layer, SInt j_min, SInt j_max, SInt log_block_size,
    LogBinomCache& cache) {
    const long double log_expected =
        log_block_size + LogProbabilityForDominant(dominant, layer) + log_edge_scaling_by_size_.at(k);

    const long double smallest = std::log(static_cast<long double>(std::numeric_limits<double>::denorm_min()));

    if (log_expected <= smallest) {
        return;
    }

    if (log_expected > std::log(static_cast<long double>(std::numeric_limits<SInt>::max()))) {
        throw ConfigurationError("CIGAM approximate block expected count exceeds SInt");
    }

    const SInt count_seed = BlockCountSeed(k, dominant, layer);

    const SInt local_m = rng_.GeneratePoisson(count_seed, static_cast<double>(std::exp(log_expected)));

    if (local_m == 0) {
        return;
    }

    mersenne_.RandomInit(EdgeSeed(count_seed));

    std::vector<SInt> pins;
    pins.reserve(k);

    FloydScratchSet scratch;
    scratch.reserve(k);

    const long double log_probability = LogProbabilityForDominant(dominant, layer) + log_edge_scaling_by_size_.at(k);

    for (SInt e = 0; e < local_m; ++e) {
        const auto event_start = std::chrono::steady_clock::now();

        pins.clear();
        pins.push_back(dominant);

        const SInt endpoint = SampleEndpoint(k, dominant, j_min, j_max, log_block_size, cache);

        if (k > 2) {
            FloydSampleInto(dominant + 1, endpoint - dominant - 1, k - 2, mersenne_, pins, scratch, 1);
        }

        pins.push_back(endpoint);

        PushCheckedHyperedge(pins, layer);

        if (debug_logger_) {
            const auto duration_ns =
                std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - event_start)
                    .count();

            debug_logger_->LogHyperedge({
                .hyperedge_id         = next_debug_hyperedge_id_++,
                .mode                 = static_cast<SInt>(config_.cigam_mode),
                .hyperedge_size       = k,
                .dominant_vertex      = dominant,
                .layer                = layer,
                .endpoint_vertex      = endpoint,
                .block_j_min          = j_min,
                .block_j_max          = j_max,
                .block_log_population = static_cast<long double>(log_block_size),
                .log_probability      = log_probability,
                .sampling_attempts    = 1,
                .duplicate_rejections = 0,
                .duration_ns          = duration_ns,
            });
        }
    }
}

// Recursively partition the dominant-vertex interval until each subproblem
// has either sufficiently small expected mass or a sufficiently small
// vertex range. Leaf ranges sample their local Poisson edge count
// directly and then draw dominant vertices according to their relative
// contribution within the leaf

// #### Math/ Layer helpers ####
template <typename BigInt>
long double HyperCIGAM<BigInt>::RankQuantile(const SInt position) const {
    if (position < 0 || position >= config_.n) {
        throw ConfigurationError("CIGAM rank position is outside [0, n)");
    }

    const long double quantile =
        1.0L - ((static_cast<long double>(position) + 0.5L) / static_cast<long double>(config_.n));

    /*
     * The midpoint expression is mathematically inside (0, 1).
     * Clamp only to protect against floating-point rounding for
     * extremely large n.
     */
    return std::clamp(quantile, std::nextafter(0.0L, 1.0L), std::nextafter(1.0L, 0.0L));
}

template <typename BigInt>
long double HyperCIGAM<BigInt>::RankValue(const SInt position) const {
    return InverseTruncatedExpCDF(RankQuantile(position));
}
template <typename BigInt>
long double HyperCIGAM<BigInt>::InverseTruncatedExpCDF(const long double u) const {
    const long double lambda = static_cast<long double>(config_.cigam_lambda);

    if (!std::isfinite(lambda) || lambda <= 0.0L) {
        throw ConfigurationError("CIGAM lambda must be finite and positive");
    }

    if (!std::isfinite(u) || u < 0.0L || u > 1.0L) {
        throw ConfigurationError("CIGAM quantile must be finite and lie in [0, 1]");
    }

    const long double normalization = -std::expm1l(-lambda);

    const long double rank = -std::log1pl(-u * normalization) / lambda;

    return std::clamp(rank, 0.0L, 1.0L);
}

template <typename BigInt>
void HyperCIGAM<BigInt>::InitLayerBounds() {
    const SInt n          = config_.n;
    const SInt num_layers = NumLayers();

    layer_begin_.assign(static_cast<std::size_t>(num_layers), 0);

    layer_end_.assign(static_cast<std::size_t>(num_layers), 0);

    SInt previous_end = 0;

    for (SInt layer = 0; layer < num_layers; ++layer) {
        const SInt end = layer == num_layers - 1
                             ? n
                             : FirstPositionAfterBreakpoint(static_cast<long double>(config_.cigam_breakpoints[layer]));

        if (end < previous_end || end > n) {
            throw ConfigurationError("CIGAM produced inconsistent layer boundaries");
        }

        layer_begin_[static_cast<std::size_t>(layer)] = previous_end;

        layer_end_[static_cast<std::size_t>(layer)] = end;

        previous_end = end;
    }

    if (previous_end != n) {
        throw ConfigurationError("CIGAM layers do not cover all rank positions");
    }
}

template <typename BigInt>
SInt HyperCIGAM<BigInt>::FirstPositionAfterBreakpoint(const long double breakpoint) const {
    const SInt n = config_.n;

    if (breakpoint < 0.0L || breakpoint > 1.0L || !std::isfinite(breakpoint)) {
        throw ConfigurationError("CIGAM breakpoint must be finite and lie in [0, 1]");
    }

    if (breakpoint >= 1.0L) {
        return n;
    }

    const long double lambda = static_cast<long double>(config_.cigam_lambda);

    const long double target_rank = 1.0L - breakpoint;

    const long double cdf = (-std::expm1l(-lambda * target_rank)) / (-std::expm1l(-lambda));

    /*
     * Initial inversion of:
     *
     *   u_i = 1 - (i + 1/2) / n < cdf.
     */
    const long double estimate = (static_cast<long double>(n) * (1.0L - cdf)) - 0.5L;

    SInt candidate;

    if (estimate < 0.0L) {
        candidate = 0;
    } else if (estimate >= static_cast<long double>(n)) {
        candidate = n;
    } else {
        candidate = static_cast<SInt>(std::floor(estimate)) + 1;
        candidate = std::clamp<SInt>(candidate, 0, n);
    }

    auto is_after = [&](const SInt position) {
        return position < n && (1.0L - RankValue(position)) > breakpoint;
    };

    /*
     * Correct any rounding error near the analytic boundary.
     * In normal operation this performs zero or one iteration.
     */
    while (candidate > 0 && is_after(candidate - 1)) {
        --candidate;
    }

    while (candidate < n && !is_after(candidate)) {
        ++candidate;
    }

    return candidate;
}

template <typename BigInt>
std::pair<SInt, SInt> HyperCIGAM<BigInt>::LayerEndpointRange(const SInt i, const SInt layer) const {
    const SInt j_min = std::max<SInt>(i + 1, layer_begin_[layer]);
    const SInt j_end = layer_end_[layer];

    if (j_min >= j_end) {
        return {1, 0};
    }

    return {j_min, j_end - 1};
}

template <typename BigInt>
long double HyperCIGAM<BigInt>::LogProbabilityForDominant(const SInt i, const SInt layer) const {
    if (i < 0 || i >= config_.n) {
        throw ConfigurationError("CIGAM dominant position is outside [0, n)");
    }

    if (layer < 0 || layer >= NumLayers()) {
        throw ConfigurationError("CIGAM layer is outside its valid range");
    }

#ifndef NDEBUG
    if (i % size_ != rank_) {
        throw ConfigurationError("CIGAM accessed a dominant not owned by this PE");
    }
#endif

    const SInt local_index = (i - rank_) / size_;

    const long double log_probability =
        dominant_exponent_[static_cast<std::size_t>(local_index)] * log_c_[static_cast<std::size_t>(layer)];

    if (!std::isfinite(log_probability) || log_probability >= 0.0L) {
        throw ConfigurationError("CIGAM generated an invalid hyperedge probability");
    }

    return log_probability;
}

template <typename BigInt>
long double HyperCIGAM<BigInt>::LogBlockSize(
    const SInt k, const SInt i, const SInt j_min, const SInt j_max, LogBinomCache& binom_cache) const {
    if (k == 2) {
        return std::log(static_cast<long double>(j_max - j_min + 1));
    }

    const SInt high = j_max - i;
    const SInt low  = j_min - i - 1;

    const long double log_high = binom_cache.Get(high, k - 1);
    const long double log_low =
        low >= k - 1 ? binom_cache.Get(low, k - 1) : -std::numeric_limits<long double>::infinity();

    return LogDifferenceOfExponentials(log_high, log_low);
}

// #### Endpoint/hyperedge sampling ####
template <typename BigInt>
void HyperCIGAM<BigInt>::SampleBoundedHyperedge(
    const SInt k, const SInt dominant, const SInt j_min, const SInt j_max, const long double log_block_size,
    LogBinomCache& binom_cache, std::vector<SInt>& pins, FloydScratchSet& scratch) {
    pins.clear();
    pins.push_back(dominant);

    const SInt endpoint = SampleEndpoint(k, dominant, j_min, j_max, log_block_size, binom_cache);

    if (k > 2) {
        FloydSampleInto(dominant + 1, endpoint - dominant - 1, k - 2, mersenne_, pins, scratch, 1);
    }

    pins.push_back(endpoint);
}

template <typename BigInt>
SInt HyperCIGAM<BigInt>::SampleEndpoint(
    const SInt k, const SInt i, const SInt j_min, const SInt j_max, const long double log_block_size,
    LogBinomCache& binom_cache) {
    assert(k >= 2);
    assert(i < j_min);
    assert(j_min <= j_max);

    const long double u = std::min<long double>(
        static_cast<long double>(rng_.GenerateCanonicalDoubleStream()), std::nextafter(1.0L, 0.0L));

    if (k == 2) {
        return j_min + static_cast<SInt>(u * static_cast<long double>(j_max - j_min + 1));
    }

    if (u <= 0.0L) {
        return j_min;
    }

    const SInt q = k - 1;

    const long double log_target = std::log(u) + log_block_size;

    const SInt low = j_min - i - 1;

    const long double log_low = low >= q ? binom_cache.Get(low, q) : -std::numeric_limits<long double>::infinity();

    const bool use_linear_prefix = q <= 8 && (j_max - i) < (SInt{1} << 30);

    const long double low_mass = use_linear_prefix && low >= q ? BinomSmallExactLD(low, q) : 0.0L;

    const auto log_prefix = [&](const SInt j) {
        if (use_linear_prefix) {
            const long double prefix = BinomSmallExactLD(j - i, q) - low_mass;

            return prefix > 0.0L ? std::log(prefix) : -std::numeric_limits<long double>::infinity();
        }

        return LogDifferenceOfExponentials(binom_cache.Get(j - i, q), log_low);
    };

    SInt j = EstimateEndpointInitialGuess(k, i, j_min, j_max, log_target, binom_cache);

    constexpr SInt kMaxCorrectionSteps = 64;
    SInt           steps               = 0;

    while (j > j_min && log_prefix(j - 1) >= log_target) {
        --j;

        if (++steps > kMaxCorrectionSteps) {
            return SampleEndpointBinarySearch(k, i, j_min, j_max, log_target, binom_cache);
        }
    }

    while (j < j_max && log_prefix(j) < log_target) {
        ++j;

        if (++steps > kMaxCorrectionSteps) {
            return SampleEndpointBinarySearch(k, i, j_min, j_max, log_target, binom_cache);
        }
    }

    return j;
}

template <typename BigInt>
SInt HyperCIGAM<BigInt>::EstimateEndpointInitialGuess(
    const SInt k, const SInt i, const SInt j_min, const SInt j_max, const long double log_target_prefix,
    LogBinomCache& binom_cache) {
    assert(k > 2);
    assert(i < j_min);
    assert(j_min <= j_max);

    const SInt q   = k - 1;
    const SInt low = j_min - i - 1;

    const long double log_low_mass = low >= q ? binom_cache.Get(low, q) : -std::numeric_limits<long double>::infinity();

    const long double log_absolute_target = log_low_mass == -std::numeric_limits<long double>::infinity()
                                                ? log_target_prefix
                                                : LogAdd(log_low_mass, log_target_prefix);

    if (!std::isfinite(log_absolute_target)) {
        return j_min;
    }

    const long double log_q_factorial = std::lgammal(static_cast<long double>(q) + 1.0L);

    const long double root = std::exp((log_absolute_target + log_q_factorial) / static_cast<long double>(q));

    const long double x_estimate = root + (0.5L * static_cast<long double>(q - 1));

    const long double j_estimate = static_cast<long double>(i) + x_estimate;

    if (!std::isfinite(j_estimate) || j_estimate <= static_cast<long double>(j_min)) {
        return j_min;
    }

    if (j_estimate >= static_cast<long double>(j_max)) {
        return j_max;
    }

    return std::clamp<SInt>(static_cast<SInt>(std::llround(j_estimate)), j_min, j_max);
}

template <typename BigInt>
SInt HyperCIGAM<BigInt>::SampleEndpointBinarySearch(
    const SInt k, const SInt i, const SInt j_min, const SInt j_max, const long double log_target,
    LogBinomCache& binom_cache) {
    const SInt q = k - 1;

    const SInt low = j_min - i - 1;

    const long double log_low = low >= q ? binom_cache.Get(low, q) : -std::numeric_limits<long double>::infinity();

    SInt left  = j_min;
    SInt right = j_max;

    while (left < right) {
        const SInt mid = left + ((right - left) / 2);

        const long double log_prefix = LogDifferenceOfExponentials(binom_cache.Get(mid - i, q), log_low);

        if (log_prefix >= log_target) {
            right = mid;
        } else {
            left = mid + 1;
        }
    }

    return left;
}

// #### Shared sampling / graph helpers
template <typename BigInt>
SInt HyperCIGAM<BigInt>::BlockCountSeed(const SInt k, const SInt i, const SInt layer) {
    return sampling::Spooky::hash(
        static_cast<unsigned long long>(config_.seed) + (static_cast<unsigned long long>(k) * kCountSeedMultiplier)
        + (static_cast<unsigned long long>(i) * kEdgeSeedMultiplier)
        + (static_cast<unsigned long long>(layer) * kRankSeedMultiplier));
}

template <typename BigInt>
SInt HyperCIGAM<BigInt>::EdgeSeed(const SInt count_seed) {
    return sampling::Spooky::hash(count_seed + kEdgeRankSeedMultiplier);
}

template <typename BigInt>
void HyperCIGAM<BigInt>::PushCheckedHyperedge(const std::vector<SInt>& pins, const SInt layer) {
    AssertHyperedgeInvariants(pins, layer);

#ifndef NDEBUG
    ++debug_edges_per_layer_[layer];
#endif

    PushUncompressedHyperedge(graph_, memory_stats_, pins);
}

template <typename BigInt>
void HyperCIGAM<BigInt>::AssertHyperedgeInvariants(const std::vector<SInt>& pins, const SInt layer) const {
#ifndef NDEBUG
    const SInt j = pins.back();

    if (j < layer_begin_[layer] || j >= layer_end_[layer]) {
        throw ConfigurationError("CIGAM layer invariant violated: endpoint outside layer");
    }

    for (SInt pos = 1; pos < static_cast<SInt>(pins.size()); ++pos) {
        if (pins[pos - 1] >= pins[pos]) {
            throw ConfigurationError("CIGAM invariant violated: pins not strictly increasing before permutation");
        }
    }
#endif
}

template <typename BigInt>
long double HyperCIGAM<BigInt>::CounterUniform01(const std::uint64_t stream_tag, const std::uint64_t index) const {
    const std::uint64_t key = static_cast<std::uint64_t>(config_.seed) ^ (stream_tag * 0x9e3779b97f4a7c15ULL)
                              ^ (index * 0xd1b54a32d192ed03ULL);

    const std::uint64_t random_word =
        static_cast<std::uint64_t>(sampling::Spooky::hash(static_cast<unsigned long long>(key)));

    // Construct a binary64-compatible uniform value in [0, 1).
    constexpr long double kInverseTwoTo53 = 1.0L / 9007199254740992.0L;

    return static_cast<long double>(random_word >> 11) * kInverseTwoTo53;
}

template <typename BigInt>
typename HyperCIGAM<BigInt>::PythonReferenceRanks HyperCIGAM<BigInt>::GeneratePythonReferenceRanks() const {
    constexpr std::uint64_t kRankStreamTag = 0x434947414d52414eULL; // "CIGAMRAN"

    const SInt n = config_.n;
    const SInt L = NumLayers();

    struct RankedVertex {
        long double rank;
        SInt        vertex;
    };

    std::vector<RankedVertex> sampled;
    sampled.reserve(static_cast<std::size_t>(n));

    const long double lambda = static_cast<long double>(config_.cigam_lambda);

    const long double truncated_mass = 1.0L - std::exp(-lambda);

    for (SInt vertex = 0; vertex < n; ++vertex) {
        const long double u = CounterUniform01(kRankStreamTag, static_cast<std::uint64_t>(vertex));

        // Inverse CDF of the truncated exponential distribution:
        //
        // F(r) = (1 - exp(-lambda r)) / (1 - exp(-lambda)).
        const long double rank = -std::log1pl(-u * truncated_mass) / lambda;

        sampled.push_back({
            .rank   = rank,
            .vertex = vertex,
        });
    }

    // Dominant position 0 has the greatest prestige.
    std::sort(sampled.begin(), sampled.end(), [](const RankedVertex& lhs, const RankedVertex& rhs) {
        if (lhs.rank != rhs.rank) {
            return lhs.rank > rhs.rank;
        }

        // Deterministic tie-breaking, although exact ties should be rare.
        return lhs.vertex < rhs.vertex;
    });

    PythonReferenceRanks result;

    result.original_vertex.resize(static_cast<std::size_t>(n));
    result.rank.resize(static_cast<std::size_t>(n));

    result.layer_begin.assign(static_cast<std::size_t>(L), n);
    result.layer_end.assign(static_cast<std::size_t>(L), n);

    for (SInt position = 0; position < n; ++position) {
        result.original_vertex[position] = sampled[position].vertex;

        result.rank[position] = sampled[position].rank;

        // The paper defines the layer from 1 - minimum rank.
        // Since ranks are decreasing by position, this value increases.
        const long double layer_coordinate = 1.0L - sampled[position].rank;

        const auto it =
            std::lower_bound(config_.cigam_breakpoints.begin(), config_.cigam_breakpoints.end(), layer_coordinate);

        if (it == config_.cigam_breakpoints.end()) {
            throw ConfigurationError("CIGAM sampled rank is not covered by a breakpoint");
        }

        const SInt layer = static_cast<SInt>(std::distance(config_.cigam_breakpoints.begin(), it));

        if (result.layer_begin[layer] == n) {
            result.layer_begin[layer] = position;
        }

        result.layer_end[layer] = position + 1;
    }

    // Represent empty layers as empty intervals.
    for (SInt layer = 0; layer < L; ++layer) {
        if (result.layer_begin[layer] == n) {
            result.layer_begin[layer] = 0;
            result.layer_end[layer]   = 0;
        }
    }

    return result;
}

template <typename BigInt>
PEID HyperCIGAM<BigInt>::BlockOwner(const SInt k, const SInt dominant, const SInt layer) const {
    const std::uint64_t key =
        static_cast<std::uint64_t>(config_.seed) ^ (static_cast<std::uint64_t>(k) * 0x9e3779b97f4a7c15ULL)
        ^ (static_cast<std::uint64_t>(dominant) * 0xd1b54a32d192ed03ULL)
        ^ (static_cast<std::uint64_t>(layer) * 0x94d049bb133111ebULL) ^ 0x424c4f434b4f574eULL; // "BLOCKOWN"

    const std::uint64_t hash = static_cast<std::uint64_t>(sampling::Spooky::hash(static_cast<unsigned long long>(key)));

    return static_cast<PEID>(hash % static_cast<std::uint64_t>(size_));
}

template <typename BigInt>
SInt HyperCIGAM<BigInt>::PythonBlockCountSeed(const SInt k, const SInt dominant, const SInt layer) const {
    const std::uint64_t key =
        static_cast<std::uint64_t>(config_.seed) + (static_cast<std::uint64_t>(k) * kCountSeedMultiplier)
        + (static_cast<std::uint64_t>(dominant) * kEdgeSeedMultiplier)
        + (static_cast<std::uint64_t>(layer) * kRankSeedMultiplier) + 0x434f554e54535452ULL; // "COUNTSTR"

    return static_cast<SInt>(sampling::Spooky::hash(static_cast<unsigned long long>(key)));
}

template <typename BigInt>
SInt HyperCIGAM<BigInt>::PythonBlockEdgeSeed(const SInt k, const SInt dominant, const SInt layer) const {
    const std::uint64_t key =
        static_cast<std::uint64_t>(PythonBlockCountSeed(k, dominant, layer)) ^ 0x454447455354524dULL; // "EDGESTRM"

    return static_cast<SInt>(sampling::Spooky::hash(static_cast<unsigned long long>(key)));
}

template <typename BigInt>
CountInt
HyperCIGAM<BigInt>::PythonBlockPopulation(const SInt k, const SInt dominant, const SInt j_min, const SInt j_max) const {
    if (j_min > j_max) {
        return 0;
    }

    if (k < 2) {
        throw ConfigurationError("CIGAM requires hyperedge size at least two");
    }

    const SInt q = k - 1;

    const SInt high = j_max - dominant;
    const SInt low  = j_min - dominant - 1;

    if (high < q) {
        return 0;
    }

    const CountInt upper = BinomialExact(high, q);

    const CountInt lower = low >= q ? BinomialExact(low, q) : CountInt{0};

    const CountInt population = upper - lower;

    if (population < 0) {
        throw ConfigurationError("Negative CIGAM block population");
    }

    return population;
}

template <typename BigInt>
std::vector<SInt> HyperCIGAM<BigInt>::SamplePythonBlockHyperedge(
    const SInt k, const SInt dominant, const SInt j_min, const SInt j_max, SInt& sampling_attempts,
    SInt& block_rejections) {
    if (j_min > j_max) {
        throw ConfigurationError("Cannot sample from an empty CIGAM block");
    }

    const SInt suffix_universe_size = j_max - dominant;
    const SInt suffix_sample_size   = k - 1;

    if (suffix_sample_size > suffix_universe_size) {
        throw ConfigurationError("CIGAM block does not contain enough vertices");
    }

    while (true) {
        ++sampling_attempts;

        auto suffix = FloydSample(dominant + 1, suffix_universe_size, suffix_sample_size, mersenne_);

        if (suffix.empty() || suffix.back() < j_min) {
            ++block_rejections;
            continue;
        }

        std::vector<SInt> pins;
        pins.reserve(static_cast<std::size_t>(k));

        pins.push_back(dominant);
        pins.insert(pins.end(), suffix.begin(), suffix.end());

        return pins;
    }
}

template <typename BigInt>
void HyperCIGAM<BigInt>::GeneratePythonBlock(
    const SInt k, const SInt dominant, const SInt layer, const SInt j_min, const SInt j_max,
    const PythonReferenceRanks& ranks) {
    const CountInt population = PythonBlockPopulation(k, dominant, j_min, j_max);

    if (population == 0) {
        return;
    }

    const long double log_probability =
        (-2.0L + ranks.rank[dominant]) * std::log(static_cast<long double>(config_.cigam_c[layer]));

    const SInt count_seed = PythonBlockCountSeed(k, dominant, layer);

    const SInt block_count =
        GenerateBinomialHybrid(population, log_probability, rng_, count_seed, "Python-reference CIGAM");

    if (block_count == 0) {
        return;
    }

    if (CountInt(block_count) > population) {
        throw ConfigurationError("CIGAM block count exceeds block population");
    }

    mersenne_.RandomInit(PythonBlockEdgeSeed(k, dominant, layer));

    struct PendingPaperDebugData {
        SInt sampling_attempts    = 0;
        SInt duplicate_rejections = 0;

        std::chrono::steady_clock::duration duration{};
    };

    std::map<std::vector<SInt>, PendingPaperDebugData> selected_rank_positions;

    SInt pending_sampling_attempts    = 0;
    SInt pending_duplicate_rejections = 0;
    auto pending_event_start          = std::chrono::steady_clock::now();

    while (static_cast<SInt>(selected_rank_positions.size()) < block_count) {
        SInt candidate_sampling_attempts = 0;
        SInt candidate_block_rejections  = 0;

        auto pins = SamplePythonBlockHyperedge(
            k, dominant, j_min, j_max, candidate_sampling_attempts, candidate_block_rejections);

        pending_sampling_attempts += candidate_sampling_attempts;

        const auto [it, inserted] = selected_rank_positions.try_emplace(std::move(pins));

        if (!inserted) {
            ++pending_duplicate_rejections;
            continue;
        }

        it->second = {
            .sampling_attempts    = pending_sampling_attempts,
            .duplicate_rejections = pending_duplicate_rejections,
            .duration             = std::chrono::steady_clock::now() - pending_event_start,
        };

        pending_sampling_attempts    = 0;
        pending_duplicate_rejections = 0;
        pending_event_start          = std::chrono::steady_clock::now();
    }

    const long double block_log_population = std::log(population.convert_to<long double>());

    for (const auto& [rank_positions, debug_data]: selected_rank_positions) {
        std::vector<SInt> original_pins;
        original_pins.reserve(rank_positions.size());

        for (const SInt position: rank_positions) {
            original_pins.push_back(ranks.original_vertex[position]);
        }

        std::sort(original_pins.begin(), original_pins.end());

        PushUncompressedHyperedge(graph_, memory_stats_, original_pins);

        if (debug_logger_) {
            const auto duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(debug_data.duration).count();

            debug_logger_->LogHyperedge({
                .hyperedge_id         = next_debug_hyperedge_id_++,
                .mode                 = static_cast<SInt>(config_.cigam_mode),
                .hyperedge_size       = k,
                .dominant_vertex      = dominant,
                .layer                = layer,
                .endpoint_vertex      = rank_positions.back(),
                .block_j_min          = j_min,
                .block_j_max          = j_max,
                .block_log_population = block_log_population,
                .log_probability      = log_probability,
                .sampling_attempts    = debug_data.sampling_attempts,
                .duplicate_rejections = debug_data.duplicate_rejections,
                .duration_ns          = duration_ns,
            });
        }
    }
}

template <typename BigInt>
void HyperCIGAM<BigInt>::GeneratePythonDistributedCSR() {
    ValidateCIGAMConfig(config_);

    if (config_.edge_budget > 0.0) {
        throw ConfigurationError(
            "The Python-reference CIGAM implementation does not "
            "support edge-budget scaling");
    }

    if (config_.size_decay > 0.0 && config_.size_decay < 1.0) {
        throw ConfigurationError(
            "The Python-reference CIGAM implementation does not "
            "support KaGen size-decay weighting");
    }

    if (config_.allow_duplicates) {
        throw ConfigurationError(
            "The CIGAM Bernoulli model generates a simple hypergraph "
            "and does not allow duplicate hyperedges");
    }

    pins_are_final_vertex_ids_ = true;

    const PythonReferenceRanks ranks = GeneratePythonReferenceRanks();

    const std::vector<SInt> sizes = HyperedgeSizes();

    const SInt n = config_.n;
    const SInt L = NumLayers();

#ifndef NDEBUG
    debug_edges_per_layer_.assign(static_cast<std::size_t>(L), 0);
#endif

    for (const SInt k: sizes) {
        if (k < 2 || k > n) {
            throw ConfigurationError("Invalid CIGAM hyperedge size");
        }

        // The dominant vertex must leave room for k - 1 later positions.
        const SInt dominant_end = n - k + 1;

        for (SInt dominant = 0; dominant < dominant_end; ++dominant) {
            for (SInt layer = 0; layer < L; ++layer) {
                if (BlockOwner(k, dominant, layer) != rank_) {
                    continue;
                }

                const SInt layer_begin = ranks.layer_begin[layer];

                const SInt layer_end = ranks.layer_end[layer];

                if (layer_begin >= layer_end) {
                    continue;
                }

                const SInt j_min = std::max<SInt>(dominant + 1, layer_begin);

                const SInt j_max = layer_end - 1;

                if (j_min > j_max) {
                    continue;
                }

                if (j_max - dominant < k - 1) {
                    continue;
                }

                const std::size_t offsets_before = graph_.hyperedge_offsets.size();

                GeneratePythonBlock(k, dominant, layer, j_min, j_max, ranks);

#ifndef NDEBUG
                const std::size_t offsets_after = graph_.hyperedge_offsets.size();

                debug_edges_per_layer_[layer] += static_cast<SInt>(offsets_after - offsets_before);
#endif
            }
        }
    }

    // This describes vertex-data ownership, not block ownership.
    const auto [local_begin, local_end] = BalancedVertexRange(config_.n, rank_, size_);

    SetVertexRange(local_begin, local_end);

#ifndef NDEBUG
    std::cerr << "PE " << rank_ << " generated " << graph_.hyperedge_offsets.size() - 1
              << " Python-reference CIGAM hyperedges\n";
#endif
}

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
template <typename BigInt>
void HyperCIGAM<BigInt>::AccumulateLogBinomStats(const LogBinomCacheStats& stats) {
    log_binom_stats_.map_calls += stats.map_calls;
    log_binom_stats_.map_hits += stats.map_hits;
    log_binom_stats_.map_misses += stats.map_misses;
    log_binom_stats_.map_inserts += stats.map_inserts;

    log_binom_stats_.candidate_calls += stats.candidate_calls;
    log_binom_stats_.candidate_exact_hits += stats.candidate_exact_hits;
    log_binom_stats_.candidate_recurrence_hits += stats.candidate_recurrence_hits;
    log_binom_stats_.candidate_direct_evals += stats.candidate_direct_evals;
    log_binom_stats_.candidate_below_range += stats.candidate_below_range;

    log_binom_stats_.candidate_distance_sum += stats.candidate_distance_sum;
    log_binom_stats_.candidate_distance_observations += stats.candidate_distance_observations;
    log_binom_stats_.candidate_max_distance =
        std::max(log_binom_stats_.candidate_max_distance, stats.candidate_max_distance);

    log_binom_stats_.forward_steps += stats.forward_steps;
    log_binom_stats_.backward_steps += stats.backward_steps;

    log_binom_stats_.min_key  = std::min(log_binom_stats_.min_key, stats.min_key);
    log_binom_stats_.max_key  = std::max(log_binom_stats_.max_key, stats.max_key);
    log_binom_stats_.max_size = std::max(log_binom_stats_.max_size, stats.max_size);
}
#endif

template <typename BigInt>
CountInt
HyperCIGAM<BigInt>::ExactBlockPopulation(const SInt k, const SInt i, const SInt j_min, const SInt j_max) const {
    if (j_min > j_max || j_max - i < k - 1) {
        return 0;
    }

    const SInt q = k - 1;

    const SInt high = j_max - i;

    const SInt low = j_min - i - 1;

    const CountInt upper = BinomialExact(high, q);

    const CountInt lower = low >= q ? BinomialExact(low, q) : CountInt{0};

    const CountInt population = upper - lower;

    if (population < 0) {
        throw ConfigurationError("CIGAM block population is negative");
    }

    return population;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
template class HyperCIGAM<__int128>;
template class HyperCIGAM<SInt>;
#pragma GCC diagnostic pop

} // namespace kagen
