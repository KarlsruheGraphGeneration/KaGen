#include "kagen/generators/hyper/h_erdos/cigam.h"

#include "kagen/generators/generator.h"
#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"
#include "kagen/tools/rng_wrapper.h"

#include <memory>
#include <tuple>
namespace kagen {
std::unique_ptr<Generator>
HyperCIGAMFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    if (config.n <= (1ull << 31)) {
        return std::make_unique<HyperCIGAMSmall>(config, rank, size);
    }
    return std::make_unique<HyperCIGAMBig>(config, rank, size);
}

void NormalizeParams(const PGeneratorConfig& config) {
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
void HyperCIGAM<BigInt>::GenerateExactCSR(SInt local_begin, SInt local_end) {
    SetVertexRange(local_begin, local_end);

    for (const SInt k: HyperedgeSizes()) {
        LogBinomCache binom_cache(k - 1);

        for (SInt i = local_begin; i < local_end; ++i) {
            if (config_.n - i - 1 < k - 1) {
                continue;
            }

            for (SInt layer = 0; layer < NumLayers(); ++layer) {
                GenerateBoundedBlock(k, i, layer, binom_cache);
            }
        }
    }
}

template <typename BigInt>
void HyperCIGAM<BigInt>::GenerateApproxCSR(SInt local_begin, SInt local_end) {
    for (const SInt k: HyperedgeSizes()) {
        LogBinomCache cache(k - 1, std::min<SInt>(config_.n, SInt{1} << 20));

        for (SInt layer = 0; layer < NumLayers(); ++layer) {
            GenerateApproxRange(k, local_begin, local_end, layer, 0, cache);
        }
    }
}

template <typename BigInt>
SInt HyperCIGAM<BigInt>::FindApproxMassSplit(
    const SInt k, const SInt i_begin, const SInt i_end, const SInt layer, const long double log_total_mass,
    LogBinomCache& cache) {
    if (i_end - i_begin <= 2) {
        return i_begin + 1;
    }
    SInt left  = i_begin + 1;
    SInt right = i_end - 1;
    SInt best  = i_begin + ((i_end - i_begin) / 2);

    const long double log_half = log_total_mass + std::log(0.5L);

    while (left <= right) {
        const SInt mid = left + ((right - left) / 2);

        const long double log_left_mass = LogExpectedRangeMass(k, i_begin, mid, layer, cache);

        if (log_left_mass <= log_half) {
            best = mid;
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }

    return std::clamp<SInt>(best, i_begin + 1, i_end - 1);
}

template <typename BigInt>
void HyperCIGAM<BigInt>::InitProbabilityConstants() {
    const SInt L = NumLayers();

    log_c_.resize(L);
    log_c_over_lambda_.resize(L);

    const long double lambda = static_cast<long double>(config_.cigam_lambda);
    lambda_exp_term_         = 1.0L - expl(-lambda);

    for (SInt layer = 0; layer < L; ++layer) {
        log_c_[layer]             = std::log(static_cast<long double>(config_.cigam_c[layer]));
        log_c_over_lambda_[layer] = log_c_[layer] / lambda;
    }
}

template <typename BigInt>
void HyperCIGAM<BigInt>::InitSizeWeights() {
    log_size_weight_.clear();

    const auto sizes = HyperedgeSizes();

    if (config_.size_decay <= 0.0 || config_.size_decay >= 1.0) {
        for (const SInt k: sizes) {
            log_size_weight_[k] = 0.0L;
        }
        return;
    }

    const SInt        k_min = sizes.front();
    const long double x     = static_cast<long double>(config_.size_decay);

    long double log_z = -std::numeric_limits<long double>::infinity();

    for (const SInt k: sizes) {
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

    for (const SInt k: HyperedgeSizes()) {
        LogBinomCache cache(k - 1);

        long double              log_sum = -std::numeric_limits<long double>::infinity();
        std::vector<long double> layer_mass(NumLayers());

        for (SInt layer = 0; layer < NumLayers(); ++layer) {
            layer_mass[layer] = LogExpectedRangeMass(k, 0, config_.n, layer, cache);
            log_sum           = LogAdd(log_sum, layer_mass[layer]);
        }

        log_mass_by_size_[k] = log_sum;
    }
}
template <typename BigInt>
long double HyperCIGAM<BigInt>::LogSizeWeight(const SInt k) const {
    return log_size_weight_.at(k);
}

template <typename BigInt>
void HyperCIGAM<BigInt>::GenerateCSR() {
    graph_.hyperedge_offsets.push_back(0);
    NormalizeParams(config_);
    InitLayerBounds();
    InitProbabilityConstants();
    InitSizeWeights();
    InitMassCache();
    const SInt n = config_.n;
#ifndef NDEBUG
    debug_edges_per_layer_.assign(NumLayers(), 0);
#endif

    if (config_.edge_budget > 0.0) {
        for (const SInt k: HyperedgeSizes()) {
            const long double log_Z_k = log_mass_by_size_.at(k);

            if (log_Z_k == -std::numeric_limits<long double>::infinity()) {
                throw ConfigurationError("CIGAM edge budget requested, but expected mass for size is zero");
            }

            log_edge_scaling_by_size_[k] =
                std::log(static_cast<long double>(config_.edge_budget)) + LogSizeWeight(k) - log_Z_k;
        }
    } else {
        for (const SInt k: HyperedgeSizes()) {
            log_edge_scaling_by_size_[k] = 0.0L;
        }
    }

    const SInt vertices_per_pe = n / size_;
    const SInt remainder       = n % size_;
    const SInt local_begin     = (rank_ * vertices_per_pe) + std::min<SInt>(rank_, remainder);
    const SInt local_end       = local_begin + vertices_per_pe + static_cast<SInt>(rank_ < remainder);

    SetVertexRange(local_begin, local_end);
    if (config_.approx || config_.n > (SInt{1} << 31)) {
        GenerateApproxCSR(local_begin, local_end);
    } else {
        GenerateExactCSR(local_begin, local_end);
    }
#ifndef NDEBUG
    std::cerr << "CIGAM layer statistics:\n";

    for (SInt layer = 0; layer < NumLayers(); ++layer) {
        std::cerr << "  layer " << layer << " [" << layer_begin_[layer] << ", " << layer_end_[layer] << ") -> "
                  << debug_edges_per_layer_[layer] << " edges\n";
    }
#endif
}

template <typename BigInt>
void HyperCIGAM<BigInt>::GenerateBoundedBlock(
    const SInt k, const SInt i, const SInt layer, LogBinomCache& binom_cache) {
    const auto [j_min, j_max] = LayerEndpointRange(i, layer);

    if (j_min > j_max || j_max - i < k - 1) {
        return;
    }

    const long double log_block_size = LogBlockSize(k, i, j_min, j_max, binom_cache);

    const long double log_p = LogProbabilityForDominant(i, layer) + log_edge_scaling_by_size_.at(k);

    const SInt seed = sampling::Spooky::hash(
        static_cast<unsigned long long>(config_.seed) + (static_cast<unsigned long long>(k) * kCountSeedMultiplier)
        + (static_cast<unsigned long long>(i) * kEdgeSeedMultiplier)
        + (static_cast<unsigned long long>(layer) * kRankSeedMultiplier));

    const SInt local_m = SampleBlockCountFromLogSize(log_block_size, log_p, rng_, seed, "CIGAM");

    if (local_m == 0) {
        return;
    }

    SInt edge_seed = sampling::Spooky::hash(seed + kEdgeRankSeedMultiplier);
    mersenne_.RandomInit(edge_seed);

    std::unordered_set<std::vector<SInt>, VectorHash> local_seen;

    if (!config_.allow_duplicates) {
        local_seen.max_load_factor(0.5);
        local_seen.reserve(static_cast<std::size_t>(local_m));
    }
    SInt           generated    = 0;
    CountInt       attempts     = 0;
    const CountInt max_attempts = std::max(CountInt(local_m) * 10, CountInt(1000));

    if (!config_.allow_duplicates) {
        const long double density = static_cast<long double>(local_m) / expl(log_block_size);

        if (density > 0.25L) {
            throw ConfigurationError("Dense exact CIGAM block not implemented yet; use --fast");
        }
    }

    while (generated < local_m) {
        const SInt j = SampleEndpoint(k, i, j_min, j_max, log_block_size, binom_cache);

        std::vector<SInt> pins;
        pins.reserve(k);

        pins.push_back(i);

        if (k > 2) {
            auto middle = FloydSample(i + 1, j - i - 1, k - 2, mersenne_);
            pins.insert(pins.end(), middle.begin(), middle.end());
        }

        pins.push_back(j);

        if (config_.allow_duplicates || local_seen.insert(pins).second) {
            AssertHyperedgeInvariants(pins, layer);
#ifndef NDEBUG
            ++debug_edges_per_layer_[layer];
#endif
            PushHyperedge(pins);
            ++generated;
        }

        if (++attempts > max_attempts) {
            throw ConfigurationError("CIGAM rejection sampling exceeded attempt limit; use --fast");
        }
    }
}

template <typename BigInt>
long double HyperCIGAM<BigInt>::RankValue(const SInt i) {
    const long double n = static_cast<long double>(config_.n);
    const long double u = 1.0L - ((static_cast<long double>(i) + 0.5L) / n);

    return InverseTruncatedExpCDF(u);
}

template <typename BigInt>
long double HyperCIGAM<BigInt>::InverseTruncatedExpCDF(const long double u) const {
    const long double lambda = static_cast<long double>(config_.cigam_lambda);
    return -std::log1pl(-u * (1.0L - expl(-lambda))) / lambda;
}

template <typename BigInt>
void HyperCIGAM<BigInt>::InitLayerBounds() {
    const SInt n = config_.n;
    const SInt L = NumLayers();

    layer_begin_.assign(L, 0);
    layer_end_.assign(L, 0);

    auto first_pos_with_x_greater_than = [&](const long double threshold) {
        SInt left  = 0;
        SInt right = n;

        while (left < right) {
            const SInt        mid = left + ((right - left) / 2);
            const long double x   = 1.0L - RankValue(mid);

            if (x > threshold) {
                right = mid;
            } else {
                left = mid + 1;
            }
        }

        return left;
    };

    for (SInt layer = 0; layer < L; ++layer) {
        layer_begin_[layer] = layer == 0 ? 0 : first_pos_with_x_greater_than(config_.cigam_breakpoints[layer - 1]);

        layer_end_[layer] = layer + 1 == L ? n : first_pos_with_x_greater_than(config_.cigam_breakpoints[layer]);
    }
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

template <typename BigInt>
SInt HyperCIGAM<BigInt>::EstimateEndpointInitialGuess(
    const SInt k, const SInt i, const SInt j_min, const SInt j_max, const long double log_target_prefix,
    LogBinomCache& binom_cache) {
    const SInt q = k - 1;

    const SInt low = j_min - i - 1;

    const long double log_low_mass = low >= q ? binom_cache.Get(low, q) : -std::numeric_limits<long double>::infinity();

    const long double log_abs_target = log_low_mass == -std::numeric_limits<long double>::infinity()
                                           ? log_target_prefix
                                           : LogAdd(log_low_mass, log_target_prefix);

    const long double log_q_factorial = std::lgammal(static_cast<long double>(q) + 1.0L);

    const long double x = expl((log_abs_target + log_q_factorial) / static_cast<long double>(q));

    SInt j = i + static_cast<SInt>(std::llround(x));

    return std::clamp<SInt>(j, j_min, j_max);
}

template <typename BigInt>
SInt HyperCIGAM<BigInt>::SampleEndpointBinarySearch(
    const SInt k, const SInt i, const SInt j_min, const SInt j_max, const long double log_target,
    LogBinomCache& binom_cache) {
    const SInt q = k - 1;

    const SInt        low     = j_min - i - 1;
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

template <typename BigInt>
SInt HyperCIGAM<BigInt>::SampleEndpoint(
    const SInt k, const SInt i, const SInt j_min, const SInt j_max, const long double log_block_size,
    LogBinomCache& binom_cache) {
    const long double u = std::min<long double>(
        static_cast<long double>(rng_.GenerateCanonicalDoubleStream()), std::nextafter(1.0L, 0.0L));

    if (k == 2) {
        return j_min + static_cast<SInt>(u * static_cast<long double>(j_max - j_min + 1));
    }

    const SInt q = k - 1;

    const long double log_target = std::log(u) + log_block_size;

    const SInt low = j_min - i - 1;

    const long double log_low = low >= q ? binom_cache.Get(low, q) : -std::numeric_limits<long double>::infinity();

    const bool use_linear_prefix = q <= 8 && (j_max - i) < (SInt{1} << 30);

    const long double low_mass = use_linear_prefix && low >= q ? BinomSmallExactLD(low, q) : 0.0L;

    auto prefix_ge_target = [&](const SInt j) {
        if (use_linear_prefix) {
            const long double prefix = BinomSmallExactLD(j - i, q) - low_mass;
            return std::log(prefix) >= log_target;
        }

        return LogDifferenceOfExponentials(binom_cache.Get(j - i, q), log_low) >= log_target;
    };

    auto prefix_lt_target = [&](const SInt j) {
        if (use_linear_prefix) {
            const long double prefix = BinomSmallExactLD(j - i, q) - low_mass;
            return std::log(prefix) < log_target;
        }

        return LogDifferenceOfExponentials(binom_cache.Get(j - i, q), log_low) < log_target;
    };

    SInt j = EstimateEndpointInitialGuess(k, i, j_min, j_max, log_target, binom_cache);

    constexpr SInt kMaxCorrectionSteps = 64;
    SInt           steps               = 0;

    while (j > j_min && prefix_ge_target(j - 1)) {
        --j;

        if (++steps > kMaxCorrectionSteps) {
            return SampleEndpointBinarySearch(k, i, j_min, j_max, log_target, binom_cache);
        }
    }

    while (j < j_max && prefix_lt_target(j)) {
        ++j;

        if (++steps > kMaxCorrectionSteps) {
            return SampleEndpointBinarySearch(k, i, j_min, j_max, log_target, binom_cache);
        }
    }

    return j;
}

template <typename BigInt>
long double HyperCIGAM<BigInt>::LogExpectedRangeMass(
    const SInt k, const SInt i_begin, const SInt i_end, const SInt layer, LogBinomCache& cache) {
    long double log_sum = -std::numeric_limits<long double>::infinity();

    for (SInt i = i_begin; i < i_end; ++i) {
        if (config_.n - i - 1 < k - 1) {
            continue;
        }

        const auto [j_min, j_max] = LayerEndpointRange(i, layer);

        if (j_min > j_max || j_max - i < k - 1) {
            continue;
        }

        const long double log_block_size = LogBlockSize(k, i, j_min, j_max, cache);
        const long double log_p          = LogProbabilityForDominant(i, layer);
        log_sum                          = LogAdd(log_sum, log_block_size + log_p);
    }
    return log_sum;
}

template <typename BigInt>
long double HyperCIGAM<BigInt>::LogProbabilityForDominant(const SInt i, const SInt layer) const {
    const long double n = static_cast<long double>(config_.n);
    const long double u = 1.0L - ((static_cast<long double>(i) + 0.5L) / n);

    return (-2.0L * log_c_[layer]) - (log_c_over_lambda_[layer] * std::log1pl(-u * lambda_exp_term_));
}

template <typename BigInt>
void HyperCIGAM<BigInt>::GenerateApproxRange(
    const SInt k, const SInt i_begin, const SInt i_end, const SInt layer, const SInt level, LogBinomCache& cache) {
    if (i_begin >= i_end) {
        return;
    }

    const long double log_expected =
        LogExpectedRangeMass(k, i_begin, i_end, layer, cache) + log_edge_scaling_by_size_.at(k);

    if (log_expected == -std::numeric_limits<long double>::infinity()) {
        return;
    }

    if (log_expected > std::log(static_cast<long double>(std::numeric_limits<SInt>::max()))) {
        throw ConfigurationError("CIGAM approximate range expected count exceeds SInt");
    }

    const long double expected = expl(log_expected);

    if (expected <= 0.0L) {
        return;
    }

    constexpr SInt kTargetEdgesPerBlock = 16384;
    constexpr SInt kTargetRangeWidth    = SInt{1} << 20;

    const SInt width = i_end - i_begin;

    const bool split_by_edges = expected > static_cast<long double>(kTargetEdgesPerBlock);

    const bool split_by_width = width > kTargetRangeWidth;

    if ((split_by_edges || split_by_width) && width > 1) {
        const SInt mid =
            split_by_width
                ? i_begin + (width / 2)
                : FindApproxMassSplit(k, i_begin, i_end, layer, log_expected - log_edge_scaling_by_size_.at(k), cache);

        GenerateApproxRange(k, i_begin, mid, layer, level + 1, cache);
        GenerateApproxRange(k, mid, i_end, layer, level + 1, cache);

        return;
    }
    const SInt count_seed = sampling::Spooky::hash(
        static_cast<unsigned long long>(config_.seed) + (static_cast<unsigned long long>(k) * kCountSeedMultiplier)
        + (static_cast<unsigned long long>(layer) * kRankSeedMultiplier)
        + (static_cast<unsigned long long>(i_begin) * kEdgeSeedMultiplier)
        + (static_cast<unsigned long long>(level) * kEdgeRankSeedMultiplier));

    const SInt local_m = rng_.GeneratePoisson(count_seed, static_cast<double>(expected));

    if (local_m == 0) {
        return;
    }

    std::vector<SInt>        candidates;
    std::vector<long double> cdf;

    const std::size_t reserve_size = static_cast<std::size_t>(std::min<SInt>(width, std::max<SInt>(1024, local_m * 4)));

    candidates.reserve(reserve_size);
    cdf.reserve(reserve_size);

    long double total = 0.0L;

    for (SInt i = i_begin; i < i_end; ++i) {
        if (config_.n - i - 1 < k - 1) {
            continue;
        }

        const auto [j_min, j_max] = LayerEndpointRange(i, layer);

        if (j_min > j_max || j_max - i < k - 1) {
            continue;
        }

        const long double log_block_size = LogBlockSize(k, i, j_min, j_max, cache);

        const long double log_p = LogProbabilityForDominant(i, layer);

        const long double weight = expl(log_block_size + log_p);

        if (weight <= 0.0L || !std::isfinite(static_cast<double>(weight))) {
            continue;
        }

        total += weight;
        candidates.push_back(i);
        cdf.push_back(total);
    }

    if (candidates.empty() || total <= 0.0L) {
        return;
    }

    for (auto& x: cdf) {
        x /= total;
    }
    cdf.back() = 1.0L;

    SInt edge_seed = sampling::Spooky::hash(static_cast<unsigned long long>(count_seed) + kEdgeRankSeedMultiplier);

    mersenne_.RandomInit(edge_seed);

    for (SInt e = 0; e < local_m; ++e) {
        const long double u = std::min<long double>(
            static_cast<long double>(rng_.GenerateCanonicalDoubleStream()), std::nextafter(1.0L, 0.0L));

        const auto it  = std::lower_bound(cdf.begin(), cdf.end(), u);
        const SInt idx = static_cast<SInt>(it - cdf.begin());
        const SInt i   = candidates[idx];

        const auto [j_min, j_max] = LayerEndpointRange(i, layer);

        const long double log_block_size = LogBlockSize(k, i, j_min, j_max, cache);

        const SInt j = SampleEndpoint(k, i, j_min, j_max, log_block_size, cache);

        std::vector<SInt> pins;
        pins.reserve(k);

        pins.push_back(i);

        if (k > 2) {
            auto middle = FloydSample(i + 1, j - i - 1, k - 2, mersenne_);
            pins.insert(pins.end(), middle.begin(), middle.end());
        }

        pins.push_back(j);
        AssertHyperedgeInvariants(pins, layer);
#ifndef NDEBUG
        ++debug_edges_per_layer_[layer];
#endif
        PushHyperedge(pins);
    }
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
HyperCIGAM<BigInt>::HyperCIGAM(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size),
      vertex_permutation_(
          random_permutation::FeistelPseudoRandomPermutation::buildPermutation(
              static_cast<std::uint64_t>(config.n - 1), static_cast<std::uint64_t>(config.seed))),
      rng_(config) {}

template <typename BigInt>
void HyperCIGAM<BigInt>::FinalizeCSR(MPI_Comm /*comm*/) {
    for (SInt& pin: graph_.hyperedge_pins) {
        pin = static_cast<SInt>(vertex_permutation_.f(static_cast<std::uint64_t>(pin)));
    }
}

template <typename BigInt>
void HyperCIGAM<BigInt>::PushHyperedge(const std::vector<SInt>& pins) {
    graph_.hyperedge_pins.insert(graph_.hyperedge_pins.end(), pins.begin(), pins.end());
    graph_.hyperedge_offsets.push_back(static_cast<SInt>(graph_.hyperedge_pins.size()));
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
template class HyperCIGAM<__int128>;
template class HyperCIGAM<SInt>;
#pragma GCC diagnostic pop

} // namespace kagen
