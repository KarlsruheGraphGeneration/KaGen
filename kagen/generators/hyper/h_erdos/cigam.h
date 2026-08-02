
#include "kagen/generators/generator.h"
#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"
#include "kagen/hypergraph/debug_logger_cigam.h"
#include "kagen/kagen.h"
#include "kagen/tools/random_permutation.h"
#include "kagen/tools/rng_wrapper.h"

namespace kagen {
class HyperCIGAMFactory : public GeneratorFactory {
public:
    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const override;
};

template <typename BigInt>
class HyperCIGAM : public virtual Generator, private CSROnlyGenerator {
public:
    HyperCIGAM(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void GenerateCSR() final;
    void FinalizeCSR(MPI_Comm comm) final;

private:
    struct BlockInfo {
        SInt        layer;
        SInt        j_min;
        SInt        j_max;
        long double log_block_size;
    };

    void                  InitQuantileGenerationState();
    std::pair<SInt, SInt> ComputeLocalVertexRange();
    void                  InitEdgeBudgetScaling();
    SInt                  EstimateEndpointInitialGuess(
        SInt k, SInt i, SInt j_min, SInt j_max, long double log_target_prefix, LogBinomCache& binom_cache);
    SInt BlockCountSeed(SInt k, SInt i, SInt layer);
    SInt EdgeSeed(SInt count_seed);

    void ValidateExactBlockDensity(SInt local_m, long double log_block_size);
    void SampleBoundedHyperedge(
        SInt k, SInt dominant, SInt j_min, SInt j_max, long double log_block_size, LogBinomCache& binom_cache,
        std::vector<SInt>& pins, FloydScratchSet& scratch);
    void PushCheckedHyperedge(const std::vector<SInt>& pins, SInt layer);
    void GenerateBoundedBlock(SInt k, SInt i, SInt layer, LogBinomCache& binom_cache);
    SInt SampleEndpoint(SInt k, SInt i, SInt j_min, SInt j_max, long double log_block_size, LogBinomCache& binom_cache);
    SInt SampleEndpointBinarySearch(
        SInt k, SInt i, SInt j_min, SInt j_max, long double log_target, LogBinomCache& binom_cache);
    long double           LogBlockSize(SInt k, SInt i, SInt j_min, SInt j_max, LogBinomCache& binom_cache) const;
    void                  InitProbabilityConstants();
    std::pair<SInt, SInt> LayerEndpointRange(SInt i, SInt layer) const;
    std::vector<SInt>     HyperedgeSizes() const;
    std::vector<double>   BuildDominantMassPrefix(SInt k, SInt layer) const;
    SInt                  FindDominantBoundaryInPrefix(const std::vector<double>& prefix, PEID rank, PEID size) const;
    void                  GenerateApproxCSR();
    void                  GenerateExactCSR();

    SInt NumLayers() const {
        return static_cast<SInt>(config_.cigam_c.size());
    }

    std::size_t PrefixIndex(const std::size_t k_idx, const SInt layer) const {
        return (k_idx * NumLayers()) + layer;
    }

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    void AccumulateLogBinomStats(const LogBinomCacheStats& stats);
#endif

    // Position i is zero-based and ordered by decreasing prestige.
    //
    // The deterministic prestige is the midpoint quantile
    //
    //   u_i = 1 - (i + 1/2) / n,
    //   r_i = F^{-1}(u_i),
    //
    // where F is the truncated exponential CDF.
    long double RankQuantile(SInt position) const;
    long double RankValue(SInt position) const;

    // Numerically stable inverse CDF of the exponential distribution
    // truncated to [0, 1].
    long double InverseTruncatedExpCDF(long double u) const;

    // First sorted position i for which
    //
    //   1 - RankValue(i) > breakpoint.
    //
    // Returned value lies in [0, n].
    SInt FirstPositionAfterBreakpoint(long double breakpoint) const;

    void        InitLayerBounds();
    long double LogProbabilityForDominant(SInt i, SInt layer) const;

    PGeneratorConfig                                   config_;
    PEID                                               rank_;
    PEID                                               size_;
    std::vector<SInt>                                  layer_begin_;
    std::vector<SInt>                                  layer_end_;
    random_permutation::FeistelPseudoRandomPermutation vertex_permutation_;

    RNGWrapper<> rng_;
    Mersenne     mersenne_;

    // Caching
    std::vector<long double>              log_c_;
    std::unordered_map<SInt, long double> log_size_weight_;
    std::unordered_map<SInt, long double> log_mass_by_size_;
    std::vector<std::vector<double>>      dominant_mass_prefix_;
    std::vector<SInt>                     cigam_sizes_;
    std::vector<long double>              dominant_exponent_;

    std::optional<CIGAMHypergraphDebugLogger> debug_logger_;
    SInt                                      next_debug_hyperedge_id_ = 0;

    HypergraphMemoryStats memory_stats_;

    void                                  InitMassCache();
    void                                  InitSizeWeights();
    long double                           LogSizeWeight(SInt k) const;
    std::unordered_map<SInt, long double> log_edge_scaling_by_size_;

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION

    LogBinomCacheStats log_binom_stats_;
#endif
    // #### Paper-accurate implementation ####

    struct PythonReferenceRanks {
        // Sorted position -> original vertex ID.
        std::vector<SInt> original_vertex;

        // Sorted position -> sampled prestige, decreasing.
        std::vector<long double> rank;

        // Half-open sorted-position interval for each layer.
        std::vector<SInt> layer_begin;
        std::vector<SInt> layer_end;
    };

    void GeneratePythonDistributedCSR();

    PythonReferenceRanks GeneratePythonReferenceRanks() const;

    PEID BlockOwner(SInt k, SInt dominant, SInt layer) const;

    SInt PythonBlockCountSeed(SInt k, SInt dominant, SInt layer) const;

    SInt PythonBlockEdgeSeed(SInt k, SInt dominant, SInt layer) const;

    long double CounterUniform01(std::uint64_t stream_tag, std::uint64_t index) const;

    CountInt PythonBlockPopulation(SInt k, SInt dominant, SInt j_min, SInt j_max) const;

    std::vector<SInt> SamplePythonBlockHyperedge(
        SInt k, SInt dominant, SInt j_min, SInt j_max, SInt& sampling_attempts, SInt& block_rejections);

    void
    GeneratePythonBlock(SInt k, SInt dominant, SInt layer, SInt j_min, SInt j_max, const PythonReferenceRanks& ranks);

    CountInt ExactBlockPopulation(SInt k, SInt i, SInt j_min, SInt j_max) const;

    void GenerateApproxBlock(
        SInt k, SInt dominant, SInt layer, SInt j_min, SInt j_max, SInt log_block_size, LogBinomCache& cache);

    bool pins_are_final_vertex_ids_ = false;

#ifndef NDEBUG
    std::vector<SInt> debug_edges_per_layer_;
#endif

    void AssertHyperedgeInvariants(const std::vector<SInt>& pins, SInt layer) const;
};

using HyperCIGAMBig   = HyperCIGAM<__int128>;
using HyperCIGAMSmall = HyperCIGAM<SInt>;

} // namespace kagen