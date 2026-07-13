
#include "kagen/generators/generator.h"
#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"
#include "kagen/kagen.h"
#include "kagen/tools/random_permutation.h"
#include "kagen/tools/rng_wrapper.h"

#include <optional>

namespace kagen {
struct ApproxRangeStatsCIGAM {
    long double log_expected;
    long double expected;
};
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
    void                  InitQuantileGenerationState();
    std::pair<SInt, SInt> ComputeLocalVertexRange();
    void                  InitEdgeBudgetScaling();
    SInt                  EstimateEndpointInitialGuess(
        SInt k, SInt i, SInt j_min, SInt j_max, long double log_target_prefix, LogBinomCache& binom_cache);
    SInt BlockCountSeed(SInt k, SInt i, SInt layer);
    SInt EdgeSeed(SInt count_seed);

    void              ValidateExactBlockDensity(SInt local_m, long double log_block_size);
    std::vector<SInt> SampleBoundedHyperedge(
        SInt k, SInt i, SInt j_min, SInt j_max, long double log_block_size, LogBinomCache& binom_cache);
    void GenerateApproxLeafRange(
        SInt k, SInt i_begin, SInt i_end, SInt layer, SInt level, long double expected, LogBinomCache& cache,
        const std::vector<double>& prefix);

    SInt ApproxRangeCountSeed(SInt k, SInt i_begin, SInt layer, SInt level) const;
    bool ShouldSplitApproxRange(SInt width, long double expected) const;
    SInt ChooseApproxRangeSplit(SInt i_begin, SInt i_end, const std::vector<double>& prefix);

    std::optional<ApproxRangeStatsCIGAM> ComputeApproxRangeStats(
        SInt k, SInt i_begin, SInt i_end, SInt layer, LogBinomCache& cache, const std::vector<double>& prefix) const;
    void PushCheckedHyperedge(const std::vector<SInt>& pins, SInt layer);
    void GenerateBoundedBlock(SInt k, SInt i, SInt layer, LogBinomCache& binom_cache);
    SInt SampleEndpoint(SInt k, SInt i, SInt j_min, SInt j_max, long double log_block_size, LogBinomCache& binom_cache);
    SInt SampleEndpointBinarySearch(
        SInt k, SInt i, SInt j_min, SInt j_max, long double log_target, LogBinomCache& binom_cache);
    SInt                  SampleDominantVertexFromPrefix(SInt i_begin, SInt i_end, const std::vector<double>& prefix);
    long double           LogBlockSize(SInt k, SInt i, SInt j_min, SInt j_max, LogBinomCache& binom_cache) const;
    void                  InitProbabilityConstants();
    std::pair<SInt, SInt> LayerEndpointRange(SInt i, SInt layer) const;
    std::vector<SInt>     HyperedgeSizes() const;
    std::vector<double>   BuildDominantMassPrefix(SInt k, SInt layer) const;
    SInt                  FindDominantBoundaryInPrefix(const std::vector<double>& prefix, PEID rank, PEID size) const;
    void                  GenerateApproxCSR();
    void                  GenerateExactCSR();
    void                  GenerateApproxRange(
        SInt k, SInt i_begin, SInt i_end, SInt layer, SInt level, LogBinomCache& cache,
        const std::vector<double>& prefix);

    SInt NumLayers() const {
        return static_cast<SInt>(config_.cigam_c.size());
    }

    std::size_t PrefixIndex(const std::size_t k_idx, const SInt layer) const {
        return (k_idx * NumLayers()) + layer;
    }

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

    HypergraphMemoryStats memory_stats_;

    void                                  InitMassCache();
    void                                  InitSizeWeights();
    long double                           LogSizeWeight(SInt k) const;
    std::unordered_map<SInt, long double> log_edge_scaling_by_size_;
    SInt FindApproxMassSplit(SInt i_begin, SInt i_end, const std::vector<double>& prefix);

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

    PEID PythonBlockOwner(SInt k, SInt dominant, SInt layer) const;

    SInt PythonBlockCountSeed(SInt k, SInt dominant, SInt layer) const;

    SInt PythonBlockEdgeSeed(SInt k, SInt dominant, SInt layer) const;

    long double CounterUniform01(std::uint64_t stream_tag, std::uint64_t index) const;

    CountInt PythonBlockPopulation(SInt k, SInt dominant, SInt j_min, SInt j_max) const;

    std::vector<SInt> SamplePythonBlockHyperedge(SInt k, SInt dominant, SInt j_min, SInt j_max);

    void
    GeneratePythonBlock(SInt k, SInt dominant, SInt layer, SInt j_min, SInt j_max, const PythonReferenceRanks& ranks);

    CountInt ExactBlockPopulation(SInt k, SInt i, SInt j_min, SInt j_max) const;

    bool pins_are_final_vertex_ids_ = false;

#ifndef NDEBUG
    std::vector<SInt> debug_edges_per_layer_;
#endif

    void AssertHyperedgeInvariants(const std::vector<SInt>& pins, SInt layer) const;
};

using HyperCIGAMBig   = HyperCIGAM<__int128>;
using HyperCIGAMSmall = HyperCIGAM<SInt>;

} // namespace kagen