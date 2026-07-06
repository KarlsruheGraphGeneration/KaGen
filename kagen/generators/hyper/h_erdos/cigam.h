
#include "kagen/generators/generator.h"
#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"
#include "kagen/hypergraph/debug_logger_erdos.h"
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
    void                  InitGenerationState();
    std::pair<SInt, SInt> ComputeLocalVertexRange();
    void                  InitEdgeBudgetScaling();

    SInt BlockCountSeed(SInt k, SInt i, SInt layer);
    SInt EdgeSeed(SInt count_seed);

    std::unordered_set<std::vector<SInt>, VectorHash> MakeLocalSeenSet(SInt local_edge_count) const;
    void              ValidateExactBlockDensity(SInt local_m, long double log_block_size);
    std::vector<SInt> SampleBoundedHyperedge(
        SInt k, SInt i, SInt j_min, SInt j_max, long double log_block_size, LogBinomCache& binom_cache);
    void GenerateApproxLeafRange(
        SInt k, SInt i_begin, SInt i_end, SInt layer, SInt level, long double expected, LogBinomCache& cache);

    SInt ApproxRangeCountSeed(SInt k, SInt i_begin, SInt layer, SInt level) const;
    bool ShouldSplitApproxRange(SInt width, long double expected) const;
    SInt ChooseApproxRangeSplit(SInt i_begin, SInt i_end, const std::vector<double>& prefix);
    bool BuildDominantVertexCDF(
        SInt k, SInt i_begin, SInt i_end, SInt layer, SInt local_m, LogBinomCache& cache, std::vector<SInt>& candidates,
        std::vector<long double>& cdf);

    std::optional<ApproxRangeStatsCIGAM> ComputeApproxRangeStats(
        SInt k, SInt i_begin, SInt i_end, SInt layer, LogBinomCache& cache, const std::vector<double>& prefix) const;
    void PushCheckedHyperedge(const std::vector<SInt>& pins, SInt layer);
    void GenerateBoundedBlock(SInt k, SInt i, SInt layer, LogBinomCache& binom_cache);
    SInt SampleEndpoint(SInt k, SInt i, SInt j_min, SInt j_max, long double log_block_size, LogBinomCache& binom_cache);
    SInt SampleEndpointBinarySearch(
        SInt k, SInt i, SInt j_min, SInt j_max, long double log_target, LogBinomCache& binom_cache);
    SInt        SampleDominantVertex(const std::vector<SInt>& candidates, const std::vector<long double>& cdf);
    long double LogBlockSize(SInt k, SInt i, SInt j_min, SInt j_max, LogBinomCache& binom_cache) const;
    long double RankValue(SInt i) const;
    long double InverseTruncatedExpCDF(long double u) const;
    void        InitLayerBounds();
    void        InitProbabilityConstants();
    std::pair<SInt, SInt> LayerEndpointRange(SInt i, SInt layer) const;
    void                  PushHyperedge(const std::vector<SInt>& pins);
    std::vector<SInt>     HyperedgeSizes() const;
    std::pair<SInt, SInt> LocalDominantOwnerRange(SInt k, SInt layer, LogBinomCache& cache) const;
    std::vector<double>   BuildDominantMassPrefix(SInt k, SInt layer, LogBinomCache& cache) const;
    SInt                  FindDominantBoundaryInPrefix(const std::vector<double>& prefix, PEID rank, PEID size) const;
    void                  GenerateApproxCSR();
    void                  GenerateExactCSR();
    void                  GenerateApproxRange(
        SInt k, SInt i_begin, SInt i_end, SInt layer, SInt level, LogBinomCache& cache,
        const std::vector<double>& prefix);
    long double LogExpectedRangeMass(SInt k, SInt i_begin, SInt i_end, SInt layer, LogBinomCache& cache) const;
    SInt        EstimateEndpointInitialGuess(
        SInt k, SInt i, SInt j_min, SInt j_max, long double log_target_prefix, LogBinomCache& cache);
    long double LogProbabilityForDominant(SInt i, SInt layer) const;

    SInt NumLayers() const {
        return static_cast<SInt>(config_.cigam_c.size());
    }

    std::size_t PrefixIndex(const std::size_t k_idx, const SInt layer) const {
        return (k_idx * NumLayers()) + layer;
    }

    PGeneratorConfig                                   config_;
    PEID                                               rank_;
    PEID                                               size_;
    std::vector<SInt>                                  layer_begin_;
    std::vector<SInt>                                  layer_end_;
    random_permutation::FeistelPseudoRandomPermutation vertex_permutation_;

    RNGWrapper<>                              rng_;
    Mersenne                                  mersenne_;
    Mersenne                                  count_mersenne_;
    std::optional<ErdosHypergraphDebugLogger> debug_logger_;

    // Caching
    std::vector<long double>              log_c_;
    std::vector<long double>              log_c_over_lambda_;
    long double                           lambda_exp_term_ = 0.0L;
    std::unordered_map<SInt, long double> log_size_weight_;
    std::unordered_map<SInt, long double> log_mass_by_size_;
    std::vector<std::vector<double>>      dominant_mass_prefix_;
    std::vector<SInt>                     cigam_sizes_;

    void                                  InitMassCache();
    void                                  InitSizeWeights();
    long double                           LogSizeWeight(SInt k) const;
    std::unordered_map<SInt, long double> log_edge_scaling_by_size_;
    SInt FindApproxMassSplit(SInt i_begin, SInt i_end, const std::vector<double>& prefix);

#ifndef NDEBUG
    std::vector<SInt> debug_edges_per_layer_;
#endif

    void AssertHyperedgeInvariants(const std::vector<SInt>& pins, SInt layer) const;
};

using HyperCIGAMBig   = HyperCIGAM<__int128>;
using HyperCIGAMSmall = HyperCIGAM<SInt>;

} // namespace kagen