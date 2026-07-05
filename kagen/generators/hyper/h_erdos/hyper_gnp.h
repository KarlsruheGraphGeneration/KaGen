#pragma once

#include "kagen/generators/generator.h"
#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"
#include "kagen/hypergraph/debug_logger_erdos.h"
#include "kagen/kagen.h"
#include "kagen/tools/rng_wrapper.h"

#include <optional>
#include <vector>

namespace kagen {

enum class ProbabilityMode : std::uint8_t {
    ExplicitProbabilities,
    ExplicitExpectedCounts,
    EdgeBudget,
    GlobalProbability
};

struct SizeGenerationParameters {
    double probability    = 0.0;
    double expected_count = -1.0;
};
struct ApproxRangeStats {
    long double local_mass = 0.0L;
    long double log_total  = 0.0L;
    long double expected   = 0.0L;
};

class HyperGNPFactory : public GeneratorFactory {
public:
    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const override;
};

template <typename BigInt>
class HyperGNP final : public virtual Generator, private CSROnlyGenerator {
public:
    HyperGNP(const PGeneratorConfig& config, PEID rank, PEID size);

    void GenerateCSR() final;

    void FinalizeCSR(MPI_Comm comm) final;

private:
    void GenerateHyperedgesForAllSizes(SInt lower_bound, SInt upper_bound);
    bool GenerateHyperedgesForSizeIfNeeded(SInt hyperedge_size, SInt lower_bound);
    void GenerateHyperedgesOfSize(SInt hyperedge_size, double p);

    void                     LoadSizeDistributionInputs();
    void                     LoadSizeProbabilitiesFromFile();
    void                     LoadExpectedSizeCountsFromFile();
    void                     SelectProbabilityMode();
    SizeGenerationParameters GetSizeGenerationParameters(SInt hyperedge_size, SInt lower_bound);
    std::optional<double> ResolveProbabilityForSize(SInt hyperedge_size, const SizeGenerationParameters& params) const;
    SInt                  EffectiveUpperBound(SInt hard_upper_bound, SInt lower_bound);
    bool                  UsesExpectedCount() const;
    bool                  ExpectedCountSkipsRemainingSizes(const SizeGenerationParameters& params) const;
    bool                  ExpectedCountSkipsCurrentSize(const SizeGenerationParameters& params) const;
    void                  ValidateProbability(double probability) const;
    bool                  ShouldSkipSizeGeneration(SInt hyperedge_size, double probability) const;
    bool                  ShouldUseApproxGeneration(SInt hyperedge_size) const;

    void                  SetLocalVertexRange();
    std::pair<SInt, SInt> LocalMinOwnerRange(SInt hyperedge_size) const;
    bool                  ShouldSkipLocalRange(SInt local_min_begin, SInt local_min_end, double probability) const;

    void GenerateHyperedgesOfSizeApprox(SInt hyperedge_size, double p);
    void GenerateApproxLocalRange(
        SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, double probability, SInt level,
        LogBinomCache& count_cache);
    std::optional<ApproxRangeStats> ComputeApproxRangeStats(
        SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, double probability,
        LogBinomCache& count_cache) const;
    bool ShouldSplitApproxRange(SInt local_min_begin, SInt local_min_end, long double expected) const;
    SInt ChooseApproxRangeSplit(
        SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, long double local_mass,
        LogBinomCache& count_cache);
    SInt FindApproxMassSplit(SInt begin, SInt end, SInt n, SInt k, long double total_mass, LogBinomCache& cache);
    void GenerateApproxLeafRange(
        SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, double probability, SInt level,
        const ApproxRangeStats& stats);
    void GenerateLocalHyperedges(
        SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, SInt local_edge_count,
        LogBinomCache& log_binom_cache);
    std::vector<SInt> SampleLocalHyperedge(
        SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, SInt& edge_seed, LogBinomCache& log_binom_cache);
    SInt LocalEdgeSeed(SInt hyperedge_size) const;

    void GenerateHyperedgesOfSizeExact(SInt hyperedge_size, double p);
    void GenerateExactOwnedPrefixRange(
        SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, double probability, SInt level);
    void GenerateExactPrefixRange(
        SInt hyperedge_size, std::vector<SInt> prefix, SInt candidate_begin, SInt candidate_end, SInt remaining_pins,
        double probability, SInt level);
    CountInt PrefixRangePopulation(SInt candidate_begin, SInt candidate_end, SInt remaining_pins) const;
    SInt
    FindExactPrefixPopulationSplit(SInt begin, SInt end, SInt remaining_pins, const CountInt& total_population) const;
    void GenerateExactPrefixLeafRange(
        SInt hyperedge_size, std::vector<SInt> prefix, SInt candidate_begin, SInt candidate_end, SInt remaining_pins,
        double probability, SInt level, const CountInt& population);
    void GenerateSampledPrefixHyperedges(
        const std::vector<SInt>& prefix, SInt candidate_begin, SInt candidate_end, SInt remaining_pins, SInt local_m,
        Mersenne& mersenne);
    std::vector<SInt> SamplePrefixHyperedge(
        const std::vector<SInt>& prefix, SInt candidate_begin, SInt candidate_end, SInt remaining_pins,
        Mersenne& mersenne);
    SInt
    SampleNextPrefixVertex(SInt candidate_begin, SInt candidate_end, SInt remaining_pins, Mersenne& mersenne) const;
    void ValidateExactSparseDensity(SInt local_m, const CountInt& population) const;
    void ValidateExactMinOwnerPartition(SInt hyperedge_size);

    void                                              ValidateDuplicateCheckingIsFeasible(SInt local_edge_count) const;
    std::unordered_set<std::vector<SInt>, VectorHash> MakeLocalSeenSet(SInt local_edge_count) const;
    SInt CalculateCountSeed(SInt hyperedge_size, SInt candidate_begin, SInt level) const;

    void LogSizeProbability(SInt hyperedge_size, double probability);
    void PushHyperedge(const std::vector<SInt>& pins);

    PGeneratorConfig config_;

    PEID rank_;
    PEID size_;

    RNGWrapper<>    rng_;
    Mersenne        mersenne_;
    ProbabilityMode probs_type_ = ProbabilityMode::GlobalProbability;

    std::optional<ErdosHypergraphDebugLogger> debug_logger_;
};
using HyperGNPSmall = HyperGNP<SInt>;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
using HyperGNPBig = HyperGNP<__int128>;
#pragma GCC diagnostic pop

} // namespace kagen