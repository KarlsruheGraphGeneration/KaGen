#pragma once

#include "kagen/generators/generator.h"
#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"
#include "kagen/hypergraph/debug_logger_erdos.h"
#include "kagen/kagen.h"
#include "kagen/tools/rng_wrapper.h"

#include <cstddef>
#include <limits>
#include <optional>
#include <vector>

namespace kagen {

enum class ProbabilityMode : std::uint8_t {
    ExplicitProbabilities,
    ExplicitExpectedCounts,
    EdgeBudget,
    EdgeAndPinBudget,
    GlobalProbability
};

struct SizeGenerationParameters {
    double probability    = 0.0;
    double expected_count = -1.0;
};

struct HGNPLocalGenerationRange {
    SInt begin = 0;
    SInt end   = 0;

    // Number of hyperedges generated locally for this size.
    SInt local_m = 0;
};

struct HGNPSizePlan {
    SInt                     hyperedge_size = 0;
    SInt                     partition_id   = 0;
    HGNPLocalGenerationRange range;
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
    SInt                      SampleExactEdgeCountNative(SInt population, double probability, SInt seed);
    SInt                      SampleExactEdgeCountHuge(const CountInt& exact_population, double probability, SInt seed);
    std::vector<HGNPSizePlan> BuildGenerationPlan(SInt lower_bound, SInt upper_bound);
    SInt                      SampleExactEdgeCount(const CountInt& population, double probability, SInt seed);
    bool         AppendSizePlanIfNeeded(SInt hyperedge_size, SInt lower_bound, std::vector<HGNPSizePlan>& plan);
    HGNPSizePlan PrepareSizePlan(SInt hyperedge_size, double probability, SInt partition_id);
    HGNPLocalGenerationRange PrepareApproxLocalRange(SInt hyperedge_size, double probability, SInt partition_id);
    SInt                     LocalCountSeed(SInt hyperedge_size, SInt partition_id) const;
    void                     PrepareSampledExactPlan(HGNPSizePlan& entry, double probability);
    long double
    LogBinomialPoissonRatioRelativeToMode(SInt value, SInt mode, long double population, long double probability) const;
    void                     ReserveCSRForPlan(const std::vector<HGNPSizePlan>& plan);
    void                     GenerateHyperedgesFromPlan(const std::vector<HGNPSizePlan>& plan);
    void                     GenerateHyperedgesFromSizePlan(const HGNPSizePlan& entry);
    void                     LoadSizeDistributionInputs();
    void                     LoadSizeProbabilitiesFromFile();
    void                     LoadExpectedSizeCountsFromFile();
    void                     SelectProbabilityMode();
    SizeGenerationParameters GetSizeGenerationParameters(SInt hyperedge_size, SInt lower_bound);
    std::optional<double> ResolveProbabilityForSize(SInt hyperedge_size, const SizeGenerationParameters& params) const;
    void                  GenerateBudgetSizeCounts(SInt lower_bound, SInt upper_bound);
    SInt                  EffectiveUpperBound(SInt hard_upper_bound, SInt lower_bound);
    bool                  UsesExpectedCount() const;
    bool                  ExpectedCountSkipsRemainingSizes(const SizeGenerationParameters& params) const;
    bool                  ExpectedCountSkipsCurrentSize(const SizeGenerationParameters& params) const;
    void                  ValidateProbability(double probability) const;
    bool                  ShouldSkipSizeGeneration(SInt hyperedge_size, double probability) const;
    void                  SetLocalVertexRange();
    std::pair<SInt, SInt> LocalMinOwnerRange(SInt hyperedge_size, SInt partition_id) const;
    void                  GenerateLocalHyperedges(
        const HGNPSizePlan& entry, SInt& edge_seed, LogBinomCache& cache, std::vector<SInt>& pins,
        HyperedgeSeenSet& seen);
    void SampleLocalHyperedgeInto(
        SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, SInt& edge_seed, LogBinomCache& log_binom_cache,
        std::vector<SInt>& pins, std::uint64_t& minimum_search_steps, std::uint64_t& minimum_cache_gets);
    SInt             LocalEdgeSeed(SInt hyperedge_size, SInt partition_id) const;
    void             GenerateApproxHyperedgesFromPlan(const HGNPSizePlan& entry);
    PGeneratorConfig config_;

    PEID                  rank_;
    PEID                  size_;
    HypergraphMemoryStats memory_stats_;
    FloydScratchSet       floyd_scratch_;

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    mutable HGNPInstrumentation instrumentation_;
#endif
    RNGWrapper<>                   rng_;
    Mersenne                       mersenne_;
    ProbabilityMode                probs_type_ = ProbabilityMode::GlobalProbability;
    std::unordered_map<SInt, SInt> budget_size_counts_;

    std::optional<ErdosHypergraphDebugLogger> debug_logger_;
    SInt                                      next_debug_hyperedge_id_ = 0;
};
using HyperGNPSmall = HyperGNP<SInt>;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
using HyperGNPBig = HyperGNP<__int128>;
#pragma GCC diagnostic pop

} // namespace kagen