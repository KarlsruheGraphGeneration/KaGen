#pragma once

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"
#include "kagen/hypergraph/debug_logger_erdos.h"
#include "kagen/kagen.h"
#include "kagen/tools/rng_wrapper.h"

#include <memory>
#include <optional>
#include <vector>

namespace kagen {

class HyperGNMFactory : public GeneratorFactory {
public:
    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const override;
};
void ValidateHGNMConfig(const PGeneratorConfig& config);
struct HGNMLocalGenerationRange {
    SInt min_begin;
    SInt min_end;
    SInt local_m;
    bool use_approx;
};

struct HGNMSizePlan {
    SInt                     hyperedge_size;
    SInt                     m_k;
    HGNMLocalGenerationRange range;
};

template <typename BigInt>
class HyperGNM : public virtual Generator, private CSROnlyGenerator {
public:
    HyperGNM(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void GenerateCSR() final;
    void FinalizeCSR(MPI_Comm comm) final;

private:
    void GenerateHyperedgesFromPlan(const std::vector<HGNMSizePlan>& plan);

    void ValidateDuplicateCheckingIsFeasible(SInt local_m) const;
    void ValidateExactSparseDensity(SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, SInt local_m) const;
    HGNMLocalGenerationRange PrepareLocalGenerationRange(SInt hyperedge_size, SInt m_k);

    std::pair<SInt, SInt> ResolveSizeRange() const;
    double                ValidateAndGetSizeAlpha() const;
    void                  GenerateSizeCountsForConfig(
        SInt lower_bound, SInt upper_bound, double alpha, std::unordered_map<SInt, SInt>& size_counts);
    void                  LogSizeCounts(const std::unordered_map<SInt, SInt>& size_counts) const;
    std::pair<SInt, SInt> ComputeLocalVertexRange() const;
    SInt                  EdgeSeed(SInt hyperedge_size) const;

    void GenerateSampledHyperedges(
        SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, SInt local_m, SInt& edge_seed,
        LogBinomCache& log_binom_cache);

    void GenerateHyperedgesOfSize(SInt hyperedge_size, SInt m_k, const HGNMLocalGenerationRange& data);

    bool TryPushHyperedge(const std::vector<SInt>& pins, HyperedgeSeenSet& local_seen);

    void SampleHyperedgeInto(SInt minimum_vertex, SInt hyperedge_size, std::vector<SInt>& pins);

    SInt ApproximateLocalHyperedgeCount(SInt hyperedge_size, SInt m_k);

    SInt ApproximateLocalHyperedgeCountRecursive(
        SInt hyperedge_size, SInt rank_begin, SInt rank_end, SInt target_rank, long double population_mass, SInt draws,
        SInt level, LogBinomCache& log_binom_cache);

    void GenerateRandomGeometricSizeCounts(
        SInt lower, SInt upper, double alpha, std::unordered_map<SInt, SInt>& size_counts);

    void GenerateDeterministicGeometricSizeCounts(
        SInt lower, SInt upper, double decay, std::unordered_map<SInt, SInt>& size_counts);

    SInt RankSplitSeed(SInt hyperedge_size, SInt rank_begin, SInt rank_end, SInt level) const;

    void GenerateBoltzmannPinBudgetSizeCounts(
        SInt lower_bound, SInt upper_bound, SInt pin_budget, std::unordered_map<SInt, SInt>& size_counts);

    void AccumulateCacheStats(const LogBinomCache& cache, std::size_t requested_capacity);

    const PGeneratorConfig& config_;
    PEID                    rank_;
    PEID                    size_;
    HypergraphMemoryStats   memory_stats_;
    FloydScratchSet         floyd_scratch_;
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    mutable HGNPInstrumentation instrumentation_;
#endif

    RNGWrapper<>                              rng_;
    Mersenne                                  mersenne_;
    std::optional<ErdosHypergraphDebugLogger> debug_logger_;
};

using HyperGNMSmall = HyperGNM<SInt>;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
using HyperGNMBig = HyperGNM<__int128>;
#pragma GCC diagnostic pop

} // namespace kagen