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

template <typename BigInt>
class HyperGNM : public virtual Generator, private CSROnlyGenerator {
public:
    HyperGNM(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void GenerateCSR() final;
    void FinalizeCSR(MPI_Comm comm) final;

private:
    void GenerateHyperedgesFromSizeCounts(const std::unordered_map<SInt, SInt>& size_counts);

    std::unordered_set<std::vector<SInt>, VectorHash> MakeLocalSeenSet(SInt local_m) const;
    void                                              ValidateDuplicateCheckingIsFeasible(SInt local_m) const;
    SInt ComputeLocalHyperedgeCount(SInt hyperedge_size, SInt m_k, bool use_approx);
    void ValidateExactSparseDensity(SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, SInt local_m) const;
    std::size_t ComputeCacheSize(SInt local_m, SInt local_min_begin, SInt local_min_end) const;

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

    void GenerateHyperedgesOfSize(SInt hyperedge_size, SInt m_k);

    bool TryPushHyperedge(const std::vector<SInt>& pins, std::unordered_set<std::vector<SInt>, VectorHash>& local_seen);

    std::vector<SInt> SampleHyperedge(SInt minimum_vertex, SInt hyperedge_size);

    void PushHyperedge(const std::vector<SInt>& pins);

    SInt ApproximateLocalHyperedgeCount(SInt hyperedge_size, SInt m_k);

    SInt ApproximateLocalHyperedgeCountRecursive(
        SInt hyperedge_size, SInt rank_begin, SInt rank_end, SInt target_rank, long double population_mass, SInt draws,
        SInt level, LogBinomCache& log_binom_cache);

    SInt ExactLocalHyperedgeCount(SInt hyperedge_size, SInt m_k);

    void GenerateRandomGeometricSizeCounts(
        SInt lower, SInt upper, double alpha, std::unordered_map<SInt, SInt>& size_counts);

    void GenerateDeterministicGeometricSizeCounts(
        SInt lower, SInt upper, double decay, std::unordered_map<SInt, SInt>& size_counts);

    SInt RankSplitSeed(SInt hyperedge_size, SInt rank_begin, SInt rank_end, SInt level) const;

    void GenerateBoltzmannPinBudgetSizeCounts(
        SInt lower_bound, SInt upper_bound, SInt pin_budget, std::unordered_map<SInt, SInt>& size_counts);

    SInt ExactLocalHyperedgeCountRecursive(
        SInt hyperedge_size, SInt rank_begin, SInt rank_end, SInt target_rank, CountInt population, SInt draws,
        SInt level);

    const PGeneratorConfig& config_;
    PEID                    rank_;
    PEID                    size_;

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