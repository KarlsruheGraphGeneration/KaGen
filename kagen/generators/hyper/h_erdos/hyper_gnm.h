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

template <typename BigInt>
class HyperGNM : public virtual Generator, private CSROnlyGenerator {
public:
    HyperGNM(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void GenerateCSR() final;
    void FinalizeCSR(MPI_Comm comm) final;

private:
    void GenerateHyperedgesOfSize(SInt hyperedge_size, SInt m_k);

    void PushHyperedge(const std::vector<SInt>& pins);

    SInt SampleHyperedgeSize(SInt lower, SInt upper, double alpha, RNGWrapper<>& rng, SInt& seed);

    SInt ApproximateLocalHyperedgeCount(SInt hyperedge_size, SInt m_k);

    SInt ApproximateLocalHyperedgeCountRecursive(
        SInt hyperedge_size, SInt rank_begin, SInt rank_end, SInt target_rank, long double population_mass, SInt draws,
        SInt level, LogBinomCache& log_binom_cache);

    SInt ExactLocalHyperedgeCount(SInt hyperedge_size, SInt m_k);

    void GenerateSizeCounts(SInt lower, SInt upper, double alpha, std::unordered_map<SInt, SInt>& size_counts);

    void GenerateDeterministicDecaySizeCounts(
        SInt lower, SInt upper, double decay, std::unordered_map<SInt, SInt>& size_counts);

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