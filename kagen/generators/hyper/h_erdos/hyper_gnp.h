#pragma once

#include "kagen/generators/generator.h"
#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"
#include "kagen/hypergraph/debug_logger_erdos.h"
#include "kagen/kagen.h"
#include "kagen/tools/rng_wrapper.h"

#include <optional>
#include <vector>

namespace kagen {

enum ProbabilityType : std::uint8_t { EXPLICIT_PROBS, EXPLICIT_EXPECTED, BUDGET_MODE, GLOBAL_PROBABILITY };

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
    void GenerateHyperedgesOfSize(SInt hyperedge_size, double p);

    void GenerateHyperedgesOfSizeApprox(SInt hyperedge_size, double p);

    void GenerateHyperedgesOfSizeExact(SInt hyperedge_size, double p);

    void GenerateLocalHyperedges(
        SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, SInt local_m, LogBinomCache& log_binom_cache);

    void PushHyperedge(const std::vector<SInt>& pins);

    std::pair<double, double> CalculateProbability(SInt hyperedge_size, SInt lower_bound);

    void SetProbability();

    SInt SetUpperBound(SInt hard_upper_bound, SInt lower_bound);

    void ValidateExactMinOwnerPartition(SInt hyperedge_size);

    void GenerateExactLocalRange(SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, double p, SInt level);

    SInt FindExactPopulationSplit(SInt begin, SInt end, SInt n, SInt k, const CountInt& total_population);

    void GenerateApproxLocalRange(
        SInt hyperedge_size, SInt local_min_begin, SInt local_min_end, double p, SInt level,
        LogBinomCache& count_cache);

    SInt FindApproxMassSplit(SInt begin, SInt end, SInt n, SInt k, long double total_mass, LogBinomCache& cache);

    PGeneratorConfig config_;

    PEID rank_;
    PEID size_;

    RNGWrapper<>    rng_;
    Mersenne        mersenne_;
    ProbabilityType probs_type_ = ProbabilityType::GLOBAL_PROBABILITY;

    std::optional<ErdosHypergraphDebugLogger> debug_logger_;
};
using HyperGNPSmall = HyperGNP<SInt>;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
using HyperGNPBig = HyperGNP<__int128>;
#pragma GCC diagnostic pop

} // namespace kagen