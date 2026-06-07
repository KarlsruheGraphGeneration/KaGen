#pragma once

#include "kagen/generators/generator.h"
#include "kagen/tools/rng_wrapper.h"

#include <boost/multiprecision/cpp_int.hpp>

#include <memory>
#include <vector>

namespace kagen {

using CountInt = boost::multiprecision::cpp_int;

class HyperGNPFactory : public GeneratorFactory {
public:
    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const;
};

template <typename BigInt>
class HyperGNP final : public virtual Generator, private CSROnlyGenerator {
public:
    HyperGNP(const PGeneratorConfig& config, PEID rank, PEID size);

    void GenerateCSR() final;
    void FinalizeCSR(MPI_Comm comm) final;

private:
    void GenerateHyperedgesOfSize(SInt hyperedge_size, double p);

    void GenerateHyperedgesOfSizeGNP(SInt hyperedge_size, const CountInt& universe, double p);

    void QueryPrefix(std::vector<SInt>& prefix, SInt remaining_k, SInt lo, SInt hi, SInt m, SInt level);

    void QueryPrefixGNP(std::vector<SInt>& prefix, SInt remaining_k, SInt lo, SInt hi, double p, SInt level);

    void
    QueryFirstVertexRangeGNP(SInt hyperedge_size, SInt lo, SInt hi, SInt query_lo, SInt query_hi, double p, SInt level);

    CountInt GenerateGeometricSkip(SInt seed, CountInt counter, double p) const;

    CountInt Binomial(SInt n, SInt k) const;

    CountInt CountFirstVertexRange(SInt lo, SInt hi, SInt n, SInt hyperedge_size) const;

    SInt FirstVertexBegin(SInt hyperedge_size) const;
    SInt FirstVertexEnd(SInt hyperedge_size) const;

    BigInt CheckedCastCount(const CountInt& value) const;

    std::vector<SInt> UnrankFirstVertexRange(BigInt index, SInt lo, SInt hi, SInt n, SInt k) const;

    BigInt BinomialNative(SInt n, SInt k) const;

    std::vector<SInt> UnrankCombination(BigInt index, SInt n, SInt k) const;

    void PushHyperedge(const std::vector<SInt>& pins);

private:
    PGeneratorConfig config_;
    PEID             rank_;
    PEID             size_;
    RNGWrapper<>     rng_;
};

using HyperGNPBig = HyperGNP<__int128>;

} // namespace kagen