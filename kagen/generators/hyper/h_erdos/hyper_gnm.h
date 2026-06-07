#pragma once

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/kagen.h"
#include "kagen/tools/rng_wrapper.h"

#include <boost/multiprecision/cpp_int.hpp>

#include <vector>

namespace kagen {
using CountInt = boost::multiprecision::cpp_int;
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
    const PGeneratorConfig& config_;
    PEID                    rank_;
    PEID                    size_;

    RNGWrapper<BigInt> rng_;

    BigInt LocalBegin(BigInt universe_size) const;
    BigInt LocalEnd(BigInt universe_size) const;

    SInt HyperedgesOfSize(SInt k) const;

    CountInt Binomial(SInt n, SInt k) const;

    std::vector<SInt> UnrankCombination(BigInt index, SInt n, SInt k) const;

    void GenerateHyperedgesOfSize(SInt k, SInt m_k);

    void QueryCombinationRange(
        SInt k, SInt m, BigInt global_begin, BigInt global_end, BigInt query_begin, BigInt query_end, SInt level);

    void GenerateRange(SInt k, SInt m, BigInt global_begin, BigInt global_end, SInt level);

    void PushHyperedge(const std::vector<SInt>& pins);

    BigInt SaturationLimit(SInt m) const;

    CountInt CountFirstVertexRange(SInt lo, SInt hi, SInt n, SInt hyperedge_size) const;

    void QueryFirstVertexRange(SInt hyperedge_size, SInt m, SInt lo, SInt hi, SInt query_lo, SInt query_hi, SInt level);

    BigInt CheckedCastCount(const CountInt& value) const;

    SInt FirstVertexBegin(SInt hyperedge_size) const;

    SInt FirstVertexEnd(SInt hyperedge_size) const;

    void QueryPrefix(std::vector<SInt>& prefix, SInt remaining_k, SInt lo, SInt hi, SInt m, SInt level);

    std::vector<SInt> UnrankFirstVertexRange(BigInt index, SInt lo, SInt hi, SInt n, SInt k) const;

    BigInt BinomialNative(SInt n, SInt k) const;

    std::string ToString(__int128 value);
};

using HyperGNMSmall = HyperGNM<SInt>;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
using HyperGNMBig = HyperGNM<__int128>;
#pragma GCC diagnostic pop

} // namespace kagen