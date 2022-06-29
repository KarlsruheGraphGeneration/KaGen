/* Copyright (C) 2009-2010 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#pragma once

#include "hash.hpp"
#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/generator.h"

// forward declaration
struct mrg_state;

namespace kagen {
class KroneckerFactory : public GeneratorFactory {
public:
    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const override;
};

class Kronecker : public Generator {
public:
    Kronecker(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void GenerateImpl() final;

private:
    // Config
    const PGeneratorConfig& config_;
    PEID                    size_, rank_;

    // Constants and variables
    int     log_n_;
    SInt    from_, to_;
    SInt    num_edges_;
    int64_t scramble1_, scramble2_;
    SInt    edges_per_pe_;

    int Bernoulli(mrg_state* st, int level, int nlevels);

    /* Reverse bits in a number; this should be optimized for performance
     * (including using bit- or byte-reverse intrinsics if your platform has them).
     * */
    uint64_t bitreverse(uint64_t x);

    /* Apply a permutation to scramble vertex numbers; a randomly generated
     * permutation is not used because applying it at scale is too expensive. */
    int64_t Scramble(int64_t v0);

    /* Make a single graph edge using a pre-set MRG state. */
    void GenerateEdge(int64_t n, int level, mrg_state* st);
};
} // namespace kagen
