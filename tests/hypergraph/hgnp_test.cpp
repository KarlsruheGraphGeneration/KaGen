#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/hyper/h_erdos/hyper_gnp.h"

#include <gtest/gtest.h>
#include <mpi.h>

#include "hypergraph/utils.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <vector>

using namespace kagen;

namespace {

Graph GenerateHGNP(PGeneratorConfig config) {
    HyperGNPFactory factory;

    PEID rank;
    PEID size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    config = factory.NormalizeParameters(config, rank, size, false);

    auto generator = factory.Create(config, rank, size);

    generator->Generate(GraphRepresentation::CSR);

    generator->Finalize(MPI_COMM_WORLD);

    return generator->Take();
}

bool IsRoot() {
    PEID rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    return rank == 0;
}

PGeneratorConfig BaseConfig() {
    PGeneratorConfig config;

    config.n = 64;

    //
    // Virtual minimum-owner partitions.
    //
    config.k = 8;

    config.seed = 1;

    config.size_dist_lower_bound = 2;
    config.size_dist_upper_bound = 4;

    config.allow_duplicates = false;
    config.approx           = false;
    config.quiet            = true;

    return config;
}

kagen::testing::ExpandedHypergraph Canonicalize(kagen::testing::ExpandedHypergraph hypergraph) {
    for (auto& edge: hypergraph) {
        std::sort(edge.begin(), edge.end());
    }

    std::sort(hypergraph.begin(), hypergraph.end());

    return hypergraph;
}

std::map<SInt, SInt> CountHyperedgesBySize(const kagen::testing::ExpandedHypergraph& hypergraph) {
    std::map<SInt, SInt> counts;

    for (const auto& edge: hypergraph) {
        ++counts[static_cast<SInt>(edge.size())];
    }

    return counts;
}

SInt CountOfSize(const std::map<SInt, SInt>& counts, const SInt size) {
    const auto it = counts.find(size);

    if (it == counts.end()) {
        return 0;
    }

    return it->second;
}

SInt NumberOfPins(const kagen::testing::ExpandedHypergraph& hypergraph) {
    SInt pins = 0;

    for (const auto& edge: hypergraph) {
        pins += static_cast<SInt>(edge.size());
    }

    return pins;
}

SInt BinomialCoefficient(const SInt n, const SInt k) {
    if (k < 0 || k > n) {
        return 0;
    }

    SInt result = 1;

    const SInt effective_k = std::min(k, n - k);

    for (SInt i = 1; i <= effective_k; ++i) {
        result = result * (n - effective_k + i) / i;
    }

    return result;
}

void ExpectBinomialCountNear(const SInt actual, const SInt population, const double probability) {
    const double expected = static_cast<double>(population) * probability;

    const double variance = static_cast<double>(population) * probability * (1.0 - probability);

    const double sigma = std::sqrt(variance);

    //
    // Six sigma makes these tests extremely unlikely to fail
    // because of ordinary sampling noise.
    //
    // The small additive margin also handles very small
    // expected counts without making the test brittle.
    //
    const double tolerance = std::max(3.0, 6.0 * sigma);

    EXPECT_NEAR(static_cast<double>(actual), expected, tolerance)
        << "population=" << population << " p=" << probability << " expected=" << expected << " sigma=" << sigma;
}

void ExpectValidPins(const kagen::testing::ExpandedHypergraph& hypergraph, const SInt n) {
    for (std::size_t e = 0; e < hypergraph.size(); ++e) {
        const auto& edge = hypergraph[e];

        for (const SInt pin: edge) {
            EXPECT_GE(pin, 0) << "negative pin in hyperedge " << e;

            EXPECT_LT(pin, n) << "invalid pin in hyperedge " << e;
        }
    }
}

void ExpectNoDuplicatePins(const kagen::testing::ExpandedHypergraph& hypergraph) {
    for (std::size_t e = 0; e < hypergraph.size(); ++e) {
        const auto& edge = hypergraph[e];

        for (std::size_t i = 1; i < edge.size(); ++i) {
            EXPECT_LT(edge[i - 1], edge[i]) << "pins are not strictly increasing "
                                            << "in hyperedge " << e;
        }
    }
}

void ExpectNoDuplicateHyperedges(kagen::testing::ExpandedHypergraph hypergraph) {
    hypergraph = Canonicalize(std::move(hypergraph));

    for (std::size_t i = 1; i < hypergraph.size(); ++i) {
        EXPECT_NE(hypergraph[i - 1], hypergraph[i]) << "duplicate hyperedge at canonical index " << i;
    }
}

} // namespace

//
// Global probability mode
//

TEST(HGNPTest, ZeroProbabilityProducesNoHyperedges) {
    auto config = BaseConfig();

    config.p = 0.0;

    const Graph local_graph = GenerateHGNP(config);

    const auto global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    EXPECT_TRUE(global_hypergraph.empty());
}

TEST(HGNPTest, GlobalProbabilityProducesExpectedSparseCounts) {
    auto config = BaseConfig();

    config.n = 64;

    config.size_dist_lower_bound = 2;
    config.size_dist_upper_bound = 4;

    //
    // Sparse enough that duplicate rejection remains cheap.
    //
    config.p = 0.002;

    const Graph local_graph = GenerateHGNP(config);

    const auto global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    const auto counts = CountHyperedgesBySize(global_hypergraph);

    for (SInt size = 2; size <= 4; ++size) {
        ExpectBinomialCountNear(CountOfSize(counts, size), BinomialCoefficient(config.n, size), config.p);
    }

    ExpectValidPins(global_hypergraph, config.n);

    ExpectNoDuplicatePins(global_hypergraph);

    ExpectNoDuplicateHyperedges(global_hypergraph);
}

//
// Explicit probability mode
//

TEST(HGNPTest, ExplicitProbabilitiesApplyToCorrespondingSizes) {
    auto config = BaseConfig();

    config.n = 64;

    config.size_dist_lower_bound = 2;

    //
    // Index zero corresponds to size_dist_lower_bound.
    //
    // All probabilities remain sparse.
    //
    config.size_probabilities = {
        0.01L,
        0.002L,
        0.0002L,
    };

    const Graph local_graph = GenerateHGNP(config);

    const auto global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    const auto counts = CountHyperedgesBySize(global_hypergraph);

    for (std::size_t i = 0; i < config.size_probabilities.size(); ++i) {
        const SInt size = config.size_dist_lower_bound + static_cast<SInt>(i);

        const double probability = static_cast<double>(config.size_probabilities[i]);

        ExpectBinomialCountNear(CountOfSize(counts, size), BinomialCoefficient(config.n, size), probability);
    }
}

//
// Explicit expected-count mode
//

TEST(HGNPTest, ExplicitExpectedCountsProduceExpectedSparseCounts) {
    auto config = BaseConfig();

    config.n = 64;

    config.size_dist_lower_bound = 2;

    //
    // These are expectations, not exact realized counts.
    //
    config.size_expected_counts = {
        20.0L,
        40.0L,
        60.0L,
    };

    const Graph local_graph = GenerateHGNP(config);

    const auto global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    const auto counts = CountHyperedgesBySize(global_hypergraph);

    for (std::size_t i = 0; i < config.size_expected_counts.size(); ++i) {
        const SInt size = config.size_dist_lower_bound + static_cast<SInt>(i);

        const SInt population = BinomialCoefficient(config.n, size);

        const double expected_count = static_cast<double>(config.size_expected_counts[i]);

        const double probability = expected_count / static_cast<double>(population);

        ExpectBinomialCountNear(CountOfSize(counts, size), population, probability);
    }
}

//
// Edge-budget mode
//

TEST(HGNPTest, SparseEdgeBudgetProducesExpectedTotalCount) {
    auto config = BaseConfig();

    config.n = 64;

    config.size_dist_lower_bound = 2;
    config.size_dist_upper_bound = 6;

    //
    // Geometric expected-count distribution.
    //
    config.edge_budget = 200.0;
    config.size_decay  = 0.5;

    const Graph local_graph = GenerateHGNP(config);

    const auto global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    //
    // edge_budget is an expectation.
    //
    // For independent Bernoulli generation across the size
    // classes the total variance is approximately the sum of
    // the per-class expected counts in the sparse regime.
    //
    const double expected = config.edge_budget;

    const double tolerance = 6.0 * std::sqrt(expected);

    EXPECT_NEAR(static_cast<double>(global_hypergraph.size()), expected, tolerance);

    for (const auto& edge: global_hypergraph) {
        EXPECT_GE(edge.size(), std::size_t{2});

        EXPECT_LE(edge.size(), std::size_t{6});
    }
}

//
// Edge + pin budget mode
//

TEST(HGNPTest, EdgeAndPinBudgetMatchesRequestedExpectations) {
    auto config = BaseConfig();

    config.n = 128;

    config.size_dist_lower_bound = 2;
    config.size_dist_upper_bound = 8;

    //
    // Expected budgets. These remain very sparse compared with
    // the available combinatorial populations.
    //
    config.edge_budget          = 500.0;
    config.size_dist_pin_budget = 2000;

    const Graph local_graph = GenerateHGNP(config);

    const auto global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    const double actual_edges = static_cast<double>(global_hypergraph.size());

    const double actual_pins = static_cast<double>(NumberOfPins(global_hypergraph));

    //
    // These are intentionally broad. The purpose is to detect
    // a broken budget interpretation, not to make a statistical
    // goodness-of-fit test out of a single realization.
    //
    kagen::testing::ExpectNearRelative(actual_edges, config.edge_budget, 0.20);

    kagen::testing::ExpectNearRelative(actual_pins, static_cast<double>(config.size_dist_pin_budget), 0.20);

    for (const auto& edge: global_hypergraph) {
        EXPECT_GE(edge.size(), std::size_t{2});

        EXPECT_LE(edge.size(), std::size_t{8});
    }
}

//
// General model invariants
//

TEST(HGNPTest, HyperedgesContainNoDuplicatePins) {
    auto config = BaseConfig();

    config.n = 64;

    config.size_dist_lower_bound = 3;
    config.size_dist_upper_bound = 6;

    config.p = 0.0002;

    const Graph local_graph = GenerateHGNP(config);

    const auto global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    ASSERT_FALSE(global_hypergraph.empty());

    ExpectNoDuplicatePins(global_hypergraph);

    ExpectValidPins(global_hypergraph, config.n);
}

TEST(HGNPTest, HyperedgesAreUniqueWhenDuplicatesAreDisabled) {
    auto config = BaseConfig();

    config.n = 64;

    config.size_dist_lower_bound = 2;
    config.size_dist_upper_bound = 4;

    config.p = 0.001;

    config.allow_duplicates = false;

    const Graph local_graph = GenerateHGNP(config);

    auto global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    ASSERT_FALSE(global_hypergraph.empty());

    ExpectNoDuplicateHyperedges(std::move(global_hypergraph));
}

//
// Reproducibility
//

TEST(HGNPTest, SameSeedProducesSameHypergraph) {
    auto config = BaseConfig();

    config.n = 64;

    config.size_dist_lower_bound = 2;
    config.size_dist_upper_bound = 4;

    config.p = 0.001;

    config.seed = 42;

    const Graph first_local = GenerateHGNP(config);

    const auto first_global = kagen::testing::GatherHypergraph(first_local);

    const Graph second_local = GenerateHGNP(config);

    const auto second_global = kagen::testing::GatherHypergraph(second_local);

    if (!IsRoot()) {
        return;
    }

    EXPECT_EQ(Canonicalize(first_global), Canonicalize(second_global));
}

TEST(HGNPTest, DifferentSeedsUsuallyProduceDifferentHypergraphs) {
    auto first_config = BaseConfig();

    auto second_config = BaseConfig();

    first_config.n = second_config.n = 64;

    first_config.size_dist_lower_bound = second_config.size_dist_lower_bound = 2;

    first_config.size_dist_upper_bound = second_config.size_dist_upper_bound = 4;

    first_config.p = second_config.p = 0.001;

    first_config.seed  = 1;
    second_config.seed = 42;

    const Graph first_local = GenerateHGNP(first_config);

    const auto first_global = kagen::testing::GatherHypergraph(first_local);

    const Graph second_local = GenerateHGNP(second_config);

    const auto second_global = kagen::testing::GatherHypergraph(second_local);

    if (!IsRoot()) {
        return;
    }

    ASSERT_FALSE(first_global.empty());

    ASSERT_FALSE(second_global.empty());

    EXPECT_NE(Canonicalize(first_global), Canonicalize(second_global));
}

//
// Approximate generation
//

TEST(HGNPTest, ApproximateModeRespectsModelInvariants) {
    auto config = BaseConfig();

    config.n = 128;

    config.size_dist_lower_bound = 3;
    config.size_dist_upper_bound = 6;

    //
    // Keep the realization sparse. Approximate mode is tested
    // here for structural/model invariants, not exact equality
    // with the exact generator.
    //
    config.p = 0.00001;

    config.approx           = true;
    config.allow_duplicates = false;

    const Graph local_graph = GenerateHGNP(config);

    auto global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    ASSERT_FALSE(global_hypergraph.empty());

    for (const auto& edge: global_hypergraph) {
        EXPECT_GE(edge.size(), std::size_t{3});

        EXPECT_LE(edge.size(), std::size_t{6});
    }

    ExpectValidPins(global_hypergraph, config.n);

    ExpectNoDuplicatePins(global_hypergraph);

    ExpectNoDuplicateHyperedges(std::move(global_hypergraph));
}