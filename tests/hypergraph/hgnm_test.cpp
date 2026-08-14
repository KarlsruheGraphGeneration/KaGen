#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/hyper/h_erdos/hyper_gnm.h"

#include <gtest/gtest.h>
#include <mpi.h>

#include "hypergraph/utils.h"
#include <algorithm>
#include <map>
#include <set>
#include <vector>

using namespace kagen;

namespace {

Graph GenerateHGNM(PGeneratorConfig config) {
    HyperGNMFactory factory;

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

    config.n = 128;
    config.m = 64;

    //
    // Keep the number of virtual partitions fixed.
    //
    config.k = 8;

    config.seed = 1;

    config.size_dist_lower_bound = 2;
    config.size_dist_upper_bound = 8;

    config.allow_duplicates = false;
    config.approx           = false;

    return config;
}

std::map<SInt, SInt> CountHyperedgesBySize(const kagen::testing::ExpandedHypergraph& hypergraph) {
    std::map<SInt, SInt> counts;

    for (const auto& edge: hypergraph) {
        ++counts[static_cast<SInt>(edge.size())];
    }

    return counts;
}

void ExpectNoDuplicateHyperedges(kagen::testing::ExpandedHypergraph hypergraph) {
    //
    // ExpandLocalHypergraph() already sorts the pins inside
    // each hyperedge. Sorting the outer vector therefore gives
    // us a canonical representation.
    //
    std::sort(hypergraph.begin(), hypergraph.end());

    for (std::size_t i = 1; i < hypergraph.size(); ++i) {
        EXPECT_NE(hypergraph[i - 1], hypergraph[i]) << "duplicate hyperedge at canonical index " << i;
    }
}

kagen::testing::ExpandedHypergraph Canonicalize(kagen::testing::ExpandedHypergraph hypergraph) {
    for (auto& edge: hypergraph) {
        std::sort(edge.begin(), edge.end());
    }

    std::sort(hypergraph.begin(), hypergraph.end());

    return hypergraph;
}

} // namespace

TEST(HGNMTest, GeneratesExactlyMHyperedges) {
    auto config = BaseConfig();

    const Graph local_graph = GenerateHGNM(config);

    const auto global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    EXPECT_EQ(global_hypergraph.size(), static_cast<std::size_t>(config.m));
}

TEST(HGNMTest, FixedSizeProducesOnlyRequestedSize) {
    auto config = BaseConfig();

    config.size_dist_lower_bound = 5;
    config.size_dist_upper_bound = 5;

    const Graph local_graph = GenerateHGNM(config);

    const auto global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    ASSERT_EQ(global_hypergraph.size(), static_cast<std::size_t>(config.m));

    for (const auto& edge: global_hypergraph) {
        EXPECT_EQ(edge.size(), std::size_t{5});
    }
}

TEST(HGNMTest, HyperedgeSizesStayInsideConfiguredBounds) {
    auto config = BaseConfig();

    config.size_dist_lower_bound = 3;
    config.size_dist_upper_bound = 9;

    const Graph local_graph = GenerateHGNM(config);

    const auto global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    for (const auto& edge: global_hypergraph) {
        EXPECT_GE(edge.size(), std::size_t{3});

        EXPECT_LE(edge.size(), std::size_t{9});
    }
}

TEST(HGNMTest, HyperedgesContainNoDuplicatePins) {
    auto config = BaseConfig();

    const Graph local_graph = GenerateHGNM(config);

    const auto global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    for (std::size_t e = 0; e < global_hypergraph.size(); ++e) {
        const auto& edge = global_hypergraph[e];

        for (std::size_t i = 1; i < edge.size(); ++i) {
            EXPECT_NE(edge[i - 1], edge[i]) << "duplicate pin " << edge[i] << " in hyperedge " << e;
        }
    }
}

TEST(HGNMTest, HyperedgesAreUniqueByDefault) {
    auto config = BaseConfig();

    config.allow_duplicates = false;

    const Graph local_graph = GenerateHGNM(config);

    auto global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    ASSERT_EQ(global_hypergraph.size(), static_cast<std::size_t>(config.m));

    ExpectNoDuplicateHyperedges(std::move(global_hypergraph));
}

TEST(HGNMTest, PinBudgetIsSatisfiedExactly) {
    auto config = BaseConfig();

    config.m = 100;

    config.size_dist_lower_bound = 2;
    config.size_dist_upper_bound = 8;

    //
    // Feasible because:
    //
    //     100 * 2 <= 500 <= 100 * 8
    //
    config.size_dist_pin_budget = 500;

    const Graph local_graph = GenerateHGNM(config);

    const auto global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    ASSERT_EQ(global_hypergraph.size(), std::size_t{100});

    SInt total_pins = 0;

    for (const auto& edge: global_hypergraph) {
        total_pins += static_cast<SInt>(edge.size());

        EXPECT_GE(edge.size(), std::size_t{2});

        EXPECT_LE(edge.size(), std::size_t{8});
    }

    EXPECT_EQ(total_pins, config.size_dist_pin_budget);
}

TEST(HGNMTest, DeterministicSizeDistributionHasExpectedCounts) {
    auto config = BaseConfig();

    config.m = 100;

    config.size_dist_lower_bound = 2;
    config.size_dist_upper_bound = 4;

    config.size_dist_deterministic = true;
    config.size_decay              = 0.5;

    const Graph local_graph = GenerateHGNM(config);

    const auto global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    ASSERT_EQ(global_hypergraph.size(), std::size_t{100});

    const auto counts = CountHyperedgesBySize(global_hypergraph);

    //
    // HGNM's deterministic distribution performs:
    //
    // size 2: floor(100 * 0.5)       = 50
    // size 3: floor(100 * 0.5 * 0.5) = 25
    // size 4: remaining              = 25
    //
    EXPECT_EQ(counts.at(2), 50);

    EXPECT_EQ(counts.at(3), 25);

    EXPECT_EQ(counts.at(4), 25);
}

TEST(HGNMTest, SameSeedProducesSameHypergraph) {
    auto config = BaseConfig();

    config.seed = 42;

    const Graph first_local = GenerateHGNM(config);

    const auto first_global = kagen::testing::GatherHypergraph(first_local);

    const Graph second_local = GenerateHGNM(config);

    const auto second_global = kagen::testing::GatherHypergraph(second_local);

    if (!IsRoot()) {
        return;
    }

    EXPECT_EQ(Canonicalize(first_global), Canonicalize(second_global));
}

TEST(HGNMTest, ApproximateModePreservesBasicModelConstraints) {
    auto config = BaseConfig();

    config.approx = true;

    config.n = 256;
    config.m = 64;

    config.size_dist_lower_bound = 3;
    config.size_dist_upper_bound = 8;

    const Graph local_graph = GenerateHGNM(config);

    auto global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    EXPECT_EQ(global_hypergraph.size(), static_cast<std::size_t>(config.m));

    for (const auto& edge: global_hypergraph) {
        EXPECT_GE(edge.size(), std::size_t{3});

        EXPECT_LE(edge.size(), std::size_t{8});
    }

    ExpectNoDuplicateHyperedges(std::move(global_hypergraph));
}