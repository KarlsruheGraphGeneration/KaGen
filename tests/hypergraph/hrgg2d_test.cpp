#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/hyper/h_geometric/h_rgg.h"

#include <gtest/gtest.h>
#include <mpi.h>

#include "hypergraph/utils.h"
#include <algorithm>
#include <cstddef>
#include <unordered_set>
#include <utility>

using namespace kagen;

namespace {

Graph GenerateHRGG2D(PGeneratorConfig config) {
    HyperRGG2DFactory factory;

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

    config.n = 1024;
    config.m = 64;

    //
    // HRGG2D requires a square power-of-two chunk layout.
    //
    config.k = 4;

    config.seed  = 1;
    config.quiet = true;
    config.debug = false;

    //
    // Constant-radius default for tests that do not specifically
    // exercise random radii.
    //
    config.random_radius = false;
    config.r             = 0.05;

    config.partial_cell_mode = PartialCellMode::GenerateAndCheck;

    return config;
}

kagen::testing::ExpandedHypergraph Canonicalize(kagen::testing::ExpandedHypergraph hypergraph) {
    for (auto& edge: hypergraph) {
        std::sort(edge.begin(), edge.end());
    }

    std::sort(hypergraph.begin(), hypergraph.end());

    return hypergraph;
}

void ExpectValidHypergraph(const kagen::testing::ExpandedHypergraph& hypergraph, const SInt num_vertices) {
    for (std::size_t e = 0; e < hypergraph.size(); ++e) {
        const auto& edge = hypergraph[e];

        std::unordered_set<SInt> seen;

        for (const SInt pin: edge) {
            EXPECT_LT(pin, num_vertices) << "invalid pin in hyperedge " << e;

            EXPECT_TRUE(seen.insert(pin).second) << "duplicate pin " << pin << " in hyperedge " << e;
        }
    }
}

kagen::testing::ExpandedHypergraph GenerateAndGather(PGeneratorConfig config) {
    const Graph local_graph = GenerateHRGG2D(std::move(config));

    return kagen::testing::GatherHypergraph(local_graph);
}

} // namespace

TEST(HRGG2DTest, GeneratesValidHypergraph) {
    const auto config = BaseConfig();

    const auto hypergraph = GenerateAndGather(config);

    if (!IsRoot()) {
        return;
    }

    ASSERT_EQ(hypergraph.size(), static_cast<std::size_t>(config.m));

    ExpectValidHypergraph(hypergraph, config.n);
}

TEST(HRGG2DTest, GeneratesExactNumberOfHyperedges) {
    auto config = BaseConfig();

    config.m = 137;

    const auto hypergraph = GenerateAndGather(config);

    if (!IsRoot()) {
        return;
    }

    EXPECT_EQ(hypergraph.size(), static_cast<std::size_t>(config.m));
}

TEST(HRGG2DTest, SameSeedProducesSameHypergraph) {
    auto config = BaseConfig();

    config.seed = 42;

    const auto first = GenerateAndGather(config);

    const auto second = GenerateAndGather(config);

    if (!IsRoot()) {
        return;
    }

    EXPECT_EQ(Canonicalize(first), Canonicalize(second));
}

TEST(HRGG2DTest, DifferentSeedsProduceDifferentHypergraphs) {
    auto first_config = BaseConfig();

    auto second_config = BaseConfig();

    first_config.seed  = 1;
    second_config.seed = 42;

    const auto first = GenerateAndGather(first_config);

    const auto second = GenerateAndGather(second_config);

    if (!IsRoot()) {
        return;
    }

    EXPECT_NE(Canonicalize(first), Canonicalize(second));
}

TEST(HRGG2DTest, GenerateAndCheckModeWorks) {
    auto config = BaseConfig();

    config.partial_cell_mode = PartialCellMode::GenerateAndCheck;

    const auto hypergraph = GenerateAndGather(config);

    if (!IsRoot()) {
        return;
    }

    ASSERT_EQ(hypergraph.size(), static_cast<std::size_t>(config.m));

    ExpectValidHypergraph(hypergraph, config.n);
}

TEST(HRGG2DTest, CoverageRangeModeWorks) {
    auto config = BaseConfig();

    config.partial_cell_mode = PartialCellMode::EstimateByCoverageRange;

    const auto hypergraph = GenerateAndGather(config);

    if (!IsRoot()) {
        return;
    }

    ASSERT_EQ(hypergraph.size(), static_cast<std::size_t>(config.m));

    ExpectValidHypergraph(hypergraph, config.n);
}

TEST(HRGG2DTest, CoverageFloydModeWorks) {
    auto config = BaseConfig();

    config.partial_cell_mode = PartialCellMode::EstimateByCoverageFloyd;

    const auto hypergraph = GenerateAndGather(config);

    if (!IsRoot()) {
        return;
    }

    ASSERT_EQ(hypergraph.size(), static_cast<std::size_t>(config.m));

    ExpectValidHypergraph(hypergraph, config.n);
}

TEST(HRGG2DTest, ConstantRadiusWorks) {
    auto config = BaseConfig();

    config.random_radius = false;
    config.r             = 0.04;

    const auto hypergraph = GenerateAndGather(config);

    if (!IsRoot()) {
        return;
    }

    ASSERT_EQ(hypergraph.size(), static_cast<std::size_t>(config.m));

    ExpectValidHypergraph(hypergraph, config.n);
}

TEST(HRGG2DTest, RandomRadiusWorks) {
    auto config = BaseConfig();

    config.random_radius = true;

    //
    // Let normalization choose its radius used for the grid.
    //
    config.r = 0.0;

    config.size_dist_lower_bound = 2;
    config.size_dist_upper_bound = 32;

    config.hyperedge_radius_exponent = 2.5;

    const auto hypergraph = GenerateAndGather(config);

    if (!IsRoot()) {
        return;
    }

    ASSERT_EQ(hypergraph.size(), static_cast<std::size_t>(config.m));

    ExpectValidHypergraph(hypergraph, config.n);
}

TEST(HRGG2DTest, PinBudgetModeWorks) {
    auto config = BaseConfig();

    config.random_radius = true;
    config.r             = 0.0;

    config.size_dist_lower_bound = 2;
    config.size_dist_upper_bound = 32;

    //
    // m = 64, pin budget = 512:
    //
    //     expected mean hyperedge size = 8.
    //
    // This lies comfortably inside [2, 32].
    //
    config.size_dist_pin_budget = 512;

    const auto hypergraph = GenerateAndGather(config);

    if (!IsRoot()) {
        return;
    }

    ASSERT_EQ(hypergraph.size(), static_cast<std::size_t>(config.m));

    ExpectValidHypergraph(hypergraph, config.n);
}