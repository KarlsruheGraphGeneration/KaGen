#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/generator.h"
#include "kagen/generators/geometric/rgg/rgg_2d.h"
#include "kagen/kagen.h"

#include <gtest/gtest.h>
#include <mpi.h>

#include "../gather.h"
#include <algorithm>
#include <string>
#include <vector>

using namespace kagen;

static Graph MakeEdgeListGraph(const Edgelist& edges) {
    Graph g;
    g.representation = GraphRepresentation::EDGE_LIST;
    g.edges          = edges;
    return g;
}

static Edgelist GatherSortedDeduplicatedEdges(const Edgelist& local_edges) {
    Edgelist edges = kagen::testing::GatherEdgeLists(MakeEdgeListGraph(local_edges)).edges;
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
}

static Graph
GenerateRGGWithDistribution(const std::string& generator, const std::string& distribution, SInt n, double r, int seed) {
    kagen::KaGen gen(MPI_COMM_WORLD);
    gen.UseEdgeListRepresentation();
    std::string options = generator + ";n=" + std::to_string(n) + ";r=" + std::to_string(r)
                          + ";seed=" + std::to_string(seed) + ";redistribution=" + distribution;
    return gen.GenerateFromOptionString(options);
}

// ---- Parameterized tests over generator, distribution, and radius -----------

using DistributionParam = std::tuple<std::string, std::string, double>;

struct RGGDistributionFixture : public ::testing::TestWithParam<DistributionParam> {};

INSTANTIATE_TEST_SUITE_P(
    RGGDistributionTests, RGGDistributionFixture,
    ::testing::Combine(
        ::testing::Values("rgg2d", "rgg3d"),
        ::testing::Values("balance-vertices", "balance-edges", "balance-edges-strict"),
        ::testing::Values(0.05, 0.1)),
    [](const ::testing::TestParamInfo<DistributionParam>& info) {
        std::string generator = std::get<0>(info.param);
        std::string name      = std::get<1>(info.param);
        std::replace(name.begin(), name.end(), '-', '_');
        std::string radius = std::to_string(std::get<2>(info.param));
        std::replace(radius.begin(), radius.end(), '.', 'p');
        radius.erase(radius.find_last_not_of('0') + 1);
        return generator + "_" + name + "_r" + radius;
    });

TEST_P(RGGDistributionFixture, OwnershipInvariant) {
    auto [generator, distribution, r] = GetParam();

    if (distribution == "balance-edges-strict") {
        GTEST_SKIP() << "true edge-balancing intentionally splits a vertex's edges across PEs, so a single PE's "
                        "vertex range does not bound the tails of its edges";
    }

    const SInt n    = 2000;
    const int  seed = 42;

    Graph graph = GenerateRGGWithDistribution(generator, distribution, n, r, seed);

    for (const auto& edge: graph.edges) {
        auto u = edge.first;
        EXPECT_GE(u, graph.vertex_range.first);
        EXPECT_LT(u, graph.vertex_range.second);
    }
}

TEST_P(RGGDistributionFixture, NoDuplicates) {
    auto [generator, distribution, r] = GetParam();
    const SInt n                      = 2000;
    const int  seed                   = 42;

    Graph graph = GenerateRGGWithDistribution(generator, distribution, n, r, seed);

    Edgelist edges = graph.edges;
    std::sort(edges.begin(), edges.end());
    auto dup = std::adjacent_find(edges.begin(), edges.end());
    EXPECT_EQ(dup, edges.end());
}

// ---- Cross-distribution tests -----------------------------------------------

using CrossDistributionParam = std::tuple<std::string, double>;

struct RGGCrossDistributionFixture : public ::testing::TestWithParam<CrossDistributionParam> {};

INSTANTIATE_TEST_SUITE_P(
    RGGCrossDistributionTests, RGGCrossDistributionFixture,
    ::testing::Combine(::testing::Values("rgg2d", "rgg3d"), ::testing::Values(0.05, 0.1)),
    [](const ::testing::TestParamInfo<CrossDistributionParam>& info) {
        std::string radius = std::to_string(std::get<1>(info.param));
        std::replace(radius.begin(), radius.end(), '.', 'p');
        radius.erase(radius.find_last_not_of('0') + 1);
        return std::get<0>(info.param) + "_r" + radius;
    });

TEST_P(RGGCrossDistributionFixture, EdgeSetIdenticalAcrossDistributions) {
    auto [generator, r] = GetParam();
    const SInt n         = 2000;
    const int  seed      = 42;

    Graph graph_bv  = GenerateRGGWithDistribution(generator, "balance-vertices", n, r, seed);
    Graph graph_be  = GenerateRGGWithDistribution(generator, "balance-edges", n, r, seed);
    Graph graph_bes = GenerateRGGWithDistribution(generator, "balance-edges-strict", n, r, seed);

    Edgelist edges_bv  = GatherSortedDeduplicatedEdges(graph_bv.edges);
    Edgelist edges_be  = GatherSortedDeduplicatedEdges(graph_be.edges);
    Edgelist edges_bes = GatherSortedDeduplicatedEdges(graph_bes.edges);

    EXPECT_EQ(edges_bv, edges_be);
    EXPECT_EQ(edges_bv, edges_bes);
}

// ---- Negative test: coordinates + edge-balanced redistribution ---------------

// GenerateFromOptionString funnels through GenerateInMemory, which catches ConfigurationError and calls
// MPI_Abort -- not observable by gtest. Drive the low-level Factory/Generator API directly instead so the
// exception can actually be caught here.
TEST(RGGCoordinatesGuardTest, CoordinatesWithBalanceEdgesThrows) {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    PGeneratorConfig config;
    config.n           = 2000;
    config.r           = 0.05;
    config.seed        = 42;
    config.coordinates = true;

    RGG2DFactory factory;
    config.redistribution = GraphRedistribution::BALANCE_EDGES;
    config                = factory.NormalizeParameters(config, rank, size, false);

    auto generator = factory.Create(config, rank, size);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    EXPECT_THROW(generator->Finalize(MPI_COMM_WORLD), ConfigurationError);
}

TEST(RGGCoordinatesGuardTest, CoordinatesWithBalanceEdgesStrictThrows) {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    PGeneratorConfig config;
    config.n           = 2000;
    config.r           = 0.05;
    config.seed        = 42;
    config.coordinates = true;

    RGG2DFactory factory;
    config.redistribution = GraphRedistribution::BALANCE_EDGES_TRUE;
    config                = factory.NormalizeParameters(config, rank, size, false);

    auto generator = factory.Create(config, rank, size);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    EXPECT_THROW(generator->Finalize(MPI_COMM_WORLD), ConfigurationError);
}

TEST(RGGCoordinatesGuardTest, CoordinatesWithBalanceVerticesSucceeds) {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    PGeneratorConfig config;
    config.n           = 2000;
    config.r           = 0.05;
    config.seed        = 42;
    config.coordinates = true;

    RGG2DFactory factory;
    config.redistribution = GraphRedistribution::BALANCE_VERTICES;
    config                = factory.NormalizeParameters(config, rank, size, false);

    auto generator = factory.Create(config, rank, size);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    EXPECT_NO_THROW(generator->Finalize(MPI_COMM_WORLD));
}
