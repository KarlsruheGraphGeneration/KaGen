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

static Graph GenerateGNMWithDistribution(const std::string& generator, const std::string& distribution, SInt n, SInt m, int seed) {
    kagen::KaGen gen(MPI_COMM_WORLD);
    gen.UseEdgeListRepresentation();
    std::string options =
        generator + ";n=" + std::to_string(n) + ";m=" + std::to_string(m) + ";seed=" + std::to_string(seed)
        + ";distribution=" + distribution;
    return gen.GenerateFromOptionString(options);
}

// ---- Parameterized tests over generator, distribution, and edge density -----

using DistributionParam = std::tuple<std::string, std::string, double>;

struct GNMDistributionFixture : public ::testing::TestWithParam<DistributionParam> {};

INSTANTIATE_TEST_SUITE_P(
    GNMDistributionTests, GNMDistributionFixture,
    ::testing::Combine(
        ::testing::Values("gnm_directed", "gnm_undirected"),
        ::testing::Values("balance-vertices", "balance-edges", "balance-edges-strict"), ::testing::Values(0.5, 4.0)),
    [](const ::testing::TestParamInfo<DistributionParam>& info) {
        std::string generator = std::get<0>(info.param);
        std::string name      = std::get<1>(info.param);
        std::replace(name.begin(), name.end(), '-', '_');
        std::string factor = std::to_string(std::get<2>(info.param));
        std::replace(factor.begin(), factor.end(), '.', 'p');
        factor.erase(factor.find_last_not_of('0') + 1);
        return generator + "_" + name + "_m" + factor + "n";
    });

TEST_P(GNMDistributionFixture, OwnershipInvariant) {
    auto [generator, distribution, factor] = GetParam();

    if (distribution == "balance-edges-strict") {
        GTEST_SKIP() << "true edge-balancing intentionally splits a vertex's edges across PEs, so a single PE's "
                        "vertex range does not bound the tails of its edges";
    }

    const SInt n    = 1024;
    const SInt m    = static_cast<SInt>(factor * n);
    const int  seed = 42;

    Graph graph = GenerateGNMWithDistribution(generator, distribution, n, m, seed);

    for (const auto& edge: graph.edges) {
        auto u = edge.first;
        EXPECT_GE(u, graph.vertex_range.first);
        EXPECT_LT(u, graph.vertex_range.second);
    }
}

TEST_P(GNMDistributionFixture, NoDuplicates) {
    auto [generator, distribution, factor] = GetParam();
    const SInt n                           = 1024;
    const SInt m                           = static_cast<SInt>(factor * n);
    const int  seed                        = 42;

    Graph graph = GenerateGNMWithDistribution(generator, distribution, n, m, seed);

    Edgelist edges = graph.edges;
    std::sort(edges.begin(), edges.end());
    auto dup = std::adjacent_find(edges.begin(), edges.end());
    EXPECT_EQ(dup, edges.end());
}

// ---- Cross-distribution tests -----------------------------------------------

using CrossDistributionParam = std::tuple<std::string, double>;

struct GNMCrossDistributionFixture : public ::testing::TestWithParam<CrossDistributionParam> {};

INSTANTIATE_TEST_SUITE_P(
    GNMCrossDistributionTests, GNMCrossDistributionFixture,
    ::testing::Combine(::testing::Values("gnm_directed", "gnm_undirected"), ::testing::Values(0.5, 4.0)),
    [](const ::testing::TestParamInfo<CrossDistributionParam>& info) {
        std::string factor = std::to_string(std::get<1>(info.param));
        std::replace(factor.begin(), factor.end(), '.', 'p');
        factor.erase(factor.find_last_not_of('0') + 1);
        return std::get<0>(info.param) + "_m" + factor + "n";
    });

TEST_P(GNMCrossDistributionFixture, EdgeSetIdenticalAcrossDistributions) {
    auto [generator, factor] = GetParam();
    const SInt n              = 1024;
    const SInt m              = static_cast<SInt>(factor * n);
    const int  seed            = 42;

    Graph graph_bv  = GenerateGNMWithDistribution(generator, "balance-vertices", n, m, seed);
    Graph graph_be  = GenerateGNMWithDistribution(generator, "balance-edges", n, m, seed);
    Graph graph_bes = GenerateGNMWithDistribution(generator, "balance-edges-strict", n, m, seed);

    Edgelist edges_bv  = GatherSortedDeduplicatedEdges(graph_bv.edges);
    Edgelist edges_be  = GatherSortedDeduplicatedEdges(graph_be.edges);
    Edgelist edges_bes = GatherSortedDeduplicatedEdges(graph_bes.edges);

    EXPECT_EQ(edges_bv, edges_be);
    EXPECT_EQ(edges_bv, edges_bes);
}
