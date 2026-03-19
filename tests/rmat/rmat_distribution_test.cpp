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

static Graph GenerateRMATWithDistribution(const std::string& distribution, SInt n, SInt m, int seed) {
    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    std::string options = "rmat;n=" + std::to_string(n) + ";m=" + std::to_string(m) + ";seed=" + std::to_string(seed)
                          + ";distribution=" + distribution;
    return generator.GenerateFromOptionString(options);
}

// ---- Parameterized tests over distribution and edge density -----------------

using DistributionParam = std::tuple<std::string, double>;

struct RMATDistributionFixture : public ::testing::TestWithParam<DistributionParam> {};

INSTANTIATE_TEST_SUITE_P(
    RMATDistributionTests, RMATDistributionFixture,
    ::testing::Combine(::testing::Values("balance-vertices", "balance-edges"), ::testing::Values(0.5, 4.0, 16.0)),
    [](const ::testing::TestParamInfo<DistributionParam>& info) {
        std::string name = std::get<0>(info.param);
        std::replace(name.begin(), name.end(), '-', '_');
        std::string factor = std::to_string(std::get<1>(info.param));
        std::replace(factor.begin(), factor.end(), '.', 'p');
        factor.erase(factor.find_last_not_of('0') + 1);
        return name + "_m" + factor + "n";
    });

TEST_P(RMATDistributionFixture, OwnershipInvariant) {
    auto [distribution, factor] = GetParam();
    const SInt n                = 1024;
    const SInt m                = static_cast<SInt>(factor * n);
    const int  seed             = 42;

    Graph graph = GenerateRMATWithDistribution(distribution, n, m, seed);

    for (const auto& edge: graph.edges) {
        auto u = edge.first;
        EXPECT_GE(u, graph.vertex_range.first);
        EXPECT_LT(u, graph.vertex_range.second);
    }
}

TEST_P(RMATDistributionFixture, NoDuplicates) {
    auto [distribution, factor] = GetParam();
    const SInt n                = 1024;
    const SInt m                = static_cast<SInt>(factor * n);
    const int  seed             = 42;

    Graph graph = GenerateRMATWithDistribution(distribution, n, m, seed);

    Edgelist edges = graph.edges;
    std::sort(edges.begin(), edges.end());
    auto dup = std::adjacent_find(edges.begin(), edges.end());
    EXPECT_EQ(dup, edges.end());
}

// ---- Cross-distribution tests -----------------------------------------------

struct RMATCrossDistributionFixture : public ::testing::TestWithParam<double> {};

INSTANTIATE_TEST_SUITE_P(
    RMATCrossDistributionTests, RMATCrossDistributionFixture, ::testing::Values(0.5, 4.0, 16.0),
    [](const ::testing::TestParamInfo<double>& info) {
        std::string factor = std::to_string(info.param);
        std::replace(factor.begin(), factor.end(), '.', 'p');
        factor.erase(factor.find_last_not_of('0') + 1);
        return "m" + factor + "n";
    });

TEST_P(RMATCrossDistributionFixture, EdgeSetIdenticalAcrossDistributions) {
    const SInt   n      = 1024;
    const double factor = GetParam();
    const SInt   m      = static_cast<SInt>(factor * n);
    const int    seed   = 42;

    Graph graph_bv = GenerateRMATWithDistribution("balance-vertices", n, m, seed);
    Graph graph_be = GenerateRMATWithDistribution("balance-edges", n, m, seed);

    Edgelist edges_bv = GatherSortedDeduplicatedEdges(graph_bv.edges);
    Edgelist edges_be = GatherSortedDeduplicatedEdges(graph_be.edges);

    EXPECT_EQ(edges_bv, edges_be);
}
