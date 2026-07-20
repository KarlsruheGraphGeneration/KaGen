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

// NOTE: for *generators* the option-string key is `redistribution` (GraphRedistribution); `distribution` is only
// honored for file inputs (GraphDistribution) and is silently ignored here -- passing it leaves the graph at the
// default redistribution, so every "distribution" would generate the same graph.
static Graph GenerateRMATWithDistribution(const std::string& distribution, SInt n, SInt m, int seed) {
    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    std::string options = "rmat;n=" + std::to_string(n) + ";m=" + std::to_string(m) + ";seed=" + std::to_string(seed)
                          + ";redistribution=" + distribution;
    return generator.GenerateFromOptionString(options);
}

static Graph GenerateRMATCSRWithDistribution(const std::string& distribution, SInt n, SInt m, int seed) {
    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseCSRRepresentation();
    std::string options = "rmat;n=" + std::to_string(n) + ";m=" + std::to_string(m) + ";seed=" + std::to_string(seed)
                          + ";redistribution=" + distribution;
    return generator.GenerateFromOptionString(options);
}

// Reconstruct this PE's (from, to) edges from a CSR graph. Row i is vertex vertex_range.first + i -- so this
// only yields the right sources if vertex_range is the physically-present row space that xadj indexes (which,
// for a left-partial split PE, means vertex_range.first is the shared replica vertex, one below the gap-free
// range). Mirrors BuildEdgeListFromCSR.
static Edgelist EdgesFromCSR(const Graph& g) {
    Edgelist edges;
    for (std::size_t i = 0; i + 1 < g.xadj.size(); ++i) {
        const SInt from = g.vertex_range.first + static_cast<SInt>(i);
        for (SInt e = g.xadj[i]; e < g.xadj[i + 1]; ++e) {
            edges.emplace_back(from, g.adjncy[e]);
        }
    }
    return edges;
}

// ---- Parameterized tests over distribution and edge density -----------------

using DistributionParam = std::tuple<std::string, double>;

struct RMATDistributionFixture : public ::testing::TestWithParam<DistributionParam> {};

INSTANTIATE_TEST_SUITE_P(
    RMATDistributionTests, RMATDistributionFixture,
    ::testing::Combine(
        ::testing::Values("balance-vertices", "balance-edges", "balance-edges-strict"),
        ::testing::Values(0.5, 4.0, 16.0)),
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

    if (distribution == "balance-edges-strict") {
        GTEST_SKIP() << "true edge-balancing intentionally splits a vertex's edges across PEs, so a single PE's "
                        "vertex range does not bound the tails of its edges";
    }

    const SInt n    = 1024;
    const SInt m    = static_cast<SInt>(factor * n);
    const int  seed = 42;

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

// ---- Strict edge-balanced CSR generation ------------------------------------

// Regression test for the defect where RMAT + balance-edges-strict + CSR threw ("CSR representation requires
// single-PE vertex ownership ..."): EdgeListOnlyGenerator::FinalizeCSR now builds a split CSR by extending the
// row space down by one for a left-partial replica vertex, instead of underflowing/throwing.
TEST(RMATStrictEdgeBalanceCSR, ProducesConsumableSplitCSR) {
    const SInt n    = 1024;
    const SInt m    = 16 * n; // dense hub-heavy RMAT: splits are essentially guaranteed at 4 PEs
    const int  seed = 42;

    // Previously this threw; now it must succeed and yield a CSR graph.
    Graph csr;
    ASSERT_NO_THROW(csr = GenerateRMATCSRWithDistribution("balance-edges-strict", n, m, seed));
    ASSERT_EQ(csr.representation, GraphRepresentation::CSR);
    ASSERT_FALSE(csr.xadj.empty());

    // Reconstructing the edges from the CSR (row i -> vertex_range.first + i) must yield exactly the same global
    // edge set as generating the same graph directly as an edge list. Holds whether or not a vertex actually
    // splits at this PE count, and fails if the row space / vertex_range is wrong (e.g. the pre-fix gap-free
    // range, which mislabels the leading replica row's source), so it pins the fix.
    Graph    el       = GenerateRMATWithDistribution("balance-edges-strict", n, m, seed);
    Edgelist from_csr = GatherSortedDeduplicatedEdges(EdgesFromCSR(csr));
    Edgelist from_el  = GatherSortedDeduplicatedEdges(el.edges);
    EXPECT_EQ(from_csr, from_el);

    // Whether this config actually splits a vertex depends on the PE count (guaranteed at higher counts; may not
    // at 1-4 PEs, where the top RMAT hub still fits within one PE's edge budget) -- deterministic split-CSR
    // coverage lives in StrictEdgeBalanceTest. When a split does occur, its metadata must be well-formed: every
    // partial vertex sits inside the row-space range, its edge block fits, and a left-partial replica's vertex is
    // exactly vertex_range.first (row 0 of the extended row space).
    if (csr.has_split_vertices) {
        for (const auto& pv: {csr.left_partial_vertex, csr.right_partial_vertex}) {
            if (!pv) {
                continue;
            }
            EXPECT_GE(pv->vertex, csr.vertex_range.first);
            EXPECT_LT(pv->vertex, csr.vertex_range.second);
            EXPECT_LE(pv->local_offset + pv->local_count, static_cast<SInt>(csr.adjncy.size()));
        }
        if (csr.left_partial_vertex) {
            EXPECT_EQ(csr.left_partial_vertex->vertex, csr.vertex_range.first);
        }
    }
}

TEST_P(RMATCrossDistributionFixture, EdgeSetIdenticalAcrossDistributions) {
    const SInt   n      = 1024;
    const double factor = GetParam();
    const SInt   m      = static_cast<SInt>(factor * n);
    const int    seed   = 42;

    Graph graph_bv  = GenerateRMATWithDistribution("balance-vertices", n, m, seed);
    Graph graph_be  = GenerateRMATWithDistribution("balance-edges", n, m, seed);
    Graph graph_bes = GenerateRMATWithDistribution("balance-edges-strict", n, m, seed);

    Edgelist edges_bv  = GatherSortedDeduplicatedEdges(graph_bv.edges);
    Edgelist edges_be  = GatherSortedDeduplicatedEdges(graph_be.edges);
    Edgelist edges_bes = GatherSortedDeduplicatedEdges(graph_bes.edges);

    EXPECT_EQ(edges_bv, edges_be);
    // RMAT symmetrizes every directed sample independently per chunk, so it can legitimately draw the same
    // undirected edge twice on two different PEs; balance-edges-strict routes a split vertex's edges by
    // position rather than by co-location, so this specifically exercises that its dedup still produces the
    // exact same edge set as the other two distributions (see RedistributeEdgesTrueBalance's escalation path).
    EXPECT_EQ(edges_bv, edges_bes);
}
