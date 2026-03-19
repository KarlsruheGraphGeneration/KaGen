#include "kagen/kagen.h"
#include "kagen/tools/postprocessor.h"

#include <gtest/gtest.h>
#include <mpi.h>

#include "../gather.h"
#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

using namespace kagen;

using GeneratorFunc = std::function<Graph(KaGen&, SInt, SInt)>;

static Graph MakeEdgeListGraph(const Edgelist& edges) {
    Graph g;
    g.representation = GraphRepresentation::EDGE_LIST;
    g.edges          = edges;
    return g;
}

static Edgelist GatherAllEdges(const Edgelist& local_edges) {
    return kagen::testing::GatherEdgeLists(MakeEdgeListGraph(local_edges)).edges;
}

static std::vector<SInt> GatherEdgeCounts(SInt local_count) {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<SInt> counts(size);
    MPI_Allgather(&local_count, 1, KAGEN_MPI_SINT, counts.data(), 1, KAGEN_MPI_SINT, MPI_COMM_WORLD);
    return counts;
}

// Build a star graph (vertex 0 connected to all others, both directions) with all edges on PE 0.
static Edgelist BuildStarOnPE0(SInt n) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Edgelist edges;
    if (rank == 0) {
        for (SInt v = 1; v < n; ++v) {
            edges.emplace_back(0, v);
            edges.emplace_back(v, 0);
        }
    }
    return edges;
}

// ---- Parameterized tests using KaGen generators ----------------------------

struct RedistributeEdgesBalancedFixture : public ::testing::TestWithParam<std::tuple<std::string, GeneratorFunc>> {};

INSTANTIATE_TEST_SUITE_P(
    RedistributeEdgesBalancedTests, RedistributeEdgesBalancedFixture,
    ::testing::Values(
        std::make_tuple("GNM", GeneratorFunc([](KaGen& gen, SInt n, SInt m) {
                            return gen.GenerateUndirectedGNM(n, m);
                        })),
        std::make_tuple("RMAT", GeneratorFunc([](KaGen& gen, SInt n, SInt m) {
                            return gen.GenerateRMAT(n, m, 0.56, 0.19, 0.19);
                        })),
        std::make_tuple("RGG2D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRGG2D_NM(n, m); })),
        std::make_tuple("RHG", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRHG_NM(2.8, n, m); })),
        std::make_tuple("Grid2D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) {
                            return gen.GenerateGrid2D_N(n, 1.0);
                        }))),
    [](const ::testing::TestParamInfo<RedistributeEdgesBalancedFixture::ParamType>& info) {
        return std::get<0>(info.param);
    });

// Apply the same round-robin remapping to the original edges, then compare
// exact edge sets with the balanced output.
TEST_P(RedistributeEdgesBalancedFixture, PreservesEdgeSet) {
    auto [name, generate] = GetParam();
    const SInt n          = 1000;
    const SInt m          = 4 * n;

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    Graph graph = generate(generator, n, m);

    Edgelist input        = graph.edges;
    SInt     num_vertices = graph.vertex_range.second;
    MPI_Allreduce(MPI_IN_PLACE, &num_vertices, 1, KAGEN_MPI_SINT, MPI_MAX, MPI_COMM_WORLD);

    // Remap the original edges with the same round-robin mapping that
    // RedistributeEdgesBalanced applies internally.
    Edgelist expected = [&num_vertices](auto edges) {
        RoundRobinRemapping(edges, num_vertices, MPI_COMM_WORLD);
        edges = GatherAllEdges(edges);
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
        return edges;
    }(input);

    Edgelist redistributed_edges;
    RedistributeEdgesBalanced(input, redistributed_edges, num_vertices, MPI_COMM_WORLD);

    Edgelist result = GatherAllEdges(redistributed_edges);
    std::sort(result.begin(), result.end());

    EXPECT_EQ(expected, result);
}

TEST_P(RedistributeEdgesBalancedFixture, OwnershipInvariant) {
    auto [name, generate] = GetParam();
    const SInt n          = 1000;
    const SInt m          = 4 * n;

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    Graph graph = generate(generator, n, m);

    Edgelist input        = graph.edges;
    SInt     num_vertices = graph.vertex_range.second;
    MPI_Allreduce(MPI_IN_PLACE, &num_vertices, 1, KAGEN_MPI_SINT, MPI_MAX, MPI_COMM_WORLD);

    Edgelist    redistributed_edges;
    VertexRange vr = RedistributeEdgesBalanced(input, redistributed_edges, num_vertices, MPI_COMM_WORLD);

    for (const auto& [u, v]: redistributed_edges) {
        EXPECT_GE(u, vr.first);
        EXPECT_LT(u, vr.second);
    }
}

TEST_P(RedistributeEdgesBalancedFixture, NoDuplicatesInOutput) {
    auto [name, generate] = GetParam();
    const SInt n          = 1000;
    const SInt m          = 4 * n;

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    Graph graph = generate(generator, n, m);

    Edgelist input        = graph.edges;
    SInt     num_vertices = graph.vertex_range.second;
    MPI_Allreduce(MPI_IN_PLACE, &num_vertices, 1, KAGEN_MPI_SINT, MPI_MAX, MPI_COMM_WORLD);

    Edgelist redistributed_edges;
    RedistributeEdgesBalanced(input, redistributed_edges, num_vertices, MPI_COMM_WORLD);

    EXPECT_TRUE(std::is_sorted(redistributed_edges.begin(), redistributed_edges.end()));
    auto dup = std::adjacent_find(redistributed_edges.begin(), redistributed_edges.end());
    EXPECT_EQ(dup, redistributed_edges.end());
}

// ---- Star graph tests (skewed degree distribution) -------------------------

TEST(RedistributeEdgesBalanced, PreservesEdgeSet_Star) {
    const SInt n     = 100;
    Edgelist   input = BuildStarOnPE0(n);

    Edgelist reference = input;
    RoundRobinRemapping(reference, n, MPI_COMM_WORLD);
    Edgelist expected = GatherAllEdges(reference);
    std::sort(expected.begin(), expected.end());
    expected.erase(std::unique(expected.begin(), expected.end()), expected.end());

    Edgelist redistributed_edges;
    RedistributeEdgesBalanced(input, redistributed_edges, n, MPI_COMM_WORLD);

    Edgelist result = GatherAllEdges(redistributed_edges);
    std::sort(result.begin(), result.end());

    EXPECT_EQ(expected, result);
}

TEST(RedistributeEdgesBalanced, OwnershipInvariant_Star) {
    const SInt  n     = 100;
    Edgelist    input = BuildStarOnPE0(n);
    Edgelist    redistributed_edges;
    VertexRange vr = RedistributeEdgesBalanced(input, redistributed_edges, n, MPI_COMM_WORLD);

    for (const auto& [u, v]: redistributed_edges) {
        EXPECT_GE(u, vr.first);
        EXPECT_LT(u, vr.second);
    }
}

// ---- Edge cases ------------------------------------------------------------

TEST(RedistributeEdgesBalanced, EmptyInput) {
    const SInt n = 100;
    Edgelist   input;
    Edgelist   redistributed_edges;

    RedistributeEdgesBalanced(input, redistributed_edges, n, MPI_COMM_WORLD);

    EXPECT_TRUE(redistributed_edges.empty());
}

TEST(RedistributeEdgesBalanced, SingleEdge) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const SInt n = 10;
    Edgelist   input;
    if (rank == 0) {
        input.emplace_back(2, 7);
    }

    Edgelist reference = input;
    RoundRobinRemapping(reference, n, MPI_COMM_WORLD);
    Edgelist expected = GatherAllEdges(reference);
    std::sort(expected.begin(), expected.end());

    Edgelist    redistributed_edges;
    VertexRange vr = RedistributeEdgesBalanced(input, redistributed_edges, n, MPI_COMM_WORLD);

    Edgelist result = GatherAllEdges(redistributed_edges);
    std::sort(result.begin(), result.end());

    EXPECT_EQ(expected, result);

    // Ownership invariant
    for (const auto& [u, v]: redistributed_edges) {
        EXPECT_GE(u, vr.first);
        EXPECT_LT(u, vr.second);
    }
}
