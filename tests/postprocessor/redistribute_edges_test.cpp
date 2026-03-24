#ifndef KAGEN_NOMPI
    #include "kagen/comm/mpi_comm.h"
#else
    #include "kagen/comm/seq_comm.h"
#endif
#include "kagen/kagen.h"
#include "kagen/tools/postprocessor.h"

#include <gtest/gtest.h>

#include "../gather.h"
#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

using namespace kagen;

using GeneratorFunc    = std::function<Graph(KaGen&, SInt, SInt)>;
using RedistributeFunc = std::function<VertexRange(Edgelist&, Edgelist&, SInt, bool, kagen::Comm&)>;

static Graph MakeEdgeListGraph(const Edgelist& edges) {
    Graph g;
    g.representation = GraphRepresentation::EDGE_LIST;
    g.edges          = edges;
    return g;
}

static Edgelist GatherAllEdges(const Edgelist& local_edges) {
#ifndef KAGEN_NOMPI
    kagen::MPIComm comm(MPI_COMM_WORLD);
#else
    kagen::SeqComm comm;
#endif
    return kagen::testing::GatherEdgeLists(MakeEdgeListGraph(local_edges), comm).edges;
}

// Build a star graph (vertex 0 connected to all others, both directions) with all edges on PE 0.
static Edgelist BuildStarOnPE0(SInt n) {
#ifndef KAGEN_NOMPI
    kagen::MPIComm mpi_world(MPI_COMM_WORLD);
#else
    kagen::SeqComm mpi_world;
#endif
    int rank = mpi_world.Rank();
    Edgelist edges;
    if (rank == 0) {
        for (SInt v = 1; v < n; ++v) {
            edges.emplace_back(0, v);
            edges.emplace_back(v, 0);
        }
    }
    return edges;
}

// ---- Parameterized tests over generator, redistribution function, and round-robin remapping ---

using GeneratorParam    = std::tuple<std::string, GeneratorFunc>;
using RedistributeParam = std::tuple<std::string, RedistributeFunc>;
using FixtureParam      = std::tuple<GeneratorParam, RedistributeParam, bool>;

struct RedistributeEdgesFixture : public ::testing::TestWithParam<FixtureParam> {};

INSTANTIATE_TEST_SUITE_P(
    RedistributeEdgesTests, RedistributeEdgesFixture,
    ::testing::Combine(
        ::testing::Values(
            std::make_tuple("GNM", GeneratorFunc([](KaGen& gen, SInt n, SInt m) {
                                return gen.GenerateUndirectedGNM(n, m);
                            })),
            std::make_tuple("RMAT", GeneratorFunc([](KaGen& gen, SInt n, SInt m) {
                                return gen.GenerateRMAT(n, m, 0.56, 0.19, 0.19);
                            })),
            std::make_tuple("RGG2D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) {
                                return gen.GenerateRGG2D_NM(n, m);
                            })),
            std::make_tuple("RHG", GeneratorFunc([](KaGen& gen, SInt n, SInt m) {
                                return gen.GenerateRHG_NM(2.8, n, m);
                            })),
            std::make_tuple("Grid2D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) {
                                return gen.GenerateGrid2D_N(n, 1.0);
                            }))),
        ::testing::Values(
            std::make_tuple("EdgeBalanced", RedistributeFunc(RedistributeEdgesBalanced)),
            std::make_tuple("VertexBalanced", RedistributeFunc(RedistributeEdges))),
        ::testing::Values(true, false)),
    [](const ::testing::TestParamInfo<FixtureParam>& info) {
        return std::get<0>(std::get<0>(info.param)) + "_" + std::get<0>(std::get<1>(info.param))
               + (std::get<2>(info.param) ? "_remap" : "_noremap");
    });

TEST_P(RedistributeEdgesFixture, PreservesEdgeSet) {
    auto [gen_pair, redist_pair, remap_round_robin] = GetParam();
    auto generate                                   = std::get<1>(gen_pair);
    auto redistribute                               = std::get<1>(redist_pair);

    const SInt n = 1000;
    const SInt m = 4 * n;

#ifndef KAGEN_NOMPI
    kagen::MPIComm comm(MPI_COMM_WORLD);
    kagen::KaGen generator(comm.GetMPIComm());
#else
    kagen::SeqComm comm;
    kagen::KaGen   generator(MPI_COMM_WORLD);
#endif
    generator.UseEdgeListRepresentation();
    Graph graph = generate(generator, n, m);

    Edgelist input        = graph.edges;
    SInt     num_vertices = graph.vertex_range.second;
    comm.Allreduce(kagen::COMM_IN_PLACE, &num_vertices, 1, kagen::CommDatatype::UNSIGNED_LONG_LONG, kagen::CommOp::MAX);

    Edgelist expected = [&](auto edges) {
        if (remap_round_robin) {
            RoundRobinRemapping(edges, num_vertices, comm);
        }
        edges = GatherAllEdges(edges);
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
        return edges;
    }(input);

    Edgelist redistributed_edges;
    redistribute(input, redistributed_edges, num_vertices, remap_round_robin, comm);

    Edgelist result = GatherAllEdges(redistributed_edges);
    std::sort(result.begin(), result.end());

    EXPECT_EQ(expected, result);
}

TEST_P(RedistributeEdgesFixture, OwnershipInvariant) {
    auto [gen_pair, redist_pair, remap_round_robin] = GetParam();
    auto generate                                   = std::get<1>(gen_pair);
    auto [redist_name, redistribute]                = redist_pair;

    const SInt n = 1000;
    const SInt m = 4 * n;

#ifndef KAGEN_NOMPI
    kagen::MPIComm comm(MPI_COMM_WORLD);
    kagen::KaGen generator(comm.GetMPIComm());
#else
    kagen::SeqComm comm;
    kagen::KaGen   generator(MPI_COMM_WORLD);
#endif
    generator.UseEdgeListRepresentation();
    Graph graph = generate(generator, n, m);

    Edgelist input        = graph.edges;
    SInt     num_vertices = graph.vertex_range.second;
    comm.Allreduce(kagen::COMM_IN_PLACE, &num_vertices, 1, kagen::CommDatatype::UNSIGNED_LONG_LONG, kagen::CommOp::MAX);

    Edgelist    redistributed_edges;
    VertexRange vr = redistribute(input, redistributed_edges, num_vertices, remap_round_robin, comm);

    for (const auto& edge: redistributed_edges) {
        EXPECT_GE(edge.first, vr.first);
        EXPECT_LT(edge.first, vr.second);
    }
}

TEST_P(RedistributeEdgesFixture, NoDuplicatesInOutput) {
    auto [gen_pair, redist_pair, remap_round_robin] = GetParam();
    auto generate                                   = std::get<1>(gen_pair);
    auto [redist_name, redistribute]                = redist_pair;

    const SInt n = 1000;
    const SInt m = 4 * n;

#ifndef KAGEN_NOMPI
    kagen::MPIComm comm(MPI_COMM_WORLD);
    kagen::KaGen generator(comm.GetMPIComm());
#else
    kagen::SeqComm comm;
    kagen::KaGen   generator(MPI_COMM_WORLD);
#endif
    generator.UseEdgeListRepresentation();
    Graph graph = generate(generator, n, m);

    Edgelist input        = graph.edges;
    SInt     num_vertices = graph.vertex_range.second;
    comm.Allreduce(kagen::COMM_IN_PLACE, &num_vertices, 1, kagen::CommDatatype::UNSIGNED_LONG_LONG, kagen::CommOp::MAX);

    Edgelist redistributed_edges;
    redistribute(input, redistributed_edges, num_vertices, remap_round_robin, comm);

    EXPECT_TRUE(std::is_sorted(redistributed_edges.begin(), redistributed_edges.end()));
    auto dup = std::adjacent_find(redistributed_edges.begin(), redistributed_edges.end());
    EXPECT_EQ(dup, redistributed_edges.end());
}

// ---- Star graph and edge case tests -----------------------------------------

using SimpleFixtureParam = std::tuple<RedistributeParam, bool>;

struct RedistributeEdgesSimpleFixture : public ::testing::TestWithParam<SimpleFixtureParam> {};

INSTANTIATE_TEST_SUITE_P(
    RedistributeEdgesSimpleTests, RedistributeEdgesSimpleFixture,
    ::testing::Combine(
        ::testing::Values(
            std::make_tuple("EdgeBalanced", RedistributeFunc(RedistributeEdgesBalanced)),
            std::make_tuple("VertexBalanced", RedistributeFunc(RedistributeEdges))),
        ::testing::Values(true, false)),
    [](const ::testing::TestParamInfo<SimpleFixtureParam>& info) {
        return std::get<0>(std::get<0>(info.param)) + (std::get<1>(info.param) ? "_remap" : "_noremap");
    });

TEST_P(RedistributeEdgesSimpleFixture, PreservesEdgeSet_Star) {
    auto [redist_pair, remap_round_robin] = GetParam();
    auto redistribute                     = std::get<1>(redist_pair);

    const SInt n     = 100;
    Edgelist   input = BuildStarOnPE0(n);

#ifndef KAGEN_NOMPI
    kagen::MPIComm comm(MPI_COMM_WORLD);
#else
    kagen::SeqComm comm;
#endif
    Edgelist reference = input;
    if (remap_round_robin) {
        RoundRobinRemapping(reference, n, comm);
    }
    Edgelist expected = GatherAllEdges(reference);
    std::sort(expected.begin(), expected.end());
    expected.erase(std::unique(expected.begin(), expected.end()), expected.end());

    Edgelist redistributed_edges;
    redistribute(input, redistributed_edges, n, remap_round_robin, comm);

    Edgelist result = GatherAllEdges(redistributed_edges);
    std::sort(result.begin(), result.end());

    EXPECT_EQ(expected, result);
}

TEST_P(RedistributeEdgesSimpleFixture, OwnershipInvariant_Star) {
    auto [redist_pair, remap_round_robin] = GetParam();
    auto redistribute                     = std::get<1>(redist_pair);

    const SInt  n     = 100;
    Edgelist    input = BuildStarOnPE0(n);
    Edgelist    redistributed_edges;
#ifndef KAGEN_NOMPI
    kagen::MPIComm comm(MPI_COMM_WORLD);
#else
    kagen::SeqComm comm;
#endif
    VertexRange vr = redistribute(input, redistributed_edges, n, remap_round_robin, comm);

    for (const auto& edge: redistributed_edges) {
        EXPECT_GE(edge.first, vr.first);
        EXPECT_LT(edge.first, vr.second);
    }
}

TEST_P(RedistributeEdgesSimpleFixture, EmptyInput) {
    auto [redist_pair, remap_round_robin] = GetParam();
    auto redistribute                     = std::get<1>(redist_pair);

    const SInt n = 100;
    Edgelist   input;
    Edgelist   redistributed_edges;

#ifndef KAGEN_NOMPI
    kagen::MPIComm comm(MPI_COMM_WORLD);
#else
    kagen::SeqComm comm;
#endif
    redistribute(input, redistributed_edges, n, remap_round_robin, comm);

    EXPECT_TRUE(redistributed_edges.empty());
}

TEST_P(RedistributeEdgesSimpleFixture, SingleEdge) {
    auto [redist_pair, remap_round_robin] = GetParam();
    auto redistribute                     = std::get<1>(redist_pair);

#ifndef KAGEN_NOMPI
    kagen::MPIComm comm(MPI_COMM_WORLD);
#else
    kagen::SeqComm comm;
#endif
    int rank = comm.Rank();

    const SInt n = 10;
    Edgelist   input;
    if (rank == 0) {
        input.emplace_back(2, 7);
    }

    Edgelist reference = input;
    if (remap_round_robin) {
        RoundRobinRemapping(reference, n, comm);
    }
    Edgelist expected = GatherAllEdges(reference);
    std::sort(expected.begin(), expected.end());

    Edgelist    redistributed_edges;
    VertexRange vr = redistribute(input, redistributed_edges, n, remap_round_robin, comm);

    Edgelist result = GatherAllEdges(redistributed_edges);
    std::sort(result.begin(), result.end());

    EXPECT_EQ(expected, result);

    for (const auto& edge: redistributed_edges) {
        EXPECT_GE(edge.first, vr.first);
        EXPECT_LT(edge.first, vr.second);
    }
}
