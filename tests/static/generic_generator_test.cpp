#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "tests/static/common.h"

const char* EMPTY_GRAPH         = "tests/data/static/empty";
const char* REAL_WORLD_GRAPH    = "tests/data/static/144";
const char* EDGE_WEIGHTED_K3    = "tests/data/static/edge_weighted_K3";
const char* EDGE_WEIGHTED_P2    = "tests/data/static/edge_weighted_P2";
const char* VERTEX_WEIGHTED_K3  = "tests/data/static/vertex_weighted_K3";
const char* VERTEX_WEIGHTED_P2  = "tests/data/static/vertex_weighted_P2";
const char* UNWEIGHTED_K3       = "tests/data/static/unweighted_K3";
const char* WEIGHTED_K3         = "tests/data/static/weighted_K3";
const char* WEIGHTED_P2         = "tests/data/static/weighted_P2";
const char* LARGE_WEIGHTS       = "tests/data/static/large_weights";
const char* GRAPH_WITH_COMMENTS = "tests/data/static/with_comments";

using namespace kagen;

struct GenericGeneratorTestFixture
    : public ::testing::TestWithParam<std::tuple<StaticGraphFormat, StaticGraphDistribution, GraphRepresentation>> {};

INSTANTIATE_TEST_SUITE_P(
    GenericGeneratorTest, GenericGeneratorTestFixture,
    ::testing::Values(
        std::make_tuple(
            StaticGraphFormat::METIS, StaticGraphDistribution::BALANCE_VERTICES, GraphRepresentation::EDGE_LIST),
        std::make_tuple(StaticGraphFormat::METIS, StaticGraphDistribution::BALANCE_VERTICES, GraphRepresentation::CSR),
        std::make_tuple(
            StaticGraphFormat::METIS, StaticGraphDistribution::BALANCE_EDGES, GraphRepresentation::EDGE_LIST)));

TEST_P(GenericGeneratorTestFixture, reads_empty_graph) {
    using namespace kagen::testing;
    const auto [format, distribution, representation] = GetParam();

    const auto graph = ReadStaticGraph(EMPTY_GRAPH, distribution, format, representation);
    ExpectEmptyGraph(graph);
    EXPECT_EQ(graph.vertex_range.first, 0);
    EXPECT_EQ(graph.vertex_range.second, 0);
}

TEST_P(GenericGeneratorTestFixture, ignores_comments) {
    using namespace kagen::testing;
    const auto [format, distribution, representation] = GetParam();

    const auto local_graph  = ReadStaticGraph(GRAPH_WITH_COMMENTS, distribution, format, representation);
    const auto global_graph = GatherGraph(local_graph);

    switch (representation) {
        case GraphRepresentation::CSR:
            EXPECT_EQ(global_graph.xadj.size(), 3);
            EXPECT_EQ(global_graph.adjncy.size(), 2);
            break;

        case GraphRepresentation::EDGE_LIST:
            EXPECT_EQ(global_graph.edges.size(), 2);
            break;
    }
}

TEST_P(GenericGeneratorTestFixture, loads_unweighted_K3) {
    using namespace kagen::testing;
    const auto [format, distribution, representation] = GetParam();

    const auto local_graph  = ReadStaticGraph(UNWEIGHTED_K3, distribution, format, representation);
    const auto global_graph = GatherGraph(local_graph);
    ExpectK3(global_graph);
    EXPECT_TRUE(global_graph.vertex_weights.empty());
    EXPECT_TRUE(global_graph.edge_weights.empty());
}

TEST_P(GenericGeneratorTestFixture, loads_edge_weighted_K3) {
    using namespace kagen::testing;
    const auto [format, distribution, representation] = GetParam();

    const auto local_graph  = ReadStaticGraph(EDGE_WEIGHTED_K3, distribution, format, representation);
    const auto global_graph = GatherGraph(local_graph);
    ExpectK3(global_graph);
    EXPECT_TRUE(global_graph.vertex_weights.empty());
    EXPECT_FALSE(global_graph.edge_weights.empty());
}

TEST_P(GenericGeneratorTestFixture, loads_vertex_weighted_K3) {
    using namespace kagen::testing;
    const auto [format, distribution, representation] = GetParam();

    const auto local_graph  = ReadStaticGraph(VERTEX_WEIGHTED_K3, distribution, format, representation);
    const auto global_graph = GatherGraph(local_graph);
    ExpectK3(global_graph);
    EXPECT_FALSE(global_graph.vertex_weights.empty());
    EXPECT_TRUE(global_graph.edge_weights.empty());
}

TEST_P(GenericGeneratorTestFixture, loads_weighted_K3) {
    using namespace kagen::testing;
    const auto [format, distribution, representation] = GetParam();

    const auto local_graph  = ReadStaticGraph(WEIGHTED_K3, distribution, format, representation);
    const auto global_graph = GatherGraph(local_graph);
    ExpectK3(global_graph);
    EXPECT_FALSE(global_graph.vertex_weights.empty());
    EXPECT_FALSE(global_graph.edge_weights.empty());
}

TEST_P(GenericGeneratorTestFixture, loads_vertex_weighted_P2) {
    using namespace kagen::testing;
    const auto [format, distribution, representation] = GetParam();

    const auto local_graph  = ReadStaticGraph(VERTEX_WEIGHTED_P2, distribution, format, representation);
    const auto global_graph = GatherGraph(local_graph);
    ExpectP2(global_graph);
    EXPECT_FALSE(global_graph.vertex_weights.empty());
    EXPECT_TRUE(global_graph.edge_weights.empty());
}

TEST_P(GenericGeneratorTestFixture, loads_edge_weighted_P2) {
    using namespace kagen::testing;
    const auto [format, distribution, representation] = GetParam();

    const auto local_graph  = ReadStaticGraph(EDGE_WEIGHTED_P2, distribution, format, representation);
    const auto global_graph = GatherGraph(local_graph);
    ExpectP2(global_graph);
    EXPECT_TRUE(global_graph.vertex_weights.empty());
    EXPECT_FALSE(global_graph.edge_weights.empty());
}

TEST_P(GenericGeneratorTestFixture, loads_weighted_P2) {
    using namespace kagen::testing;
    const auto [format, distribution, representation] = GetParam();

    const auto local_graph  = ReadStaticGraph(WEIGHTED_P2, distribution, format, representation);
    const auto global_graph = GatherGraph(local_graph);
    ExpectP2(global_graph);
    EXPECT_FALSE(global_graph.vertex_weights.empty());
    EXPECT_FALSE(global_graph.edge_weights.empty());
}

TEST_P(GenericGeneratorTestFixture, loads_large_weights) {
    using namespace kagen::testing;
    const auto [format, distribution, representation] = GetParam();

    const auto local_graph  = ReadStaticGraph(LARGE_WEIGHTS, distribution, format, representation);
    const auto global_graph = GatherGraph(local_graph);

    switch (representation) {
        case GraphRepresentation::CSR:
            ASSERT_EQ(global_graph.xadj.size(), 3);
            ASSERT_EQ(global_graph.adjncy.size(), 2);
            break;

        case GraphRepresentation::EDGE_LIST:
            ASSERT_EQ(global_graph.edges.size(), 2);
            break;
    }

    EXPECT_FALSE(global_graph.vertex_weights.empty());
    EXPECT_TRUE(global_graph.edge_weights.empty());
    EXPECT_EQ(global_graph.vertex_weights[0], 123456789);
    EXPECT_EQ(global_graph.vertex_weights[1], 234567891);
}

TEST_P(GenericGeneratorTestFixture, loads_real_world_graph) {
    using namespace kagen::testing;
    using namespace ::testing;
    const auto [format, distribution, representation] = GetParam();

    PEID rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const auto local_graph  = ReadStaticGraph(REAL_WORLD_GRAPH, distribution, format, representation);
    auto       global_graph = GatherGraph(local_graph);

    switch (representation) {
        case GraphRepresentation::CSR:
            ASSERT_EQ(global_graph.xadj.size() - 1, 144649);
            ASSERT_EQ(global_graph.adjncy.size(), 2 * 1074393);
            break;

        case GraphRepresentation::EDGE_LIST:
            ASSERT_EQ(global_graph.edges.size(), 2 * 1074393);
            break;
    }

    EXPECT_TRUE(global_graph.vertex_weights.empty());
    EXPECT_TRUE(global_graph.edge_weights.empty());

    // Compare graph against graph read on just a single PE
    if (rank == 0) {
        auto root_graph = ReadStaticGraphOnRoot(REAL_WORLD_GRAPH, distribution, format, representation);

        ASSERT_EQ(global_graph.xadj.size(), root_graph.xadj.size());
        ASSERT_EQ(global_graph.adjncy.size(), root_graph.adjncy.size());
        ASSERT_EQ(global_graph.edges.size(), root_graph.edges.size());

        // Check every 100nth node in CSR representation
        bool              ok = true;
        std::vector<SInt> global_neighbors;
        std::vector<SInt> root_neighbors;

        for (SInt u = 0; u + 1 < global_graph.xadj.size(); u += 100) {
            const SInt global_degree = global_graph.xadj[u + 1] - global_graph.xadj[u];
            const SInt root_degree   = root_graph.xadj[u + 1] - root_graph.xadj[u];
            if (global_degree != root_degree) {
                ok = false;
                break;
            }

            global_neighbors.clear();
            root_neighbors.clear();
            for (SInt e = global_graph.xadj[u]; e < global_graph.xadj[u + 1]; ++e) {
                global_neighbors.push_back(global_graph.adjncy[e]);
                root_neighbors.push_back(root_graph.adjncy[e]);
            }

            std::sort(global_neighbors.begin(), global_neighbors.end());
            std::sort(root_neighbors.begin(), root_neighbors.end());
            if (global_neighbors != root_neighbors) {
                ok = false;
                break;
            }
        }
        ASSERT_TRUE(ok);

        // Check edge list, TODO reduce running time
        std::sort(global_graph.edges.begin(), global_graph.edges.end());
        std::sort(root_graph.edges.begin(), root_graph.edges.end());
        ASSERT_EQ(global_graph.edges, root_graph.edges);
    }

    // Check first and last node by hand
    std::vector<SInt> first_neighbors;
    std::vector<SInt> last_neighbors;

    switch (representation) {
        case GraphRepresentation::CSR: {
            const SInt first_degree = global_graph.xadj[1] - global_graph.xadj[0];
            ASSERT_EQ(first_degree, 18);
            first_neighbors = std::vector<SInt>(global_graph.adjncy.begin(), global_graph.adjncy.begin() + 18);

            const SInt last_degree = global_graph.xadj.back() - global_graph.xadj[global_graph.xadj.size() - 2];
            ASSERT_EQ(last_degree, 12);
            last_neighbors = std::vector<SInt>(global_graph.adjncy.end() - last_degree, global_graph.adjncy.end());
        } break;

        case GraphRepresentation::EDGE_LIST: {
            for (const auto& [u, v]: global_graph.edges) {
                if (u == 0) {
                    first_neighbors.push_back(v);
                } else if (u == 144648) {
                    last_neighbors.push_back(v);
                }
            }
        } break;
    }

    ASSERT_THAT(
        first_neighbors, UnorderedElementsAre(
                             104774, 6, 86994, 132495, 101, 200, 131737, 82652, 142239, 105301, 327, 37, 115494, 70,
                             264, 131468, 95890, 137));
    ASSERT_THAT(
        last_neighbors,
        UnorderedElementsAre(50177, 90587, 50681, 27591, 27593, 96339, 50315, 97561, 50835, 99400, 115158, 120458));
}
