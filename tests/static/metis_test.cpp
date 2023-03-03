#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "tests/static/common.h"

const char* EMPTY_GRAPH         = "tests/data/static/empty.metis";
const char* REAL_WORLD_GRAPH    = "tests/data/static/144.metis";
const char* EDGE_WEIGHTED_K3    = "tests/data/static/edge_weighted_K3.metis";
const char* EDGE_WEIGHTED_P2    = "tests/data/static/edge_weighted_P2.metis";
const char* VERTEX_WEIGHTED_K3  = "tests/data/static/vertex_weighted_K3.metis";
const char* VERTEX_WEIGHTED_P2  = "tests/data/static/vertex_weighted_P2.metis";
const char* UNWEIGHTED_K3       = "tests/data/static/unweighted_K3.metis";
const char* WEIGHTED_K3         = "tests/data/static/weighted_K3.metis";
const char* WEIGHTED_P2         = "tests/data/static/weighted_P2.metis";
const char* LARGE_WEIGHTS       = "tests/data/static/large_weights.metis";
const char* GRAPH_WITH_COMMENTS = "tests/data/static/with_comments.metis";

using namespace kagen;

namespace {
StaticGraphDistribution distributions[] = {StaticGraphDistribution::BALANCE_VERTICES};
}

TEST(MetisGenerator, reads_empty_graph_edge_list) {
    using namespace kagen::testing;

    for (const auto distribution: distributions) {
        const auto graph =
            ReadStaticGraph(EMPTY_GRAPH, distribution, StaticGraphFormat::METIS, GraphRepresentation::EDGE_LIST);
        ExpectEmptyGraphEdgeList(graph);
        EXPECT_EQ(graph.vertex_range.first, 0);
        EXPECT_EQ(graph.vertex_range.second, 0);
    }
}

TEST(MetisGenerator, reads_empty_graph_csr) {
    using namespace kagen::testing;

    for (const auto distribution: distributions) {
        const auto graph =
            ReadStaticGraph(EMPTY_GRAPH, distribution, StaticGraphFormat::METIS, GraphRepresentation::CSR);
        ExpectEmptyGraphCSR(graph);
        EXPECT_EQ(graph.vertex_range.first, 0);
        EXPECT_EQ(graph.vertex_range.second, 0);
    }
}

TEST(MetisGenerator, ignores_comments_edge_list) {
    using namespace kagen::testing;

    for (const auto distribution: distributions) {
        const auto local_graph = ReadStaticGraph(
            GRAPH_WITH_COMMENTS, distribution, StaticGraphFormat::METIS, GraphRepresentation::EDGE_LIST);

        const auto global_graph = GatherEdgeLists(local_graph);
        EXPECT_EQ(global_graph.edges.size(), 2);
    }
}

TEST(MetisGenerator, ignores_comments_csr) {
    using namespace kagen::testing;

    for (const auto distribution: distributions) {
        const auto local_graph =
            ReadStaticGraph(GRAPH_WITH_COMMENTS, distribution, StaticGraphFormat::METIS, GraphRepresentation::CSR);
        const auto global_graph = GatherCSR(local_graph);

        ASSERT_EQ(global_graph.xadj.size(), 3);
        EXPECT_EQ(global_graph.xadj[0], 0);
        EXPECT_EQ(global_graph.xadj[1], 1);
        EXPECT_EQ(global_graph.xadj[2], 2);

        EXPECT_EQ(global_graph.adjncy.size(), 2);
    }
}

TEST(MetisGenerator, loads_unweighted_K3_csr) {
    using namespace kagen::testing;

    for (const auto distribution: distributions) {
        const auto local_graph =
            ReadStaticGraph(UNWEIGHTED_K3, distribution, StaticGraphFormat::METIS, GraphRepresentation::CSR);
        const auto global_graph = GatherCSR(local_graph);
        ExpectK3CSR(global_graph);

        EXPECT_TRUE(global_graph.vertex_weights.empty());
        EXPECT_TRUE(global_graph.edge_weights.empty());
    }
}

TEST(MetisGenerator, load_edge_weighted_K3_csr) {
    using namespace kagen::testing;

    for (const auto distribution: distributions) {
        const auto local_graph =
            ReadStaticGraph(EDGE_WEIGHTED_K3, distribution, StaticGraphFormat::METIS, GraphRepresentation::CSR);
        const auto global_graph = GatherCSR(local_graph);
        ExpectK3CSR(global_graph);

        EXPECT_TRUE(global_graph.vertex_weights.empty());
        EXPECT_FALSE(global_graph.edge_weights.empty());
    }
}

TEST(MetisGenerator, load_vertex_weighted_K3_csr) {
    using namespace kagen::testing;

    for (const auto distribution: distributions) {
        const auto local_graph =
            ReadStaticGraph(VERTEX_WEIGHTED_K3, distribution, StaticGraphFormat::METIS, GraphRepresentation::CSR);
        const auto global_graph = GatherCSR(local_graph);
        ExpectK3CSR(global_graph);

        EXPECT_FALSE(global_graph.vertex_weights.empty());
        EXPECT_TRUE(global_graph.edge_weights.empty());
    }
}

TEST(MetisGenerator, load_weighted_K3_csr) {
    using namespace kagen::testing;

    for (const auto distribution: distributions) {
        const auto local_graph =
            ReadStaticGraph(WEIGHTED_K3, distribution, StaticGraphFormat::METIS, GraphRepresentation::CSR);
        const auto global_graph = GatherCSR(local_graph);
        ExpectK3CSR(global_graph);

        EXPECT_FALSE(global_graph.vertex_weights.empty());
        EXPECT_FALSE(global_graph.edge_weights.empty());
    }
}

TEST(MetisGenerator, load_vertex_weighted_P2_csr) {
    using namespace kagen::testing;

    for (const auto distribution: distributions) {
        const auto local_graph =
            ReadStaticGraph(VERTEX_WEIGHTED_P2, distribution, StaticGraphFormat::METIS, GraphRepresentation::CSR);
        const auto global_graph = GatherCSR(local_graph);
        ExpectP2CSR(global_graph);

        EXPECT_FALSE(global_graph.vertex_weights.empty());
        EXPECT_TRUE(global_graph.edge_weights.empty());
    }
}

TEST(MetisGenerator, load_edge_weighted_P2_csr) {
    using namespace kagen::testing;

    for (const auto distribution: distributions) {
        const auto local_graph =
            ReadStaticGraph(EDGE_WEIGHTED_P2, distribution, StaticGraphFormat::METIS, GraphRepresentation::CSR);
        const auto global_graph = GatherCSR(local_graph);
        ExpectP2CSR(global_graph);

        EXPECT_TRUE(global_graph.vertex_weights.empty());
        EXPECT_FALSE(global_graph.edge_weights.empty());
    }
}

TEST(MetisGenerator, load_weighted_P2_csr) {
    using namespace kagen::testing;

    for (const auto distribution: distributions) {
        const auto local_graph =
            ReadStaticGraph(WEIGHTED_P2, distribution, StaticGraphFormat::METIS, GraphRepresentation::CSR);
        const auto global_graph = GatherCSR(local_graph);
        ExpectP2CSR(global_graph);

        EXPECT_FALSE(global_graph.vertex_weights.empty());
        EXPECT_FALSE(global_graph.edge_weights.empty());
    }
}

TEST(MetisGenerator, load_large_weights_csr) {
    using namespace kagen::testing;

    for (const auto distribution: distributions) {
        const auto local_graph =
            ReadStaticGraph(LARGE_WEIGHTS, distribution, StaticGraphFormat::METIS, GraphRepresentation::CSR);
        const auto global_graph = GatherCSR(local_graph);

        ASSERT_EQ(global_graph.xadj.size(), 3);
        ASSERT_EQ(global_graph.adjncy.size(), 2);

        EXPECT_FALSE(global_graph.vertex_weights.empty());
        EXPECT_TRUE(global_graph.edge_weights.empty());

        EXPECT_EQ(global_graph.vertex_weights[0], 123456789);
        EXPECT_EQ(global_graph.vertex_weights[1], 234567891);
    }
}

TEST(MetisGenerator, load_real_world_graph_vertex_balanced_csr) {
    using namespace kagen::testing;

    const auto local_graph = ReadStaticGraph(
        REAL_WORLD_GRAPH, StaticGraphDistribution::BALANCE_VERTICES, StaticGraphFormat::METIS,
        GraphRepresentation::CSR);
    const auto global_graph = GatherCSR(local_graph);

    PEID size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const SInt global_n        = global_graph.xadj.size() - 1;
    const SInt vertices_per_pe = std::ceil(1.0 * global_n / size);

    EXPECT_LE(local_graph.xadj.size() - 1, vertices_per_pe);
    EXPECT_EQ(local_graph.xadj.size() - 1, local_graph.vertex_range.second - local_graph.vertex_range.first);
}

TEST(MetisGenerator, load_real_world_graph_csr) {
    using namespace kagen::testing;
    using namespace ::testing;

    PEID rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (const auto distribution: distributions) {
        const auto local_graph =
            ReadStaticGraph(REAL_WORLD_GRAPH, distribution, StaticGraphFormat::METIS, GraphRepresentation::CSR);
        const auto global_graph = GatherCSR(local_graph);

        ASSERT_EQ(global_graph.xadj.size() - 1, 144649);
        ASSERT_EQ(global_graph.adjncy.size(), 2 * 1074393);
        EXPECT_TRUE(global_graph.vertex_weights.empty());
        EXPECT_TRUE(global_graph.edge_weights.empty());

        // Compare graph against graph read on just a single PE
        if (rank == 0) {
            const auto root_graph = ReadStaticGraphOnRoot(
                REAL_WORLD_GRAPH, distribution, StaticGraphFormat::METIS, GraphRepresentation::CSR);

            ASSERT_EQ(global_graph.xadj.size(), root_graph.xadj.size());
            ASSERT_EQ(global_graph.adjncy.size(), root_graph.adjncy.size());

            bool              ok = true;
            std::vector<SInt> global_neighbors;
            std::vector<SInt> root_neighbors;

            for (SInt u = 0; u + 1 < global_graph.xadj.size(); ++u) {
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
        }

        // Check first and last node by hand
        const SInt first_degree = global_graph.xadj[1] - global_graph.xadj[0];
        ASSERT_EQ(first_degree, 18);
        std::vector<SInt> first_neighbors(global_graph.adjncy.begin(), global_graph.adjncy.begin() + 18);
        ASSERT_THAT(
            first_neighbors, UnorderedElementsAre(
                                 104774, 6, 86994, 132495, 101, 200, 131737, 82652, 142239, 105301, 327, 37, 115494, 70,
                                 264, 131468, 95890, 137));

        const SInt last_degree = global_graph.xadj.back() - global_graph.xadj[global_graph.xadj.size() - 2];
        ASSERT_EQ(last_degree, 12);
        std::vector<SInt> last_neighbors(global_graph.adjncy.end() - last_degree, global_graph.adjncy.end());
        ASSERT_THAT(
            last_neighbors,
            UnorderedElementsAre(50177, 90587, 50681, 27591, 27593, 96339, 50315, 97561, 50835, 99400, 115158, 120458));
    }
}
