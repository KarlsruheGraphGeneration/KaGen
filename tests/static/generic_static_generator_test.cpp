#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <mpi.h>

#include <numeric>
#include <utility>

#include "kagen/context.h"
#include "kagen/generators/static/static_graph.h"

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
            StaticGraphFormat::METIS, StaticGraphDistribution::BALANCE_EDGES, GraphRepresentation::EDGE_LIST),
        std::make_tuple(StaticGraphFormat::METIS, StaticGraphDistribution::BALANCE_EDGES, GraphRepresentation::CSR),
        std::make_tuple(
            StaticGraphFormat::BINARY_PARHIP, StaticGraphDistribution::BALANCE_VERTICES,
            GraphRepresentation::EDGE_LIST)));

namespace {
inline Graph ReadStaticGraph(
    const std::string& filename, const StaticGraphDistribution distribution, const StaticGraphFormat format,
    const GraphRepresentation representation) {
    PGeneratorConfig config;

    switch (format) {
        case StaticGraphFormat::METIS:
            config.static_graph.filename = filename + ".metis";
            break;

        case StaticGraphFormat::BINARY_PARHIP:
            config.static_graph.filename = filename + ".bgf";
            break;
    }

    config.static_graph.distribution = distribution;
    config.static_graph.format       = format;

    PEID size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    StaticGraph generator(config, rank, size);
    generator.Generate(representation);
    return generator.Take();
}

inline Graph ReadStaticGraphOnRoot(
    const std::string& filename, const StaticGraphDistribution distribution, const StaticGraphFormat format,
    const GraphRepresentation representation) {
    PEID rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        PGeneratorConfig config;
        switch (format) {
            case StaticGraphFormat::METIS:
                config.static_graph.filename = filename + ".metis";
                break;

            case StaticGraphFormat::BINARY_PARHIP:
                config.static_graph.filename = filename + ".bgf";
                break;
        }
        config.static_graph.distribution = distribution;
        config.static_graph.format       = format;

        StaticGraph generator(config, 0, 1);
        generator.Generate(representation);
        return generator.Take();
    } else {
        return {};
    }
}

inline void ExpectEmptyGraph(const Graph& graph) {
    switch (graph.representation) {
        case GraphRepresentation::CSR:
            EXPECT_EQ(graph.xadj.size(), 1);
            EXPECT_TRUE(graph.adjncy.empty());
            EXPECT_TRUE(graph.vertex_weights.empty());
            EXPECT_TRUE(graph.edge_weights.empty());
            break;

        case GraphRepresentation::EDGE_LIST:
            EXPECT_TRUE(graph.edges.empty());
            EXPECT_TRUE(graph.vertex_weights.empty());
            EXPECT_TRUE(graph.edge_weights.empty());
            break;
    }
}

inline void ExpectK3CSR(const Graph& graph) {
    using namespace ::testing;

    ASSERT_EQ(graph.xadj.size(), 4);
    EXPECT_EQ(graph.xadj[0], 0);
    EXPECT_EQ(graph.xadj[1], 2);
    EXPECT_EQ(graph.xadj[2], 4);
    EXPECT_EQ(graph.xadj[3], 6);

    ASSERT_EQ(graph.adjncy.size(), 6);
    EXPECT_THAT(graph.adjncy[0], AnyOf(Eq(1), Eq(2)));
    EXPECT_THAT(graph.adjncy[1], AnyOf(Eq(1), Eq(2)));
    EXPECT_NE(graph.adjncy[0], graph.adjncy[1]);
    EXPECT_THAT(graph.adjncy[2], AnyOf(Eq(0), Eq(2)));
    EXPECT_THAT(graph.adjncy[3], AnyOf(Eq(0), Eq(2)));
    EXPECT_NE(graph.adjncy[2], graph.adjncy[3]);
    EXPECT_THAT(graph.adjncy[4], AnyOf(Eq(1), Eq(0)));
    EXPECT_THAT(graph.adjncy[5], AnyOf(Eq(1), Eq(0)));
    EXPECT_NE(graph.adjncy[4], graph.adjncy[5]);
}

inline void ExpectK3EdgeList(const Graph& graph) {
    using namespace ::testing;

    auto e = [](const SInt u, const SInt v) -> std::tuple<SInt, SInt> {
        return std::make_tuple(u, v);
    };

    EXPECT_THAT(graph.edges, UnorderedElementsAre(e(0, 1), e(1, 0), e(0, 2), e(2, 0), e(1, 2), e(2, 1)));
}

inline void ExpectK3(const Graph& graph) {
    switch (graph.representation) {
        case GraphRepresentation::CSR:
            ExpectK3CSR(graph);
            break;

        case GraphRepresentation::EDGE_LIST:
            ExpectK3EdgeList(graph);
            break;
    }
}

inline void ExpectP2CSR(const Graph& graph) {
    using namespace ::testing;

    ASSERT_EQ(graph.xadj.size(), 4);
    EXPECT_EQ(graph.xadj[0], 0);
    EXPECT_EQ(graph.xadj[1], 1);
    EXPECT_EQ(graph.xadj[2], 3);
    EXPECT_EQ(graph.xadj[3], 4);

    ASSERT_EQ(graph.adjncy.size(), 4);
    EXPECT_EQ(graph.adjncy[0], 1);
    EXPECT_THAT(graph.adjncy[1], AnyOf(Eq(0), Eq(2)));
    EXPECT_THAT(graph.adjncy[2], AnyOf(Eq(0), Eq(2)));
    EXPECT_NE(graph.adjncy[1], graph.adjncy[2]);
    EXPECT_EQ(graph.adjncy[3], 1);
}

inline void ExpectP2EdgeList(const Graph& graph) {
    using namespace ::testing;

    auto e = [](const SInt u, const SInt v) -> std::tuple<SInt, SInt> {
        return std::make_tuple(u, v);
    };

    EXPECT_THAT(graph.edges, UnorderedElementsAre(e(0, 1), e(1, 0), e(1, 2), e(2, 1)));
}

inline void ExpectP2(const Graph& graph) {
    switch (graph.representation) {
        case GraphRepresentation::CSR:
            ExpectP2CSR(graph);
            break;

        case GraphRepresentation::EDGE_LIST:
            ExpectP2EdgeList(graph);
            break;
    }
}

inline void GatherWeights(const Graph& local_graph, Graph& global_graph) {
    PEID size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int has_vertex_weights = !local_graph.vertex_weights.empty();
    MPI_Allreduce(MPI_IN_PLACE, &has_vertex_weights, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (has_vertex_weights) {
        const int        num_local_vertices = local_graph.vertex_range.second - local_graph.vertex_range.first;
        std::vector<int> recvcounts(size);
        MPI_Allgather(&num_local_vertices, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
        std::vector<int> displs(size);
        std::exclusive_scan(recvcounts.begin(), recvcounts.end(), displs.begin(), 0);

        global_graph.vertex_weights.resize(recvcounts.back() + displs.back());
        MPI_Allgatherv(
            local_graph.vertex_weights.data(), num_local_vertices, KAGEN_MPI_SSINT, global_graph.vertex_weights.data(),
            recvcounts.data(), displs.data(), KAGEN_MPI_SSINT, MPI_COMM_WORLD);
    }

    int has_edge_weights = !local_graph.edge_weights.empty();
    MPI_Allreduce(MPI_IN_PLACE, &has_edge_weights, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (has_edge_weights) {
        const int        num_local_edges = std::max<int>(local_graph.edges.size(), local_graph.adjncy.size());
        std::vector<int> recvcounts(size);
        MPI_Allgather(&num_local_edges, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
        std::vector<int> displs(size);
        std::exclusive_scan(recvcounts.begin(), recvcounts.end(), displs.begin(), 0);

        global_graph.edge_weights.resize(recvcounts.back() + displs.back());
        MPI_Allgatherv(
            local_graph.edge_weights.data(), num_local_edges, KAGEN_MPI_SSINT, global_graph.edge_weights.data(),
            recvcounts.data(), displs.data(), KAGEN_MPI_SSINT, MPI_COMM_WORLD);
    }
}

inline Graph GatherEdgeLists(const Graph& local_graph) {
    PEID size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int        num_local_edges = local_graph.edges.size() * 2;
    std::vector<int> recvcounts(size);
    MPI_Allgather(&num_local_edges, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<int> displs(size);
    std::exclusive_scan(recvcounts.begin(), recvcounts.end(), displs.begin(), 0u);
    const SInt num_global_edges = displs.back() + recvcounts.back();

    Graph global_graph;
    global_graph.edges.resize(num_global_edges / 2);
    MPI_Allgatherv(
        local_graph.edges.data(), num_local_edges, KAGEN_MPI_SINT, global_graph.edges.data(), recvcounts.data(),
        displs.data(), KAGEN_MPI_SINT, MPI_COMM_WORLD);

    GatherWeights(local_graph, global_graph);
    return global_graph;
}

inline Graph GatherCSR(const Graph& local_graph) {
    EXPECT_EQ(local_graph.xadj.size() - 1, local_graph.vertex_range.second - local_graph.vertex_range.first);
    EXPECT_EQ(local_graph.adjncy.size(), local_graph.xadj.back());

    PEID size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int        num_local_vertices = local_graph.xadj.size() - 1;
    std::vector<int> degree_recvcounts(size);
    MPI_Allgather(&num_local_vertices, 1, MPI_INT, degree_recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    std::vector<int> degree_displs(size);
    std::exclusive_scan(degree_recvcounts.begin(), degree_recvcounts.end(), degree_displs.begin(), 0);

    const int        num_local_edges = local_graph.adjncy.size();
    std::vector<int> edges_recvcounts(size);
    MPI_Allgather(&num_local_edges, 1, MPI_INT, edges_recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    std::vector<int> edges_displs(size);
    std::exclusive_scan(edges_recvcounts.begin(), edges_recvcounts.end(), edges_displs.begin(), 0);

    std::vector<SInt> local_degrees(num_local_vertices);
    for (SInt u = 0; u < num_local_vertices; ++u) {
        local_degrees[u] = local_graph.xadj[u + 1] - local_graph.xadj[u];
    }

    Graph global_graph;
    global_graph.xadj.resize(degree_recvcounts.back() + degree_displs.back() + 1);
    global_graph.adjncy.resize(edges_recvcounts.back() + edges_displs.back());

    MPI_Allgatherv(
        local_degrees.data(), num_local_vertices, KAGEN_MPI_SINT, global_graph.xadj.data(), degree_recvcounts.data(),
        degree_displs.data(), KAGEN_MPI_SINT, MPI_COMM_WORLD);
    MPI_Allgatherv(
        local_graph.adjncy.data(), num_local_edges, KAGEN_MPI_SINT, global_graph.adjncy.data(), edges_recvcounts.data(),
        edges_displs.data(), KAGEN_MPI_SINT, MPI_COMM_WORLD);

    std::exclusive_scan(global_graph.xadj.begin(), global_graph.xadj.end(), global_graph.xadj.begin(), 0);

    EXPECT_EQ(global_graph.xadj.back(), global_graph.adjncy.size());

    GatherWeights(local_graph, global_graph);
    return global_graph;
}

inline Graph GatherGraph(const Graph& local_graph) {
    switch (local_graph.representation) {
        case GraphRepresentation::CSR:
            return GatherCSR(local_graph);

        case GraphRepresentation::EDGE_LIST:
            return GatherEdgeLists(local_graph);
    }

    __builtin_unreachable();
}
} // namespace

TEST_P(GenericGeneratorTestFixture, reads_empty_graph) {
    const auto [format, distribution, representation] = GetParam();

    const auto graph = ReadStaticGraph(EMPTY_GRAPH, distribution, format, representation);
    ExpectEmptyGraph(graph);
    EXPECT_EQ(graph.vertex_range.first, 0);
    EXPECT_EQ(graph.vertex_range.second, 0);
}

TEST_P(GenericGeneratorTestFixture, ignores_comments) {
    const auto [format, distribution, representation] = GetParam();

    const auto local_graph  = ReadStaticGraph(GRAPH_WITH_COMMENTS, distribution, format, representation);
    std::cout << "OK" << std::endl;
    const auto global_graph = GatherGraph(local_graph);
    std::cout << "OK" << std::endl;

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
    const auto [format, distribution, representation] = GetParam();

    const auto local_graph  = ReadStaticGraph(UNWEIGHTED_K3, distribution, format, representation);
    const auto global_graph = GatherGraph(local_graph);
    ExpectK3(global_graph);
    EXPECT_TRUE(global_graph.vertex_weights.empty());
    EXPECT_TRUE(global_graph.edge_weights.empty());
}

TEST_P(GenericGeneratorTestFixture, loads_edge_weighted_K3) {
    const auto [format, distribution, representation] = GetParam();

    const auto local_graph  = ReadStaticGraph(EDGE_WEIGHTED_K3, distribution, format, representation);
    const auto global_graph = GatherGraph(local_graph);
    ExpectK3(global_graph);
    EXPECT_TRUE(global_graph.vertex_weights.empty());
    EXPECT_FALSE(global_graph.edge_weights.empty());
}

TEST_P(GenericGeneratorTestFixture, loads_vertex_weighted_K3) {
    const auto [format, distribution, representation] = GetParam();

    const auto local_graph  = ReadStaticGraph(VERTEX_WEIGHTED_K3, distribution, format, representation);
    const auto global_graph = GatherGraph(local_graph);
    ExpectK3(global_graph);
    EXPECT_FALSE(global_graph.vertex_weights.empty());
    EXPECT_TRUE(global_graph.edge_weights.empty());
}

TEST_P(GenericGeneratorTestFixture, loads_weighted_K3) {
    const auto [format, distribution, representation] = GetParam();

    const auto local_graph  = ReadStaticGraph(WEIGHTED_K3, distribution, format, representation);
    const auto global_graph = GatherGraph(local_graph);
    ExpectK3(global_graph);
    EXPECT_FALSE(global_graph.vertex_weights.empty());
    EXPECT_FALSE(global_graph.edge_weights.empty());
}

TEST_P(GenericGeneratorTestFixture, loads_vertex_weighted_P2) {
    const auto [format, distribution, representation] = GetParam();

    const auto local_graph  = ReadStaticGraph(VERTEX_WEIGHTED_P2, distribution, format, representation);
    const auto global_graph = GatherGraph(local_graph);
    ExpectP2(global_graph);
    EXPECT_FALSE(global_graph.vertex_weights.empty());
    EXPECT_TRUE(global_graph.edge_weights.empty());
}

TEST_P(GenericGeneratorTestFixture, loads_edge_weighted_P2) {
    const auto [format, distribution, representation] = GetParam();

    const auto local_graph  = ReadStaticGraph(EDGE_WEIGHTED_P2, distribution, format, representation);
    const auto global_graph = GatherGraph(local_graph);
    ExpectP2(global_graph);
    EXPECT_TRUE(global_graph.vertex_weights.empty());
    EXPECT_FALSE(global_graph.edge_weights.empty());
}

TEST_P(GenericGeneratorTestFixture, loads_weighted_P2) {
    const auto [format, distribution, representation] = GetParam();

    const auto local_graph  = ReadStaticGraph(WEIGHTED_P2, distribution, format, representation);
    const auto global_graph = GatherGraph(local_graph);
    ExpectP2(global_graph);
    EXPECT_FALSE(global_graph.vertex_weights.empty());
    EXPECT_FALSE(global_graph.edge_weights.empty());
}

TEST_P(GenericGeneratorTestFixture, loads_large_weights) {
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

        std::vector<std::tuple<SInt, SInt>> global_sample, root_sample;
        for (const auto& [u, v]: global_graph.edges) {
            if (u % 100 == 0) {
                global_sample.emplace_back(u, v);
            }
        }
        for (const auto& [u, v]: root_graph.edges) {
            if (u % 100 == 0) {
                root_sample.emplace_back(u, v);
            }
        }

        // Check edge list
        std::sort(global_sample.begin(), global_sample.end());
        std::sort(root_sample.begin(), root_sample.end());
        ASSERT_EQ(global_sample, root_sample);
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
