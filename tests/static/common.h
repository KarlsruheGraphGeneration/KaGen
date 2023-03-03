#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "kagen/context.h"
#include "kagen/generators/static/static_graph.h"

#include <numeric>

#include <mpi.h>

namespace kagen::testing {
inline Graph ReadStaticGraph(
    const char* filename, const StaticGraphDistribution distribution, const StaticGraphFormat format,
    const GraphRepresentation representation) {
    PGeneratorConfig config;
    config.static_graph.filename     = filename;
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
    const char* filename, const StaticGraphDistribution distribution, const StaticGraphFormat format,
    const GraphRepresentation representation) {
    PEID rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        PGeneratorConfig config;
        config.static_graph.filename     = filename;
        config.static_graph.distribution = distribution;
        config.static_graph.format       = format;

        StaticGraph generator(config, 0, 1);
        generator.Generate(representation);
        return generator.Take();
    } else {
        return {};
    }
}

inline void ExpectEmptyGraphCSR(const Graph& graph) {
    EXPECT_EQ(graph.xadj.size(), 1);
    EXPECT_TRUE(graph.adjncy.empty());
    EXPECT_TRUE(graph.vertex_weights.empty());
    EXPECT_TRUE(graph.edge_weights.empty());
}

inline void ExpectEmptyGraphEdgeList(const Graph& graph) {
    EXPECT_TRUE(graph.edges.empty());
    EXPECT_TRUE(graph.vertex_weights.empty());
    EXPECT_TRUE(graph.edge_weights.empty());
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
} // namespace kagen::testing
