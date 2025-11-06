#include "kagen/edge_range.h"

#include "kagen/kagen.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <numeric>

#include "tests/gather.h"
#include "tests/utils.h"
#include "tools/converter.h"
#include "tools/geometry.h"

using namespace kagen;

void check_edge_range(const Graph& graph) {
    using ::testing::ElementsAreArray;

    Edgelist edgelist = graph.edges;
    if (graph.representation == GraphRepresentation::CSR) {
        edgelist = BuildEdgeListFromCSR(graph.vertex_range, graph.xadj, graph.adjncy);
    }
    EdgeRange edge_range(graph);

    // Collect edges from edge_range into a vector for comparison
    std::vector<EdgeRange::Edge> edges_from_range(edge_range.begin(), edge_range.end());

    // Check that edge_range produces the same edges as edgelist
    EXPECT_THAT(edges_from_range, ElementsAreArray(edgelist));

    // Collect edge indices to verify they are consecutive starting from 0
    std::vector<std::size_t> edge_indices;
    for (auto it = edge_range.begin(); it != edge_range.end(); ++it) {
        edge_indices.push_back(it.edge_index());
    }

    // Generate expected indices: [0, 1, 2, ..., n-1]
    std::vector<std::size_t> expected_indices(edge_indices.size());
    std::iota(expected_indices.begin(), expected_indices.end(), 0);

    EXPECT_THAT(edge_indices, ElementsAreArray(expected_indices));
}

void check_edge_range(KaGen& generator, SInt n, SInt m) {
    // GNM
    {
        kagen::Graph graph = generator.GenerateUndirectedGNM(n, m);
        check_edge_range(graph);
    }
    // RMAT
    {
        kagen::Graph graph = generator.GenerateRMAT(n, m, 0.56, 0.19, 0.19);
        check_edge_range(graph);
    }
    // RGG2D
    {
        kagen::Graph graph = generator.GenerateRGG2D_NM(n, m);
        check_edge_range(graph);
    }
    // RGG3D
    {
        kagen::Graph graph = generator.GenerateRGG3D_NM(n, m);
        check_edge_range(graph);
    }
    // RHG
    {
        kagen::Graph graph = generator.GenerateRHG_NM(2.6, n, m);
        check_edge_range(graph);
    }
    // GRID2D
    {
        kagen::Graph graph = generator.GenerateGrid2D_NM(n, m);
        check_edge_range(graph);
    }
    // GRID2D
    {
        kagen::Graph graph = generator.GenerateGrid3D_NM(n, m);
        check_edge_range(graph);
    }
}

TEST(EdgeRangeTest, iterate_edgelist_representation) {
    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    const SInt n = 1000;
    const SInt m = 16 * n;
    check_edge_range(generator, n, m);
}

TEST(EdgeRangeTest, iterate_sparse_edgelist_representation) {
    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    const SInt n = 1000;
    const SInt m = 2 * n;
    check_edge_range(generator, n, m);
}

TEST(EdgeRangeTest, iterate_csr_representation) {
    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseCSRRepresentation();
    const SInt n = 1000;
    const SInt m = 16 * n;
    check_edge_range(generator, n, m);
}

TEST(EdgeRangeTest, iterate_sparse_csr_representation) {
    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseCSRRepresentation();
    const SInt n = 1000;
    const SInt m = 2 * n;
    check_edge_range(generator, n, m);
}
