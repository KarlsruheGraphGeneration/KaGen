#include "kagen/edge_range.h"

#include "kagen/kagen.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "tests/gather.h"
#include "tests/utils.h"
#include "tools/converter.h"
#include "tools/geometry.h"

using namespace kagen;

TEST(EdgeRangeTest, iterate_edgelist_representation) {
    using ::testing::ElementsAreArray;

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    const SInt n = 1000;
    const SInt m = 16 * n;

    auto check = [](const Graph& graph) {
        Edgelist expected = graph.edges;
        EdgeRange edge_range(graph);

        // Check edges match and indices are consecutive
        EXPECT_THAT(std::vector(edge_range.begin(), edge_range.end()), ElementsAreArray(expected));
        
        std::size_t idx = 0;
        for (auto it = edge_range.begin(); it != edge_range.end(); ++it, ++idx) {
            EXPECT_EQ(it.edge_index(), idx);
        }
    };

    check(generator.GenerateUndirectedGNM(n, m));
    check(generator.GenerateRMAT(n, m, 0.56, 0.19, 0.19));
    check(generator.GenerateRGG2D_NM(n, m));
    check(generator.GenerateRGG3D_NM(n, m));
    check(generator.GenerateRHG_NM(2.6, n, m));
    check(generator.GenerateGrid2D_NM(n, m));
    check(generator.GenerateGrid3D_NM(n, m));
}

TEST(EdgeRangeTest, iterate_sparse_edgelist_representation) {
    using ::testing::ElementsAreArray;

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    const SInt n = 1000;
    const SInt m = 2 * n;

    auto check = [](const Graph& graph) {
        Edgelist expected = graph.edges;
        EdgeRange edge_range(graph);

        // Check edges match and indices are consecutive
        EXPECT_THAT(std::vector(edge_range.begin(), edge_range.end()), ElementsAreArray(expected));
        
        std::size_t idx = 0;
        for (auto it = edge_range.begin(); it != edge_range.end(); ++it, ++idx) {
            EXPECT_EQ(it.edge_index(), idx);
        }
    };

    check(generator.GenerateUndirectedGNM(n, m));
    check(generator.GenerateRMAT(n, m, 0.56, 0.19, 0.19));
    check(generator.GenerateRGG2D_NM(n, m));
    check(generator.GenerateRGG3D_NM(n, m));
    check(generator.GenerateRHG_NM(2.6, n, m));
    check(generator.GenerateGrid2D_NM(n, m));
    check(generator.GenerateGrid3D_NM(n, m));
}

TEST(EdgeRangeTest, iterate_csr_representation) {
    using ::testing::ElementsAreArray;

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseCSRRepresentation();
    const SInt n = 1000;
    const SInt m = 16 * n;

    auto check = [](const Graph& graph) {
        Edgelist expected = graph.edges;
        if (graph.representation == GraphRepresentation::CSR) {
            expected = BuildEdgeListFromCSR(graph.vertex_range, graph.xadj, graph.adjncy);
        }
        EdgeRange edge_range(graph);

        // Check edges match and indices are consecutive
        EXPECT_THAT(std::vector(edge_range.begin(), edge_range.end()), ElementsAreArray(expected));
        
        std::size_t idx = 0;
        for (auto it = edge_range.begin(); it != edge_range.end(); ++it, ++idx) {
            EXPECT_EQ(it.edge_index(), idx);
        }
    };

    check(generator.GenerateUndirectedGNM(n, m));
    check(generator.GenerateRMAT(n, m, 0.56, 0.19, 0.19));
    check(generator.GenerateRGG2D_NM(n, m));
    check(generator.GenerateRGG3D_NM(n, m));
    check(generator.GenerateRHG_NM(2.6, n, m));
    check(generator.GenerateGrid2D_NM(n, m));
    check(generator.GenerateGrid3D_NM(n, m));
}

TEST(EdgeRangeTest, iterate_sparse_csr_representation) {
    using ::testing::ElementsAreArray;

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseCSRRepresentation();
    const SInt n = 1000;
    const SInt m = 2 * n;

    auto check = [](const Graph& graph) {
        Edgelist expected = graph.edges;
        if (graph.representation == GraphRepresentation::CSR) {
            expected = BuildEdgeListFromCSR(graph.vertex_range, graph.xadj, graph.adjncy);
        }
        EdgeRange edge_range(graph);

        // Check edges match and indices are consecutive
        EXPECT_THAT(std::vector(edge_range.begin(), edge_range.end()), ElementsAreArray(expected));
        
        std::size_t idx = 0;
        for (auto it = edge_range.begin(); it != edge_range.end(); ++it, ++idx) {
            EXPECT_EQ(it.edge_index(), idx);
        }
    };

    check(generator.GenerateUndirectedGNM(n, m));
    check(generator.GenerateRMAT(n, m, 0.56, 0.19, 0.19));
    check(generator.GenerateRGG2D_NM(n, m));
    check(generator.GenerateRGG3D_NM(n, m));
    check(generator.GenerateRHG_NM(2.6, n, m));
    check(generator.GenerateGrid2D_NM(n, m));
    check(generator.GenerateGrid3D_NM(n, m));
}
