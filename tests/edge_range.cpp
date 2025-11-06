#include "kagen/edge_range.h"

#include "kagen/kagen.h"

#include <gtest/gtest.h>

#include "tests/gather.h"
#include "tests/utils.h"
#include "tools/converter.h"

using namespace kagen;

void check_edge_range(const Graph& graph) {
    Edgelist edgelist = graph.edges;
    if (graph.representation == GraphRepresentation::CSR) {
        edgelist = BuildEdgeListFromCSR(graph.vertex_range, graph.xadj, graph.adjncy);
    }
    EdgeRange edge_range(graph);

    {
        std::size_t expected_index = 0;
        for (auto it = edge_range.begin(); it != edge_range.end(); ++it) {
            auto edge = *it;
            EXPECT_EQ(it.edge_index(), expected_index);
            EXPECT_EQ(*it, edge);
            ++expected_index;
        }
    }

    {
        EXPECT_EQ(edge_range.size(), edgelist.size());
        for (std::size_t i = 0; auto elem: edge_range) {
            EXPECT_EQ(elem, edgelist[i]);
            ++i;
        }
    }
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
