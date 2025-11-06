#include "kagen/edge_range.h"

#include "kagen/kagen.h"

#include <gtest/gtest.h>

#include "tests/gather.h"
#include "tests/utils.h"
#include "tools/converter.h"
#include "tools/geometry.h"

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

void check_bidirectional_iteration(const Graph& graph) {
    Edgelist edgelist = graph.edges;
    if (graph.representation == GraphRepresentation::CSR) {
        edgelist = BuildEdgeListFromCSR(graph.vertex_range, graph.xadj, graph.adjncy);
    }
    EdgeRange edge_range(graph);

    if (edgelist.empty()) {
        return;
    }

    // Forward then backward iteration
    {
        std::vector<EdgeRange::Edge> forward_edges;
        for (auto it = edge_range.begin(); it != edge_range.end(); ++it) {
            forward_edges.push_back(*it);
        }

        std::vector<EdgeRange::Edge> backward_edges;
        auto it = edge_range.end();
        while (it != edge_range.begin()) {
            --it;
            backward_edges.push_back(*it);
        }

        std::reverse(backward_edges.begin(), backward_edges.end());
        EXPECT_EQ(forward_edges.size(), backward_edges.size());
        for (std::size_t i = 0; i < forward_edges.size(); ++i) {
            EXPECT_EQ(forward_edges[i], backward_edges[i]);
        }
    }

    // Mixed forward and backward iteration
    {
        auto it = edge_range.begin();
        ++it;
        ++it;
        auto edge1 = *it;
        --it;
        auto edge2 = *it;
        ++it;
        EXPECT_EQ(*it, edge1);
        EXPECT_EQ(edge2, edgelist[1]);
    }
}

TEST(EdgeRangeTest, bidirectional_iteration_edgelist) {
    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    const SInt n = 100;
    const SInt m = 10 * n;
    
    kagen::Graph graph = generator.GenerateUndirectedGNM(n, m);
    check_bidirectional_iteration(graph);
    
    graph = generator.GenerateRGG2D_NM(n, m);
    check_bidirectional_iteration(graph);
}

TEST(EdgeRangeTest, bidirectional_iteration_csr) {
    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseCSRRepresentation();
    const SInt n = 100;
    const SInt m = 10 * n;
    
    kagen::Graph graph = generator.GenerateUndirectedGNM(n, m);
    check_bidirectional_iteration(graph);
    
    graph = generator.GenerateRGG2D_NM(n, m);
    check_bidirectional_iteration(graph);
}
