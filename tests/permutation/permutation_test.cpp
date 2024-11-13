#include "kagen/kagen.h"

#include <gtest/gtest.h>
#include <kagen/tools/random_permutation.h>

#include "../gather.h"
#include "factories.h"

using namespace kagen;

using SrcDstEdgeWeight = std::tuple<SInt, SInt, SSInt>;
std::vector<SrcDstEdgeWeight> convert_to_weighted_edgelist(const Graph& graph) {
    std::vector<SrcDstEdgeWeight> weighted_edges;
    bool                          is_edge_weighted = !graph.edge_weights.empty();
    switch (graph.representation) {
        case kagen::GraphRepresentation::EDGE_LIST: {
            for (std::size_t i = 0; i < graph.edges.size(); ++i) {
                const auto& [src, dst]  = graph.edges[i];
                const SSInt edge_weight = is_edge_weighted ? graph.edge_weights[i] : 0;
                weighted_edges.emplace_back(src, dst, edge_weight);
            }
        }
        case kagen::GraphRepresentation::CSR: {
            const SInt  global_v_offset = graph.vertex_range.first;
            std::size_t j               = 0;
            for (std::size_t i = 0; i + 1 < graph.xadj.size(); ++i) {
                const SInt src = global_v_offset + i;
                for (; j < graph.xadj[i + 1]; ++j) {
                    const SInt  dst         = graph.adjncy[j];
                    const SSInt edge_weight = is_edge_weighted ? graph.edge_weights[j] : 0;
                    weighted_edges.emplace_back(src, dst, edge_weight);
                }
            }
        };
    }
    return weighted_edges;
}
template <typename Transformer>
void transform_vertices(std::vector<SrcDstEdgeWeight>& weighted_edges, Transformer&& op) {
    for (auto& [src, dst, _]: weighted_edges) {
        src = op(src);
        dst = op(dst);
    }
}

void test_distribution_of_edges(SSInt n_from_config, const Graph& permuted_graph) {
    int      size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &size);

    // vertices should be split equally
    const SInt n_local = permuted_graph.NumberOfLocalVertices();
    EXPECT_LE(n_local, n_from_config / size + 2);

    auto       permuted_graph_edges = convert_to_weighted_edgelist(permuted_graph);
    const auto vertex_range         = permuted_graph.vertex_range;
    for (const auto& [src, _, __]: permuted_graph_edges) {
        EXPECT_GE(src, vertex_range.first);
        EXPECT_LT(src, vertex_range.second);
    }
}

void test_equality_of_permuted_graph(SSInt n_from_config, const Graph& permuted_graph, const Graph& graph) {
    auto gathered_graph                = kagen::testing::GatherGraph(graph);
    auto gathered_graph_edges          = convert_to_weighted_edgelist(gathered_graph);
    auto gathered_permuted_graph       = kagen::testing::GatherGraph(permuted_graph);
    auto gathered_permuted_graph_edges = convert_to_weighted_edgelist(gathered_permuted_graph);

    // we do not want the 'random' permutation to be the ID
    EXPECT_NE(gathered_graph_edges, gathered_permuted_graph_edges);

    // apply inverse permutation to permuted graph
    auto permutator = random_permutation::FeistelPseudoRandomPermutation::buildPermutation(n_from_config - 1, 0);
    auto permute    = [&permutator](SInt v) {
        return permutator.finv(v);
    };
    transform_vertices(gathered_permuted_graph_edges, permute);
    std::sort(gathered_graph_edges.begin(), gathered_graph_edges.end());
    std::sort(gathered_permuted_graph_edges.begin(), gathered_permuted_graph_edges.end());

    EXPECT_EQ(gathered_graph_edges, gathered_permuted_graph_edges);
}

TEST(GraphPermutation, check_applied_permutation_rgg2d_edgelist) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 64'000;

    kagen::KaGen generator(comm);
    generator.UseEdgeListRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateRGG2D_NM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateRGG2D_NM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}

TEST(GraphPermutation, check_applied_permutation_rgg2d_csr) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 64'000;

    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateRGG2D_NM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateRGG2D_NM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}

TEST(GraphPermutation, check_applied_permutation_rgg2d_edgelist_sparse) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 4'000;

    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateRGG2D_NM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateRGG2D_NM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}

TEST(GraphPermutation, check_applied_permutation_rgg2d_csr_sparse) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 4'000;

    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateRGG2D_NM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateRGG2D_NM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}

TEST(GraphPermutation, check_applied_permutation_gnm_edgelist_sparse) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 4'000;

    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateUndirectedGNM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateUndirectedGNM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}

TEST(GraphPermutation, check_applied_permutation_gnm_csr_sparse) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 4'000;

    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateUndirectedGNM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateUndirectedGNM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}
