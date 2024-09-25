#include "kagen/kagen.h"

#include <gtest/gtest.h>

#include "tests/gather.h"
#include "tests/geometric/utils.h"

using namespace kagen;
using WeightRange = std::pair<SInt, SInt>;

bool is_in_weight_range(SInt weight, const WeightRange& range) {
    return range.first <= weight && weight < range.second;
}

void check_weights_range(const Graph& graph, const WeightRange& range) {
    EXPECT_EQ(graph.NumberOfLocalEdges(), graph.edge_weights.size());
    for (const auto& weight: graph.edge_weights) {
        EXPECT_TRUE(is_in_weight_range(weight, range));
    }
}

void check_weights_range(KaGen& generator, const WeightRange& weight_range) {
    const SInt n = 1000;
    const SInt m = 16 * n;
    // GNM
    {
        kagen::Graph graph = generator.GenerateUndirectedGNM(n, m);
        check_weights_range(graph, weight_range);
    }
    // RMAT
    {
        kagen::Graph graph = generator.GenerateRMAT(n, m, 0.56, 0.19, 0.19);
        check_weights_range(graph, weight_range);
    }
    // RGG2D
    {
        kagen::Graph graph = generator.GenerateRGG2D_NM(n, m);
        check_weights_range(graph, weight_range);
    }
    // RGG3D
    {
        kagen::Graph graph = generator.GenerateRGG3D_NM(n, m);
        check_weights_range(graph, weight_range);
    }
    // RHG
    {
        kagen::Graph graph = generator.GenerateRHG_NM(2.6, n, m);
        check_weights_range(graph, weight_range);
    }
    // GRID2D
    {
        kagen::Graph graph = generator.GenerateGrid2D_NM(n, m);
        check_weights_range(graph, weight_range);
    }
    // GRID2D
    {
        kagen::Graph graph = generator.GenerateGrid3D_NM(n, m);
        check_weights_range(graph, weight_range);
    }
}

TEST(EdgeWeightsTest, edge_weights_in_range_for_edge_list_representation) {
    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    const WeightRange weight_range{1, 100};
    generator.ConfigureEdgeWeightGeneration(
        kagen::EdgeWeightGeneratorType::UNIFORM_RANDOM, weight_range.first, weight_range.second);
    check_weights_range(generator, weight_range);
}

TEST(EdgeWeightsTest, edge_weights_in_range_for_CSR_representation) {
    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseCSRRepresentation();
    const WeightRange weight_range{1, 100};
    generator.ConfigureEdgeWeightGeneration(
        kagen::EdgeWeightGeneratorType::UNIFORM_RANDOM, weight_range.first, weight_range.second);
    check_weights_range(generator, weight_range);
}
