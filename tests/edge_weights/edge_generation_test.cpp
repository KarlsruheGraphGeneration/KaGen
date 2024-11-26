#include "kagen/kagen.h"

#include <gtest/gtest.h>

#include "tests/gather.h"
#include "tests/utils.h"
#include "tools/geometry.h"

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

void check_euclidean_weights(const Graph& graph, double max_distance, const WeightRange& range, bool is_2d) {
    auto global_graph       = kagen::testing::GatherGraph(graph);
    auto weighted_edge_list = kagen::testing::ConvertToWeightedEdgelist(global_graph);

    auto get_squard_distance = [&](SInt src, SInt dst) {
        if (is_2d) {
            const auto& coords = global_graph.coordinates.first;
            EXPECT_LT(src, coords.size());
            EXPECT_LT(dst, coords.size());
            return PGGeometry<LPFloat>::SquaredEuclideanDistance(coords[src], coords[dst]);
        } else {
            const auto& coords = global_graph.coordinates.second;
            EXPECT_LT(src, coords.size());
            EXPECT_LT(dst, coords.size());
            const auto& v1 = coords[src];
            const auto& v2 = coords[dst];
            LPFloat     x  = std::get<0>(v1) - std::get<0>(v2);
            LPFloat     y  = std::get<1>(v1) - std::get<1>(v2);
            LPFloat     z  = std::get<2>(v1) - std::get<2>(v2);
            return x * x + y * y + z * z;
        }
    };

    auto compute_euclidean_distance = [&](const SInt src, const SInt dst) {
        const auto squared_dist = get_squard_distance(src, dst) / (max_distance * max_distance);
        const auto expected_dist =
            range.first + static_cast<SSInt>((range.second - range.first) * std::sqrt(squared_dist));
        return expected_dist;
    };

    for (const auto& [src, dst, w]: weighted_edge_list) {
        const auto expected_weight = compute_euclidean_distance(src, dst);
        EXPECT_EQ(w, expected_weight);
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

void check_geometric_weights(KaGen& generator, const WeightRange& weight_range) {
    const SInt n = 1000;
    int        rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // RGG2D
    {
        const double radius = 0.05;
        kagen::Graph graph  = generator.GenerateRGG2D(n, radius, true);
        const bool   is_2d  = true;
        check_euclidean_weights(graph, radius, weight_range, is_2d);
        check_weights_range(graph, weight_range);
    }
    // RGG3D
    {
        const double radius = 0.05;
        kagen::Graph graph  = generator.GenerateRGG3D(n, radius, true);
        const bool   is_2d  = false;
        check_euclidean_weights(graph, radius, weight_range, is_2d);
        check_weights_range(graph, weight_range);
    }
    // RHG
    // TODO not implemented yet
    //{
    //    kagen::Graph graph = generator.GenerateRHG_NM(2.6, n, m);
    //    check_weights_range(graph, weight_range);
    //}
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

TEST(EdgeWeightsTest, euclidean_weights_edge_list_representation) {
    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    const WeightRange weight_range{1, 1'000'000'000};
    generator.ConfigureEdgeWeightGeneration(
        kagen::EdgeWeightGeneratorType::EUCLIDEAN_DISTANCE, weight_range.first, weight_range.second);
    check_geometric_weights(generator, weight_range);
}

TEST(EdgeWeightsTest, euclidean_weights_CSR_representation) {
    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseCSRRepresentation();
    const WeightRange weight_range{1, 1'000'000'000};
    generator.ConfigureEdgeWeightGeneration(
        kagen::EdgeWeightGeneratorType::EUCLIDEAN_DISTANCE, weight_range.first, weight_range.second);
    check_geometric_weights(generator, weight_range);
}
