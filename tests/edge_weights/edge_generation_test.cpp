#include "kagen/kagen.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <functional>
#include <map>

#include "tests/gather.h"
#include "tests/utils.h"
#include "tools/geometry.h"

using namespace kagen;
using WeightRange = std::pair<SInt, SInt>;
using GeneratorFunc = std::function<Graph(KaGen&, SInt, SInt)>;

// Custom matcher for checking if a weight is in range
MATCHER_P(IsInWeightRange, range, "") {
    return range.first <= arg && arg < range.second;
}

void check_weights_range(const Graph& graph, const WeightRange& range) {
    using ::testing::Each;
    EXPECT_EQ(graph.NumberOfLocalEdges(), graph.edge_weights.size());
    EXPECT_THAT(graph.edge_weights, Each(IsInWeightRange(range)));
}

void check_backedge_weight(const Graph& local_graph) {
    auto const graph = kagen::testing::GatherGraph(local_graph);
    EXPECT_EQ(graph.NumberOfLocalEdges(), graph.edge_weights.size());
    if (graph.NumberOfLocalEdges() == 0) {
        return;
    }

    auto check_or_insert_edge = [](auto& map, SInt u, SInt v, SSInt weight) {
        if (u > v) {
            std::swap(u, v);
        }
        auto it = map.find(std::make_pair(u, v));
        if (it != map.end()) {
            EXPECT_EQ(it->second, weight);
        } else {
            map.emplace(std::make_pair(u, v), weight);
        }
    };
    if (graph.representation == GraphRepresentation::EDGE_LIST) {
        std::map<std::pair<SInt, SInt>, SSInt> edge_to_weight_map;
        for (std::size_t i = 0; i < graph.NumberOfLocalEdges(); ++i) {
            auto [u, v]       = graph.edges[i];
            auto const weight = graph.edge_weights[i];
            check_or_insert_edge(edge_to_weight_map, u, v, weight);
        }
        EXPECT_EQ(edge_to_weight_map.size() * 2, graph.NumberOfLocalEdges());
    } else if (graph.representation == GraphRepresentation::CSR) {
        std::size_t const                      n = graph.xadj.size() - 1;
        std::map<std::pair<SInt, SInt>, SSInt> edge_to_weight_map;
        for (std::size_t u_idx = 0; u_idx < n; ++u_idx) {
            for (std::size_t offset = static_cast<std::size_t>(graph.xadj[u_idx]);
                 offset < static_cast<std::size_t>(graph.xadj[u_idx + 1]); ++offset) {
                SInt  u      = static_cast<SInt>(u_idx);
                SInt  v      = static_cast<SInt>(graph.adjncy[offset]);
                SSInt weight = static_cast<SSInt>(graph.edge_weights[offset]);
                check_or_insert_edge(edge_to_weight_map, u, v, weight);
            }
        }
        EXPECT_EQ(edge_to_weight_map.size() * 2, graph.NumberOfLocalEdges());
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

// Test fixture for uniform random weight tests
struct UniformWeightTestFixture : public ::testing::TestWithParam<std::tuple<std::string, GeneratorFunc, GraphRepresentation>> {};

INSTANTIATE_TEST_SUITE_P(
    UniformWeightTests, UniformWeightTestFixture,
    ::testing::Values(
        std::make_tuple("GNM", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateUndirectedGNM(n, m); }), GraphRepresentation::EDGE_LIST),
        std::make_tuple("RMAT", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRMAT(n, m, 0.56, 0.19, 0.19); }), GraphRepresentation::EDGE_LIST),
        std::make_tuple("RGG2D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRGG2D_NM(n, m); }), GraphRepresentation::EDGE_LIST),
        std::make_tuple("RGG3D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRGG3D_NM(n, m); }), GraphRepresentation::EDGE_LIST),
        std::make_tuple("RHG", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRHG_NM(2.6, n, m); }), GraphRepresentation::EDGE_LIST),
        std::make_tuple("Grid2D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateGrid2D_NM(n, m); }), GraphRepresentation::EDGE_LIST),
        std::make_tuple("Grid3D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateGrid3D_NM(n, m); }), GraphRepresentation::EDGE_LIST),
        std::make_tuple("GNM", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateUndirectedGNM(n, m); }), GraphRepresentation::CSR),
        std::make_tuple("RMAT", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRMAT(n, m, 0.56, 0.19, 0.19); }), GraphRepresentation::CSR),
        std::make_tuple("RGG2D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRGG2D_NM(n, m); }), GraphRepresentation::CSR),
        std::make_tuple("RGG3D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRGG3D_NM(n, m); }), GraphRepresentation::CSR),
        std::make_tuple("RHG", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRHG_NM(2.6, n, m); }), GraphRepresentation::CSR),
        std::make_tuple("Grid2D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateGrid2D_NM(n, m); }), GraphRepresentation::CSR),
        std::make_tuple("Grid3D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateGrid3D_NM(n, m); }), GraphRepresentation::CSR)),
    [](const ::testing::TestParamInfo<UniformWeightTestFixture::ParamType>& info) {
        std::string name = std::get<0>(info.param);
        GraphRepresentation repr = std::get<2>(info.param);
        return name + (repr == GraphRepresentation::EDGE_LIST ? "_EdgeList" : "_CSR");
    });

TEST_P(UniformWeightTestFixture, weights_in_range) {
    std::string name = std::get<0>(GetParam());
    GeneratorFunc generate = std::get<1>(GetParam());
    GraphRepresentation repr = std::get<2>(GetParam());
    const SInt n = 1000;
    const SInt m = 16 * n;
    const WeightRange weight_range{1, 100};

    kagen::KaGen generator(MPI_COMM_WORLD);
    if (repr == GraphRepresentation::EDGE_LIST) {
        generator.UseEdgeListRepresentation();
    } else {
        generator.UseCSRRepresentation();
    }
    generator.ConfigureEdgeWeightGeneration(
        kagen::EdgeWeightGeneratorType::UNIFORM_RANDOM, weight_range.first, weight_range.second);
    
    Graph graph = generate(generator, n, m);
    check_weights_range(graph, weight_range);
}

TEST_P(UniformWeightTestFixture, correct_backedge_weights) {
    std::string name = std::get<0>(GetParam());
    GeneratorFunc generate = std::get<1>(GetParam());
    GraphRepresentation repr = std::get<2>(GetParam());
    const SInt n = 1000;
    const SInt m = 16 * n;
    const WeightRange weight_range{1, 100};

    kagen::KaGen generator(MPI_COMM_WORLD);
    if (repr == GraphRepresentation::EDGE_LIST) {
        generator.UseEdgeListRepresentation();
    } else {
        generator.UseCSRRepresentation();
    }
    generator.ConfigureEdgeWeightGeneration(
        kagen::EdgeWeightGeneratorType::UNIFORM_RANDOM, weight_range.first, weight_range.second);
    
    Graph graph = generate(generator, n, m);
    check_backedge_weight(graph);
}

// Test fixture for Euclidean weight tests  
struct EuclideanWeightTestFixture : public ::testing::TestWithParam<std::tuple<std::string, std::function<Graph(KaGen&, SInt, double)>, GraphRepresentation, bool>> {};

INSTANTIATE_TEST_SUITE_P(
    EuclideanWeightTests, EuclideanWeightTestFixture,
    ::testing::Values(
        std::make_tuple("RGG2D", std::function<Graph(KaGen&, SInt, double)>([](KaGen& gen, SInt n, double r) { return gen.GenerateRGG2D(n, r, true); }), GraphRepresentation::EDGE_LIST, true),
        std::make_tuple("RGG3D", std::function<Graph(KaGen&, SInt, double)>([](KaGen& gen, SInt n, double r) { return gen.GenerateRGG3D(n, r, true); }), GraphRepresentation::EDGE_LIST, false),
        std::make_tuple("RGG2D", std::function<Graph(KaGen&, SInt, double)>([](KaGen& gen, SInt n, double r) { return gen.GenerateRGG2D(n, r, true); }), GraphRepresentation::CSR, true),
        std::make_tuple("RGG3D", std::function<Graph(KaGen&, SInt, double)>([](KaGen& gen, SInt n, double r) { return gen.GenerateRGG3D(n, r, true); }), GraphRepresentation::CSR, false)),
    [](const ::testing::TestParamInfo<EuclideanWeightTestFixture::ParamType>& info) {
        std::string name = std::get<0>(info.param);
        GraphRepresentation repr = std::get<2>(info.param);
        return name + (repr == GraphRepresentation::EDGE_LIST ? "_EdgeList" : "_CSR");
    });

TEST_P(EuclideanWeightTestFixture, euclidean_weights) {
    std::string name = std::get<0>(GetParam());
    auto generate = std::get<1>(GetParam());
    GraphRepresentation repr = std::get<2>(GetParam());
    bool is_2d = std::get<3>(GetParam());
    const SInt n = 1000;
    const double radius = 0.05;
    const WeightRange weight_range{1, 1'000'000'000};

    kagen::KaGen generator(MPI_COMM_WORLD);
    if (repr == GraphRepresentation::EDGE_LIST) {
        generator.UseEdgeListRepresentation();
    } else {
        generator.UseCSRRepresentation();
    }
    generator.ConfigureEdgeWeightGeneration(
        kagen::EdgeWeightGeneratorType::EUCLIDEAN_DISTANCE, weight_range.first, weight_range.second);
    
    Graph graph = generate(generator, n, radius);
    check_euclidean_weights(graph, radius, weight_range, is_2d);
    check_weights_range(graph, weight_range);
}
