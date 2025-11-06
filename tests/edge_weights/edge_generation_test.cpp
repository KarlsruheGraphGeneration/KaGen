#include "kagen/kagen.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "kagen/edge_range.h"
#include "tests/gather.h"
#include "tests/utils.h"
#include "tools/geometry.h"
#include <functional>
#include <map>

using namespace kagen;
using WeightRange   = std::pair<SInt, SInt>;
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
    EdgeRange  edge_range(graph);
    EXPECT_EQ(edge_range.size(), graph.edge_weights.size());
    if (graph.NumberOfLocalEdges() == 0) {
        return;
    }

    auto check_or_insert_edge = [](auto& map, auto edge, SSInt weight) {
        auto& [u, v] = edge;
        if (u > v) {
            std::swap(u, v);
        }
        auto it = map.find(edge);
        if (it != map.end()) {
            EXPECT_EQ(it->second, weight);
        } else {
            map.emplace(edge, weight);
        }
    };

    std::map<EdgeRange::Edge, SSInt> edge_to_weight_map;
    for (auto it = edge_range.begin(); it != edge_range.end(); ++it) {
        check_or_insert_edge(edge_to_weight_map, *it, graph.edge_weights[it.edge_index()]);
    }
    EXPECT_EQ(edge_to_weight_map.size() * 2, edge_range.size());
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
struct UniformWeightTestFixture : public ::testing::TestWithParam<std::tuple<std::string, GeneratorFunc>> {};

INSTANTIATE_TEST_SUITE_P(
    UniformWeightTests, UniformWeightTestFixture,
    ::testing::Values(
        std::make_tuple("GNM", GeneratorFunc([](KaGen& gen, SInt n, SInt m) {
                            return gen.GenerateUndirectedGNM(n, m);
                        })),
        std::make_tuple("RMAT", GeneratorFunc([](KaGen& gen, SInt n, SInt m) {
                            return gen.GenerateRMAT(n, m, 0.56, 0.19, 0.19);
                        })),
        std::make_tuple("RGG2D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRGG2D_NM(n, m); })),
        std::make_tuple("RGG3D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRGG3D_NM(n, m); })),
        std::make_tuple("RHG", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRHG_NM(2.6, n, m); })),
        std::make_tuple("Grid2D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) {
                            return gen.GenerateGrid2D_NM(n, m);
                        })),
        std::make_tuple("Grid3D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) {
                            return gen.GenerateGrid3D_NM(n, m);
                        }))),
    [](const ::testing::TestParamInfo<UniformWeightTestFixture::ParamType>& info) { return std::get<0>(info.param); });

TEST_P(UniformWeightTestFixture, weights_in_range_edgelist_representation) {
    std::string       name     = std::get<0>(GetParam());
    GeneratorFunc     generate = std::get<1>(GetParam());
    const SInt        n        = 1000;
    const SInt        m        = 16 * n;
    const WeightRange weight_range{1, 100};

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    generator.ConfigureEdgeWeightGeneration(
        kagen::EdgeWeightGeneratorType::UNIFORM_RANDOM, weight_range.first, weight_range.second);

    Graph graph = generate(generator, n, m);
    check_weights_range(graph, weight_range);
}

TEST_P(UniformWeightTestFixture, weights_in_range_csr_representation) {
    std::string       name     = std::get<0>(GetParam());
    GeneratorFunc     generate = std::get<1>(GetParam());
    const SInt        n        = 1000;
    const SInt        m        = 16 * n;
    const WeightRange weight_range{1, 100};

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(
        kagen::EdgeWeightGeneratorType::UNIFORM_RANDOM, weight_range.first, weight_range.second);

    Graph graph = generate(generator, n, m);
    check_weights_range(graph, weight_range);
}

TEST_P(UniformWeightTestFixture, correct_backedge_weights_edgelist_representation) {
    std::string       name     = std::get<0>(GetParam());
    GeneratorFunc     generate = std::get<1>(GetParam());
    const SInt        n        = 1000;
    const SInt        m        = 16 * n;
    const WeightRange weight_range{1, 100};

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    generator.ConfigureEdgeWeightGeneration(
        kagen::EdgeWeightGeneratorType::UNIFORM_RANDOM, weight_range.first, weight_range.second);

    Graph graph = generate(generator, n, m);
    check_backedge_weight(graph);
}

TEST_P(UniformWeightTestFixture, correct_backedge_weights_csr_representation) {
    std::string       name     = std::get<0>(GetParam());
    GeneratorFunc     generate = std::get<1>(GetParam());
    const SInt        n        = 1000;
    const SInt        m        = 16 * n;
    const WeightRange weight_range{1, 100};

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(
        kagen::EdgeWeightGeneratorType::UNIFORM_RANDOM, weight_range.first, weight_range.second);

    Graph graph = generate(generator, n, m);
    check_backedge_weight(graph);
}

// Test fixture for Euclidean weight tests
struct EuclideanWeightTestFixture
    : public ::testing::TestWithParam<
          std::tuple<std::string, std::function<Graph(KaGen&, SInt, double)>, bool>> {};

INSTANTIATE_TEST_SUITE_P(
    EuclideanWeightTests, EuclideanWeightTestFixture,
    ::testing::Values(
        std::make_tuple("RGG2D", std::function<Graph(KaGen&, SInt, double)>([](KaGen& gen, SInt n, double r) {
                            return gen.GenerateRGG2D(n, r, true);
                        }), true),
        std::make_tuple("RGG3D", std::function<Graph(KaGen&, SInt, double)>([](KaGen& gen, SInt n, double r) {
                            return gen.GenerateRGG3D(n, r, true);
                        }), false)),
    [](const ::testing::TestParamInfo<EuclideanWeightTestFixture::ParamType>& info) {
        return std::get<0>(info.param) + "is_2d_" + std::to_string(std::get<2>(info.param));
    });

TEST_P(EuclideanWeightTestFixture, euclidean_weights_edgelist_representation) {
    std::string       name     = std::get<0>(GetParam());
    auto              generate = std::get<1>(GetParam());
    bool              is_2d    = std::get<2>(GetParam());
    const SInt        n        = 1000;
    const double      radius   = 0.05;
    const WeightRange weight_range{1, 1'000'000'000};

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    generator.ConfigureEdgeWeightGeneration(
        kagen::EdgeWeightGeneratorType::EUCLIDEAN_DISTANCE, weight_range.first, weight_range.second);

    Graph graph = generate(generator, n, radius);
    check_euclidean_weights(graph, radius, weight_range, is_2d);
    check_weights_range(graph, weight_range);
}

TEST_P(EuclideanWeightTestFixture, euclidean_weights_csr_representation) {
    std::string       name     = std::get<0>(GetParam());
    auto              generate = std::get<1>(GetParam());
    bool              is_2d    = std::get<2>(GetParam());
    const SInt        n        = 1000;
    const double      radius   = 0.05;
    const WeightRange weight_range{1, 1'000'000'000};

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(
        kagen::EdgeWeightGeneratorType::EUCLIDEAN_DISTANCE, weight_range.first, weight_range.second);

    Graph graph = generate(generator, n, radius);
    check_euclidean_weights(graph, radius, weight_range, is_2d);
    check_weights_range(graph, weight_range);
}
