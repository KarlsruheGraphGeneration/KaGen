#include "kagen/edge_range.h"

#include "kagen/kagen.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <functional>
#include <numeric>
#include <string>

#include "tests/gather.h"
#include "tests/utils.h"
#include "tools/converter.h"
#include "tools/geometry.h"

using namespace kagen;

using GeneratorFunc = std::function<Graph(KaGen&, SInt, SInt)>;

MATCHER(EdgeIndexMatches, "") {
    auto [iter, expected_idx] = arg;
    return iter.edge_index() == expected_idx;
}

struct EdgeRangeTestFixture : public ::testing::TestWithParam<std::tuple<std::string, GeneratorFunc>> {};

INSTANTIATE_TEST_SUITE_P(
    EdgeRangeTests, EdgeRangeTestFixture,
    ::testing::Values(
        std::make_tuple("GNM", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateUndirectedGNM(n, m); })),
        std::make_tuple("RMAT", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRMAT(n, m, 0.56, 0.19, 0.19); })),
        std::make_tuple("RGG2D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRGG2D_NM(n, m); })),
        std::make_tuple("RGG3D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRGG3D_NM(n, m); })),
        std::make_tuple("RHG", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRHG_NM(2.6, n, m); })),
        std::make_tuple("Grid2D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateGrid2D_NM(n, m); })),
        std::make_tuple("Grid3D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateGrid3D_NM(n, m); }))),
    [](const ::testing::TestParamInfo<EdgeRangeTestFixture::ParamType>& info) {
        return std::get<0>(info.param);
    });

TEST_P(EdgeRangeTestFixture, iterate_edgelist_representation) {
    using ::testing::ElementsAreArray;
    using ::testing::Pointwise;

    auto [name, generate] = GetParam();
    const SInt n = 1000;
    const SInt m = 16 * n;

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    Graph graph = generate(generator, n, m);

    Edgelist expected = graph.edges;
    EdgeRange edge_range(graph);

    // Check edges match and indices are consecutive
    EXPECT_THAT(std::vector(edge_range.begin(), edge_range.end()), ElementsAreArray(expected));
    
    std::vector<EdgeRange::iterator> iterators;
    for (auto it = edge_range.begin(); it != edge_range.end(); ++it) {
        iterators.push_back(it);
    }
    std::vector<std::size_t> expected_indices(edge_range.size());
    std::iota(expected_indices.begin(), expected_indices.end(), 0);
    EXPECT_THAT(iterators, Pointwise(EdgeIndexMatches(), expected_indices));
}

TEST_P(EdgeRangeTestFixture, iterate_sparse_edgelist_representation) {
    using ::testing::ElementsAreArray;
    using ::testing::Pointwise;

    auto [name, generate] = GetParam();
    const SInt n = 1000;
    const SInt m = 2 * n;

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    Graph graph = generate(generator, n, m);

    Edgelist expected = graph.edges;
    EdgeRange edge_range(graph);

    // Check edges match and indices are consecutive
    EXPECT_THAT(std::vector(edge_range.begin(), edge_range.end()), ElementsAreArray(expected));
    
    std::vector<EdgeRange::iterator> iterators;
    for (auto it = edge_range.begin(); it != edge_range.end(); ++it) {
        iterators.push_back(it);
    }
    std::vector<std::size_t> expected_indices(edge_range.size());
    std::iota(expected_indices.begin(), expected_indices.end(), 0);
    EXPECT_THAT(iterators, Pointwise(EdgeIndexMatches(), expected_indices));
}

TEST_P(EdgeRangeTestFixture, iterate_csr_representation) {
    using ::testing::ElementsAreArray;
    using ::testing::Pointwise;

    auto [name, generate] = GetParam();
    const SInt n = 1000;
    const SInt m = 16 * n;

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseCSRRepresentation();
    Graph graph = generate(generator, n, m);

    Edgelist expected = BuildEdgeListFromCSR(graph.vertex_range, graph.xadj, graph.adjncy);
    EdgeRange edge_range(graph);

    // Check edges match and indices are consecutive
    EXPECT_THAT(std::vector(edge_range.begin(), edge_range.end()), ElementsAreArray(expected));
    
    std::vector<EdgeRange::iterator> iterators;
    for (auto it = edge_range.begin(); it != edge_range.end(); ++it) {
        iterators.push_back(it);
    }
    std::vector<std::size_t> expected_indices(edge_range.size());
    std::iota(expected_indices.begin(), expected_indices.end(), 0);
    EXPECT_THAT(iterators, Pointwise(EdgeIndexMatches(), expected_indices));
}

TEST_P(EdgeRangeTestFixture, iterate_sparse_csr_representation) {
    using ::testing::ElementsAreArray;
    using ::testing::Pointwise;

    auto [name, generate] = GetParam();
    const SInt n = 1000;
    const SInt m = 2 * n;

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseCSRRepresentation();
    Graph graph = generate(generator, n, m);

    Edgelist expected = BuildEdgeListFromCSR(graph.vertex_range, graph.xadj, graph.adjncy);
    EdgeRange edge_range(graph);

    // Check edges match and indices are consecutive
    EXPECT_THAT(std::vector(edge_range.begin(), edge_range.end()), ElementsAreArray(expected));
    
    std::vector<EdgeRange::iterator> iterators;
    for (auto it = edge_range.begin(); it != edge_range.end(); ++it) {
        iterators.push_back(it);
    }
    std::vector<std::size_t> expected_indices(edge_range.size());
    std::iota(expected_indices.begin(), expected_indices.end(), 0);
    EXPECT_THAT(iterators, Pointwise(EdgeIndexMatches(), expected_indices));
}
