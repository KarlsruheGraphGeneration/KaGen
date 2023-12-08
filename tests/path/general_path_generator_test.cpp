#include "kagen/context.h"
#include "kagen/generators/path/path_directed.h"

#include <gtest/gtest.h>

#include "tests/gather.h"

class PathGeneratorTestFixture : public ::testing::TestWithParam<kagen::SInt> {};
INSTANTIATE_TEST_SUITE_P(PathGenerationTests, PathGeneratorTestFixture, ::testing::Values(10, 500, 16384, 30001));

namespace {
constexpr bool debug_output = false;

auto ComputeStartEndVertexAndSuccesors(const kagen::Graph& graph, kagen::SInt n) {
    using namespace kagen;
    std::vector<SInt> in_degree(n);
    std::vector<SInt> out_degree(n);
    std::vector<SInt> successor_array(n);
    for (const auto& edge: graph.edges) {
        if constexpr (debug_output) {
            std::cout << edge.first << " " << edge.second << std::endl;
        }
        ++out_degree[edge.first];
        ++in_degree[edge.second];
        successor_array[edge.first] = edge.second;
    }
    SInt in_degree_zero_vertex  = -1;
    SInt in_degree_zero         = 0;
    SInt in_degree_one          = 0;
    SInt out_degree_zero_vertex = -1;
    SInt out_degree_zero        = 0;
    SInt out_degree_one         = 0;
    for (SInt i = 0; i < n; ++i) {
        if (in_degree[i] == 0) {
            ++in_degree_zero;
            in_degree_zero_vertex = i;
        }
        if (in_degree[i] == 1) {
            ++in_degree_one;
        }
        if (out_degree[i] == 0) {
            ++out_degree_zero;
            out_degree_zero_vertex = i;
        }
        if (out_degree[i] == 1) {
            ++out_degree_one;
        }
    }
    EXPECT_EQ(in_degree_zero, 1);
    EXPECT_EQ(out_degree_zero, 1);
    EXPECT_EQ(in_degree_one, n - 1);
    EXPECT_EQ(out_degree_one, n - 1);
    return std::make_tuple(in_degree_zero_vertex, out_degree_zero_vertex, successor_array);
}

void WalkPath(const kagen::Graph& graph, kagen::SInt n) {
    using namespace kagen;
    const auto [first_vertex, last_vertex, successor_array] = ComputeStartEndVertexAndSuccesors(graph, n);
    SInt cur_vertex                                         = first_vertex;
    for (SInt i = 0; i + 1 < n; ++i) {
        cur_vertex = successor_array[cur_vertex];
    }
    if constexpr (debug_output) {
        std::cout << "first: " << first_vertex << " n: " << n << std::endl;
    }
    EXPECT_EQ(cur_vertex, last_vertex);
}
} // namespace

TEST_P(PathGeneratorTestFixture, path_generation_without_permutation) {
    using namespace kagen;
    PEID size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    PGeneratorConfig config;
    config.n = GetParam();

    PathDirected generator(config, rank, size);
    generator.Generate(GraphRepresentation::EDGE_LIST);
    generator.Finalize(MPI_COMM_WORLD);
    const Graph graph = kagen::testing::GatherGraph(generator.Take());

    WalkPath(graph, config.n);
}

TEST_P(PathGeneratorTestFixture, path_generation_with_permutation) {
    using namespace kagen;
    PEID size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    PGeneratorConfig config;
    config.n       = GetParam();
    config.permute = true;

    PathDirected generator(config, rank, size);
    generator.Generate(GraphRepresentation::EDGE_LIST);
    generator.Finalize(MPI_COMM_WORLD);
    const Graph graph = kagen::testing::GatherGraph(generator.Take());

    WalkPath(graph, config.n);
}
