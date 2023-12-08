#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/geometric/geometric_2d.h"
#include "kagen/generators/geometric/rgg.h"

#include <gtest/gtest.h>

#include "tests/gather.h"
#include "tests/geometric/utils.h"

using namespace kagen;

namespace {
void validate_graph(const Graph& local_graph, const PGeneratorConfig& config) {
    Graph global_graph = kagen::testing::GatherEdgeLists(local_graph);

    PEID rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        EXPECT_EQ(config.n, global_graph.coordinates.second.size());

        // Creating the correct edge list as a test instance
        std::vector<std::pair<SInt, SInt>> expected_edges =
            kagen::testing::CreateExpectedRGG3DEdges(config, global_graph);

        // Sorting both lists before comparing them
        std::sort(global_graph.edges.begin(), global_graph.edges.end());
        std::sort(expected_edges.begin(), expected_edges.end());

        EXPECT_EQ(global_graph.edges, expected_edges);
    }
}

void test_configuration(const SInt n, const double radius) {
    PGeneratorConfig config;
    config.n           = n;
    config.r           = radius;
    config.coordinates = true;

    RGG3DFactory factory;
    PEID         size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    config         = factory.NormalizeParameters(config, rank, size, false);
    auto generator = factory.Create(config, rank, size);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    generator->Finalize(MPI_COMM_WORLD);
    const auto local_graph = generator->Take();

    validate_graph(local_graph, config);
}
} // namespace

TEST(Geometric3DTest, generates_graph_on_np_PE_n32_r125) {
    test_configuration(32, 0.125);
}

TEST(Geometric3DTest, generates_graph_on_np_PE_n16_r10) {
    test_configuration(16, 0.1);
}

TEST(Geometric3DTest, generates_graph_on_np_PE_n512_r01) {
    test_configuration(512, 0.01);
}
