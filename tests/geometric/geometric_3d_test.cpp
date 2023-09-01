#include <gtest/gtest.h>

#include <cmath>

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/geometric/geometric_2d.h"
#include "kagen/generators/geometric/rgg.h"
#include "tests/util/utils.h"
#include "geometric/rgg_utils.h"

using namespace kagen;

TEST(Geometric3DTest, generates_graph_on_np_PE_n32_r125) {
    PGeneratorConfig config;
    config.n           = 32;
    config.r           = 0.125;
    config.coordinates = true;

    RGG3DFactory factory;
    PEID         size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    config         = factory.NormalizeParameters(config, rank, size, false);
    auto generator = factory.Create(config, rank, size);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    generator->Finalize(MPI_COMM_WORLD);
    auto result = generator->Take();

    Graph complete_graph = kagen::testing::GatherEdgeLists(result);
    auto  global_graph   = kagen::testing::GatherCoordinates3D(result);

    if (rank == 0) {
        EXPECT_EQ(config.n, global_graph.coordinates.second.size());

        // Creating the correct edge list as a test instance
        std::vector<std::pair<SInt, SInt>> edge_List =
            kagen::testing::CreateExpectedRGG3DEdges(config, global_graph);

        // Sorting both lists before comparing them
        std::sort(complete_graph.edges.begin(), complete_graph.edges.end());
        std::sort(edge_List.begin(), edge_List.end());

        EXPECT_EQ(complete_graph.edges, edge_List);
    }
}

TEST(Geometric3DTest, generates_graph_on_np_PE_n16_r10) {
    PGeneratorConfig config;
    config.n           = 16;
    config.r           = 0.1;
    config.coordinates = true;

    RGG3DFactory factory;
    PEID         size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    config         = factory.NormalizeParameters(config, rank, size, false);
    auto generator = factory.Create(config, rank, size);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    generator->Finalize(MPI_COMM_WORLD);
    auto result = generator->Take();

    Graph complete_graph = kagen::testing::GatherEdgeLists(result);
    auto  global_graph   = kagen::testing::GatherCoordinates3D(result);

    if (rank == 0) {
        EXPECT_EQ(config.n, global_graph.coordinates.second.size());

        // Creating the correct edge list as a test instance
        std::vector<std::pair<SInt, SInt>> edge_List =
            kagen::testing::CreateExpectedRGG3DEdges(config, global_graph);

        // Sorting both lists before comparing them
        std::sort(complete_graph.edges.begin(), complete_graph.edges.end());
        std::sort(edge_List.begin(), edge_List.end());

        EXPECT_EQ(complete_graph.edges, edge_List);
    }
}
