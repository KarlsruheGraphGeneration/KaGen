#include <gtest/gtest.h>

#include <cmath>

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/geometric/geometric_2d.h"
#include "kagen/generators/geometric/rgg.h"
#include "tests/util/utils.h"
#include "geometric/rgg_utils.h"

using namespace kagen;

TEST(Geometric1PETest, generates_2D_graph) {
    PGeneratorConfig config;
    config.n           = 32;
    config.r           = 0.5;
    config.coordinates = true;

    RGG2DFactory factory;
    config         = factory.NormalizeParameters(config, 0, 1, false);
    auto generator = factory.Create(config, 0, 1);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    generator->Finalize(MPI_COMM_WORLD);
    auto result = generator->Take();

    ASSERT_EQ(config.n, result.coordinates.first.size());

    // Creating the correct edge list as a test instance
    std::vector<std::pair<SInt, SInt>> edge_List = kagen::testing::CreateExpectedRGG2DEdges(config, result);

    // Sorting both lists before comparing them
    std::sort(result.edges.begin(), result.edges.end());
    std::sort(edge_List.begin(), edge_List.end());

    ASSERT_EQ(result.edges, edge_List);
}

TEST(Geometric1PETest, generates_3D_graph) {
    PGeneratorConfig config;
    config.n           = 32;
    config.r           = 0.5;
    config.coordinates = true;

    RGG3DFactory factory;
    config         = factory.NormalizeParameters(config, 0, 1, false);
    auto generator = factory.Create(config, 0, 1);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    generator->Finalize(MPI_COMM_WORLD);
    auto result = generator->Take();

    ASSERT_EQ(config.n, result.coordinates.second.size());

    // Creating the correct edge list as a test instance
    std::vector<std::pair<SInt, SInt>> edge_List = kagen::testing::CreateExpectedRGG3DEdges(config, result);

    // Sorting both lists before comparing them
    std::sort(result.edges.begin(), result.edges.end());
    std::sort(edge_List.begin(), edge_List.end());

    ASSERT_EQ(result.edges, edge_List);
}