#include <gtest/gtest.h>

#include <cmath>

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/geometric/geometric_2d.h"
#include "kagen/generators/geometric/rgg.h"
#include "tests/util/utils.h"

using namespace kagen;

TEST(GEOMETRIC3D, generates_graph_on_one_PE) {
    PGeneratorConfig config;
    config.n           = 32;
    config.r           = 0.5;
    config.coordinates = true;

    RGG3DFactory factory;
    config         = factory.NormalizeParameters(config, 0, 1, false);
    auto generator = factory.Create(config, 0, 1);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    auto result = generator->Take();

    // Creating the correct edge list as a test instance
    std::vector<std::pair<SInt, SInt>> edgeList;
    for (SInt i = 0; i < config.n; i++) {
        auto [x1, y1, z1] = result.coordinates.second[i];
        for (SInt j = 0; j < config.n; j++) {
            auto [x2, y2, z2] = result.coordinates.second[j];
            // Comparing all coordinates
            if (i != j && std::hypot(x1 - x2, y1 - y2, z1 - z2) <= config.r) {
                edgeList.emplace_back(i, j);
            }
        }
    }
    // Sorting both lists before comparing them
    std::sort(result.edges.begin(), result.edges.end());
    std::sort(edgeList.begin(), edgeList.end());

    ASSERT_EQ(result.edges, edgeList);
}

TEST(GEOMETRIC3D, generates_graph_on_np_PE_n32_r125) {
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
    auto result = generator->Take();

    Graph complete_graph = kagen::testing::GatherEdgeLists(result);
    auto global_graph = kagen::testing::GatherCoordinates3D(result);

    if (rank == 0) {
        // Creating the correct edge list as a test instance
        std::vector<std::pair<SInt, SInt>> edgeList;
        for (SInt i = 0; i < config.n; i++) {
            auto [x1, y1, z1] = global_graph.coordinates.second[i];
            for (SInt j = 0; j < config.n; j++) {
                auto [x2, y2, z2] = global_graph.coordinates.second[j];

                // Comparing all coordinates
                if (i != j && std::hypot(x1 - x2, y1 - y2, z1 - z2) < config.r) {
                    edgeList.emplace_back(i, j);
                }
            }
        }
        // Sorting both lists before comparing them
        std::sort(complete_graph.edges.begin(), complete_graph.edges.end());
        std::sort(edgeList.begin(), edgeList.end());

        ASSERT_EQ(complete_graph.edges, edgeList);
    }
}

TEST(GEOMETRIC3D, generates_graph_on_np_PE_n16_r10) {
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
    auto result = generator->Take();

    Graph complete_graph = kagen::testing::GatherEdgeLists(result);
    auto global_graph = kagen::testing::GatherCoordinates3D(result);

    if (rank == 0) {
        // Creating the correct edge list as a test instance
        std::vector<std::pair<SInt, SInt>> edgeList;
        for (SInt i = 0; i < config.n; i++) {
            auto [x1, y1, z1] = global_graph.coordinates.second[i];
            for (SInt j = 0; j < config.n; j++) {
                auto [x2, y2, z2] = global_graph.coordinates.second[j];

                // Comparing all coordinates
                if (i != j && std::hypot(x1 - x2, y1 - y2, z1 - z2) < config.r) {
                    edgeList.emplace_back(i, j);
                }
            }
        }
        // Sorting both lists before comparing them
        std::sort(complete_graph.edges.begin(), complete_graph.edges.end());
        std::sort(edgeList.begin(), edgeList.end());

        ASSERT_EQ(complete_graph.edges, edgeList);
    }
}