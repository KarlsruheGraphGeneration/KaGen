#include <gtest/gtest.h>

#include <cmath>

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "tests/util/utils.h"
#include "edgeweight/edgeweight_utils.h"
#include "generators/geometric/rgg.h"

using namespace kagen;

TEST(GeometricEdgeWeightTest, geometric_constant) {
    PGeneratorConfig config;
    config.n = 32;
    config.m = 32;
    config.generating_edge_weights = true;
    config.edge_weights.edge_weight_type = EdgeWeightType::CONSTANT;

    RGG2DFactory factory;
    PEID size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    config = factory.NormalizeParameters(config, rank, size, false);
    auto generator = factory.Create(config, rank, size);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    generator->Finalize(MPI_COMM_WORLD);
    generator->GenerateEdgeWeights(config.edge_weights, MPI_COMM_WORLD);
    auto result = generator->Take();

    Graph gathered_edges = kagen::testing::GatherEdgeLists(result);

    if (rank == 0) {
        ASSERT_EQ(gathered_edges.edges.size(), gathered_edges.edge_weights.size());

        bool check = true;
        for (SInt i = 0; i < gathered_edges.edge_weights.size(); i++) {
            check = check && (gathered_edges.edge_weights.at(i) == 1);
        }
        ASSERT_TRUE(check);
    }
}

TEST(GeometricEdgeWeightTest, geometric_random) {
    PGeneratorConfig config;
    config.n = 32;
    config.m = 32;
    config.coordinates = true;
    config.generating_edge_weights = true;
    config.edge_weights.edge_weight_type = EdgeWeightType::RANDOM;

    RGG2DFactory factory;
    PEID size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    config = factory.NormalizeParameters(config, rank, size, false);
    auto generator = factory.Create(config, rank, size);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    generator->Finalize(MPI_COMM_WORLD);
    generator->GenerateEdgeWeights(config.edge_weights, MPI_COMM_WORLD);
    auto result = generator->Take();

    Graph gathered_edges = kagen::testing::GatherEdgeLists(result);

    if (rank == 0) {
        ASSERT_EQ(gathered_edges.edges.size(), gathered_edges.edge_weights.size());

        kagen::testing::CheckReverseEdges(gathered_edges.edges, gathered_edges.edge_weights);
    }
}

TEST(GeometricEdgeWeightTest, geometric_distance) {
    PGeneratorConfig config;
    config.n = 32;
    config.m = 32;
    config.coordinates = true;
    config.generating_edge_weights = true;
    config.edge_weights.edge_weight_type = EdgeWeightType::DISTANCE;

    RGG2DFactory factory;
    PEID size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    config = factory.NormalizeParameters(config, rank, size, false);
    auto generator = factory.Create(config, rank, size);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    generator->Finalize(MPI_COMM_WORLD);
    generator->GenerateEdgeWeights(config.edge_weights, MPI_COMM_WORLD);
    auto result = generator->Take();

    Graph gathered_edges = kagen::testing::GatherEdgeLists(result);

    if (rank == 0) {
        ASSERT_EQ(gathered_edges.edges.size(), gathered_edges.edge_weights.size());

//        for(int i = 0; i < gathered_edges.edges.size(); i++) {
//            std::cout << "Edge " << i << " (" << gathered_edges.edges[i].first << ", " << gathered_edges.edges[i].second << ") " << gathered_edges.edge_weights[i] << std::endl;
//        }
    }
}

TEST(GeometricEdgeWeightTest, geometric_hashed) {
    PGeneratorConfig config;
    config.n = 32;
    config.m = 32;
    config.generating_edge_weights = true;
    config.edge_weights.edge_weight_type = EdgeWeightType::HASHED;

    RGG2DFactory factory;
    PEID size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    config = factory.NormalizeParameters(config, rank, size, false);
    auto generator = factory.Create(config, rank, size);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    generator->Finalize(MPI_COMM_WORLD);
    generator->GenerateEdgeWeights(config.edge_weights, MPI_COMM_WORLD);
    auto result = generator->Take();

    Graph gathered_edges = kagen::testing::GatherEdgeLists(result);

    if (rank == 0) {
        ASSERT_EQ(gathered_edges.edges.size(), gathered_edges.edge_weights.size());

        for(int i = 0; i < gathered_edges.edges.size(); i++) {
            std::cout << "Edge " << i << " (" << gathered_edges.edges[i].first << ", " << gathered_edges.edges[i].second << ") " << gathered_edges.edge_weights[i] << std::endl;
        }

        kagen::testing::CheckReverseEdges(gathered_edges.edges, gathered_edges.edge_weights);
    }
}