#include <gtest/gtest.h>

#include <cmath>

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "tests/util/utils.h"
#include "generators/barabassi/barabassi.h"
#include "edgeweight/edgeweight_utils.h"

using namespace kagen;

TEST(BarabassiEdgeWeightTest, barabassi_constant) {
    PGeneratorConfig config;
    config.n = 32;
    config.m = 32;
    config.generating_edge_weights = true;
    config.edge_weights.edge_weight_type = EdgeWeightType::CONSTANT;

    BarabassiFactory factory;
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

        for(int i = 0; i < gathered_edges.edges.size(); i++) {
            std::cout << "Edge " << i << " (" << gathered_edges.edges[i].first << ", " << gathered_edges.edges[i].second << ") " << gathered_edges.edge_weights[i] << std::endl;
        }
        ASSERT_TRUE(check);
    }
}

TEST(BarabassiEdgeWeightTest, barabassi_hashed) {
    PGeneratorConfig config;
    config.n = 32;
    config.m = 32;
    config.generating_edge_weights = true;
    config.coordinates = true;
    config.edge_weights.edge_weight_type = EdgeWeightType::HASHED;

    BarabassiFactory factory;
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

        // ToDo: Check if edges are the same for both directions
    }
}

TEST(BarabassiEdgeWeightTest, barabassi_random) {
    PGeneratorConfig config;
    config.n = 32;
    config.m = 32;
    config.generating_edge_weights = true;
    config.coordinates = true;
    config.edge_weights.edge_weight_type = EdgeWeightType::RANDOM;

    BarabassiFactory factory;
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
        // ToDo: Check if edges are the same for both directions
    }
}