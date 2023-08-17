#include <gtest/gtest.h>

#include <cmath>

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/geometric/geometric_2d.h"
#include "kagen/generators/geometric/rgg.h"
#include "tests/util/utils.h"

using namespace kagen;

TEST(EdgeWeightTest, edgeWeight) {
    PGeneratorConfig config;
    config.n = 32;
    config.r = 0.125;
    config.coordinates = true;
    config.edge_weights = true;
    config.edge_weight_type = EdgeWeightType::RANDOM;

    RGG2DFactory factory;
    PEID size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    config = factory.NormalizeParameters(config, rank, size, false);
    auto generator = factory.Create(config, rank, size);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    generator->Finalize(MPI_COMM_WORLD);
    generator->GenerateEdgeWeights(config.edge_weight_type, MPI_COMM_WORLD);
    auto result = generator->Take();

    Graph gathered_edges = kagen::testing::GatherEdgeLists(result);
    Graph gathered_weights;
    kagen::testing::GatherWeights(result, gathered_weights);

    for (SInt i = 0; i < result.edge_weights.size(); i++) {
        std::cout << "Result Edge " << i << " rank " << rank << " (" << result.edges.at(i).first << ", "
                  << result.edges.at(i).second << ") " << result.edge_weights.at(i) << std::endl;
    }

    if (rank == 0) {

        for (SInt i = 0; i < gathered_weights.edge_weights.size(); i++) {
            std::cout << "Complete Edge " << i << " (" << gathered_edges.edges.at(i).first << ", "
                      << gathered_edges.edges.at(i).second << ") " << gathered_weights.edge_weights.at(i) << std::endl;
        }
    }
}

//TEST(EdgeWeightTest, edgeWeight3D) {
//    PGeneratorConfig config;
//    config.n = 32;
//    config.r = 0.125;
//    config.coordinates = true;
//    config.edge_weights = true;
//    config.edge_weight_type = EdgeWeightType::HASHED;
//
//    RGG3DFactory factory;
//    PEID size, rank;
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    config = factory.NormalizeParameters(config, rank, size, false);
//    auto generator = factory.Create(config, rank, size);
//    generator->Generate(GraphRepresentation::EDGE_LIST);
//    generator->Finalize(MPI_COMM_WORLD);
//    generator->GenerateEdgeWeights(config.edge_weight_type, MPI_COMM_WORLD);
//    auto result = generator->Take();
//
//    std::cout << "Graph complete " << std::endl;
//
//    Graph gathered_edges = kagen::testing::GatherEdgeLists(result);
//
//    std::cout << "Graph complete " << std::endl;
//
//    Graph gathered_weights;
//    kagen::testing::GatherWeights(result, gathered_weights);
//
//    std::cout << "Graph complete " << std::endl;
//
//    if (rank == 0) {
//        for (SInt i = 0; i < gathered_weights.edge_weights.size(); i++) {
//            std::cout << "Complete Edge " << i << " (" << gathered_edges.edges.at(i).first << ", "
//                      << gathered_edges.edges.at(i).second << ") " << gathered_weights.edge_weights.at(i) << std::endl;
//        }
//    }
//}