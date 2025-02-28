#include "kagen/context.h"
#include "kagen/generators/geometric/rgg.h"

#include <gtest/gtest.h>

#include "tests/gather.h"
#include "tests/geometric/utils.h"
#include "tests/utils.h"

using namespace kagen;

namespace kagen::testing {

void ValidateRGG3D(const Graph& graph, const double radius) {
    Graph global_graph = GatherEdgeLists(graph);
    if (IsRoot(MPI_COMM_WORLD)) {
        const auto& actual_edges   = global_graph.edges;
        const auto  expected_edges = GenerateLocalRGG3DEdges(global_graph, radius);
        EXPECT_EQ(actual_edges, expected_edges);
    }
}

Graph GenerateAndValidateRGG3D(const SInt n, const double radius) {
    PEID size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    PGeneratorConfig config = {
        .n           = n,
        .r           = radius,
        .coordinates = true,
    };

    RGG3DFactory factory;
    config = factory.NormalizeParameters(config, rank, size, false);

    auto generator = factory.Create(config, rank, size);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    generator->Finalize(MPI_COMM_WORLD);

    const auto graph = generator->Take();

    ValidateRGG3D(graph, radius);

    return graph;
}

TEST(RGG3D, N1) {
    for (const double radius: {0.01, 0.5, 1.0}) {
        Graph global_graph = GenerateAndValidateRGG3D(1, radius);
        EXPECT_LE(global_graph.NumberOfGlobalVertices(), 1);
        EXPECT_EQ(global_graph.NumberOfGlobalEdges(), 0);
    }
}

TEST(RGG3D, N16) {
    for (const double radius: {0.001, 0.01, 0.1, 0.5, 1.0, 2.0}) {
        Graph global_graph = GenerateAndValidateRGG3D(16, radius);
        EXPECT_EQ(global_graph.NumberOfGlobalVertices(), 16);
    }
}

TEST(RGG3D, N2000) {
    Graph global_graph = GenerateAndValidateRGG3D(2000, 0.01);
    EXPECT_EQ(global_graph.NumberOfGlobalVertices(), 2000);
}

} // namespace kagen::testing
