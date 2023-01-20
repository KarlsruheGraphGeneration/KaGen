#include <gtest/gtest.h>
#include <mpi.h>

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/facade.h"

using namespace kagen;

PGeneratorConfig CreateConfig() {
    PGeneratorConfig config;
    config.n      = 1 << 10; // 2^10 vertices
    config.m      = 1 << 12; // 2^12 edges
    config.plexp  = 3.0;
    config.p      = 0.5;
    config.grid_x = 5;
    config.grid_y = 5;
    config.grid_z = 5;
    config.rmat_a = 0.1;
    config.rmat_b = 0.1;
    config.rmat_c = 0.1;

    // Fail test if the graph is not simple
    config.validate_simple_graph = true;
    return config;
}

void TestGenerator(GeneratorType type) {
    auto config      = CreateConfig();
    config.generator = type;
    Generate(config, GraphRepresentation::EDGE_LIST, MPI_COMM_WORLD);
}

TEST(SimpleGraphsTest, gnm) { // @todo
    // TestGenerator(GeneratorType::GNM_UNDIRECTED);
}

TEST(SimpleGraphsTest, gnp) {
    TestGenerator(GeneratorType::GNP_UNDIRECTED);
}

TEST(SimpleGraphsTest, rgg2d) {
    TestGenerator(GeneratorType::RGG_2D);
}

TEST(SimpleGraphsTest, rgg3d) {
    TestGenerator(GeneratorType::RGG_3D);
}

TEST(SimpleGraphsTest, rhg) {
    TestGenerator(GeneratorType::RHG);
}

TEST(SimpleGraphsTest, grid2d) { // @todo
    // TestGenerator(GeneratorType::GRID_2D);
}

TEST(SimpleGraphsTest, grid3d) { // @todo
    // TestGenerator(GeneratorType::GRID_3D);
}

TEST(SimpleGraphsTest, ba) {
    TestGenerator(GeneratorType::BA);
}

TEST(SimpleGraphsTest, kronecker) {
    TestGenerator(GeneratorType::KRONECKER);
}

TEST(SimpleGraphsTest, rmat) {
    TestGenerator(GeneratorType::RMAT);
}
