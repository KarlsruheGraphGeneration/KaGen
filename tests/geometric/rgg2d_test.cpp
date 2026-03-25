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
        EXPECT_EQ(config.n, global_graph.coordinates.first.size());

        // Creating the correct edge list as a test instance
        std::vector<std::pair<SInt, SInt>> expected_edges =
            kagen::testing::CreateExpectedRGG2DEdges(config, global_graph);

        // Sorting both lists before comparing them
        std::sort(global_graph.edges.begin(), global_graph.edges.end());
        std::sort(expected_edges.begin(), expected_edges.end());

        EXPECT_EQ(global_graph.edges, expected_edges);
    }
}

void test_configuration(const SInt n, const double radius, const int seed = 1) {
    PGeneratorConfig config;
    config.n           = n;
    config.r           = radius;
    config.seed        = seed;
    config.coordinates = true;

    RGG2DFactory factory;
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

// All tests use:
// - n >= 8 to avoid MPI buffer aliasing in gather infrastructure
// - r <= 0.125 to stay within k <= 1/r^2 constraint for k=64 (3,5,6 PEs)

// === Small graphs ===

TEST(RGG2DTest, n8_r0125) {
    test_configuration(8, 0.125);
}

TEST(RGG2DTest, n8_r005) {
    test_configuration(8, 0.05);
}

TEST(RGG2DTest, n8_r001) {
    test_configuration(8, 0.01);
}

// === Original tests ===

TEST(RGG2DTest, n16_r01) {
    test_configuration(16, 0.1);
}

TEST(RGG2DTest, n32_r0125) {
    test_configuration(32, 0.125);
}

TEST(RGG2DTest, n512_r001) {
    test_configuration(512, 0.01);
}

// === Moderate graphs with varying radius ===

TEST(RGG2DTest, n16_r0125) {
    test_configuration(16, 0.125);
}

TEST(RGG2DTest, n16_r005) {
    test_configuration(16, 0.05);
}

TEST(RGG2DTest, n16_r002) {
    test_configuration(16, 0.02);
}

TEST(RGG2DTest, n32_r01) {
    test_configuration(32, 0.1);
}

TEST(RGG2DTest, n32_r005) {
    test_configuration(32, 0.05);
}

TEST(RGG2DTest, n32_r002) {
    test_configuration(32, 0.02);
}

TEST(RGG2DTest, n64_r01) {
    test_configuration(64, 0.1);
}

TEST(RGG2DTest, n64_r005) {
    test_configuration(64, 0.05);
}

TEST(RGG2DTest, n64_r002) {
    test_configuration(64, 0.02);
}

TEST(RGG2DTest, n100_r008) {
    test_configuration(100, 0.08);
}

TEST(RGG2DTest, n128_r005) {
    test_configuration(128, 0.05);
}

// === Larger graphs ===

TEST(RGG2DTest, n256_r003) {
    test_configuration(256, 0.03);
}

TEST(RGG2DTest, n256_r005) {
    test_configuration(256, 0.05);
}

// === Edge case: maximum allowed radius (r=0.125 boundary) ===
// With k=64: 1/r^2 = 64, so r=0.125 is the exact boundary

TEST(RGG2DTest, n32_r0125_boundary) {
    test_configuration(32, 0.125);
}

TEST(RGG2DTest, n64_r0125) {
    test_configuration(64, 0.125);
}

TEST(RGG2DTest, n16_r0124) {
    // Just below boundary
    test_configuration(16, 0.124);
}

// === Edge case: very small radius (sparse, few/zero edges) ===

TEST(RGG2DTest, n32_r0005) {
    test_configuration(32, 0.005);
}

TEST(RGG2DTest, n64_r0005) {
    test_configuration(64, 0.005);
}

// === Different seeds for reproducibility ===

TEST(RGG2DTest, n32_r01_seed42) {
    test_configuration(32, 0.1, 42);
}

TEST(RGG2DTest, n32_r01_seed123) {
    test_configuration(32, 0.1, 123);
}

TEST(RGG2DTest, n64_r005_seed7) {
    test_configuration(64, 0.05, 7);
}

TEST(RGG2DTest, n64_r005_seed999) {
    test_configuration(64, 0.05, 999);
}

// === Non-power-of-two node counts ===

TEST(RGG2DTest, n9_r01) {
    test_configuration(9, 0.1);
}

TEST(RGG2DTest, n13_r01) {
    test_configuration(13, 0.1);
}

TEST(RGG2DTest, n50_r005) {
    test_configuration(50, 0.05);
}

TEST(RGG2DTest, n99_r003) {
    test_configuration(99, 0.03);
}

TEST(RGG2DTest, n127_r004) {
    test_configuration(127, 0.04);
}

// === Larger graphs (1000+ vertices) ===

TEST(RGG2DTest, n1000_r005) {
    test_configuration(1000, 0.05);
}

TEST(RGG2DTest, n1000_r01) {
    test_configuration(1000, 0.1);
}

TEST(RGG2DTest, n1000_r0125) {
    // Max radius boundary with large n
    test_configuration(1000, 0.125);
}

TEST(RGG2DTest, n2000_r003) {
    test_configuration(2000, 0.03);
}

TEST(RGG2DTest, n2000_r008) {
    test_configuration(2000, 0.08);
}

TEST(RGG2DTest, n3000_r005) {
    test_configuration(3000, 0.05);
}

TEST(RGG2DTest, n4000_r003) {
    test_configuration(4000, 0.03);
}

TEST(RGG2DTest, n1500_r007) {
    // Non-round node count
    test_configuration(1500, 0.07);
}

TEST(RGG2DTest, n2500_r004_seed42) {
    test_configuration(2500, 0.04, 42);
}

// === cells_per_dim transitions ===
// With k=64 (3,5,6 PEs), chunk_size=0.125
// cells_per_dim = floor(0.125 / r)
// r=0.125: cells_per_dim=1, r=0.063: cells_per_dim=1, r=0.062: cells_per_dim=2
// r=0.04: cells_per_dim=3, r=0.03: cells_per_dim=4

TEST(RGG2DTest, n32_r0062) {
    // cells_per_dim=2 boundary
    test_configuration(32, 0.062);
}

TEST(RGG2DTest, n32_r004) {
    // cells_per_dim=3
    test_configuration(32, 0.04);
}

TEST(RGG2DTest, n32_r003) {
    // cells_per_dim=4
    test_configuration(32, 0.03);
}
