#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/geometric/geometric_3d.h"
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

void test_configuration(const SInt n, const double radius, const int seed = 1) {
    PGeneratorConfig config;
    config.n           = n;
    config.r           = radius;
    config.seed        = seed;
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

// All tests use:
// - n >= 8 to avoid MPI buffer aliasing in gather infrastructure
// - r <= 0.25 to stay within chunk_size/r >= 1 for k=64 (3,5,6 PEs)
// - Avoid very small r to prevent excessive cell counts

// === Small graphs ===

TEST(RGG3DTest, n8_r025) {
    test_configuration(8, 0.25);
}

TEST(RGG3DTest, n8_r015) {
    test_configuration(8, 0.15);
}

TEST(RGG3DTest, n8_r005) {
    test_configuration(8, 0.05);
}

// === Original tests ===

TEST(RGG3DTest, n16_r01) {
    test_configuration(16, 0.1);
}

TEST(RGG3DTest, n32_r0125) {
    test_configuration(32, 0.125);
}

TEST(RGG3DTest, n512_r001) {
    test_configuration(512, 0.01);
}

// === Moderate graphs with varying radius ===

TEST(RGG3DTest, n16_r025) {
    test_configuration(16, 0.25);
}

TEST(RGG3DTest, n16_r015) {
    test_configuration(16, 0.15);
}

TEST(RGG3DTest, n16_r005) {
    test_configuration(16, 0.05);
}

TEST(RGG3DTest, n32_r02) {
    test_configuration(32, 0.2);
}

TEST(RGG3DTest, n32_r01) {
    test_configuration(32, 0.1);
}

TEST(RGG3DTest, n32_r005) {
    test_configuration(32, 0.05);
}

TEST(RGG3DTest, n64_r015) {
    test_configuration(64, 0.15);
}

TEST(RGG3DTest, n64_r01) {
    test_configuration(64, 0.1);
}

TEST(RGG3DTest, n64_r005) {
    test_configuration(64, 0.05);
}

TEST(RGG3DTest, n100_r01) {
    test_configuration(100, 0.1);
}

TEST(RGG3DTest, n128_r008) {
    test_configuration(128, 0.08);
}

// === Larger graphs ===

TEST(RGG3DTest, n256_r005) {
    test_configuration(256, 0.05);
}

TEST(RGG3DTest, n256_r008) {
    test_configuration(256, 0.08);
}

// === Edge case: maximum allowed radius (r=0.25 boundary) ===
// With k=64 for 3D: chunks_per_dim=4, chunk_size=0.25
// cells_per_dim = floor(0.25 / r), needs to be >= 1

TEST(RGG3DTest, n32_r025) {
    test_configuration(32, 0.25);
}

TEST(RGG3DTest, n64_r025) {
    test_configuration(64, 0.25);
}

TEST(RGG3DTest, n16_r024) {
    // Just below boundary
    test_configuration(16, 0.24);
}

// === Edge case: very small radius (sparse, few/zero edges) ===

TEST(RGG3DTest, n32_r002) {
    test_configuration(32, 0.02);
}

TEST(RGG3DTest, n64_r002) {
    test_configuration(64, 0.02);
}

// === Different seeds ===

TEST(RGG3DTest, n32_r015_seed42) {
    test_configuration(32, 0.15, 42);
}

TEST(RGG3DTest, n32_r015_seed123) {
    test_configuration(32, 0.15, 123);
}

TEST(RGG3DTest, n64_r008_seed7) {
    test_configuration(64, 0.08, 7);
}

TEST(RGG3DTest, n64_r008_seed999) {
    test_configuration(64, 0.08, 999);
}

// === Non-power-of-two node counts ===

TEST(RGG3DTest, n9_r015) {
    test_configuration(9, 0.15);
}

TEST(RGG3DTest, n13_r01) {
    test_configuration(13, 0.1);
}

TEST(RGG3DTest, n50_r008) {
    test_configuration(50, 0.08);
}

TEST(RGG3DTest, n99_r005) {
    test_configuration(99, 0.05);
}

TEST(RGG3DTest, n127_r006) {
    test_configuration(127, 0.06);
}

// === Larger graphs (1000+ vertices) ===

TEST(RGG3DTest, n1000_r008) {
    test_configuration(1000, 0.08);
}

TEST(RGG3DTest, n1000_r015) {
    test_configuration(1000, 0.15);
}

TEST(RGG3DTest, n1000_r025) {
    // Max radius boundary with large n
    test_configuration(1000, 0.25);
}

TEST(RGG3DTest, n2000_r005) {
    test_configuration(2000, 0.05);
}

TEST(RGG3DTest, n2000_r01) {
    test_configuration(2000, 0.1);
}

TEST(RGG3DTest, n3000_r008) {
    test_configuration(3000, 0.08);
}

TEST(RGG3DTest, n4000_r005) {
    test_configuration(4000, 0.05);
}

TEST(RGG3DTest, n1500_r01) {
    // Non-round node count
    test_configuration(1500, 0.1);
}

TEST(RGG3DTest, n2500_r006_seed42) {
    test_configuration(2500, 0.06, 42);
}

// === cells_per_dim transitions ===
// With k=64 (3,5,6 PEs), chunks_per_dim=4, chunk_size=0.25
// cells_per_dim = floor(0.25 / r)
// r=0.25: cells_per_dim=1, r=0.13: cells_per_dim=1, r=0.125: cells_per_dim=2
// r=0.08: cells_per_dim=3, r=0.06: cells_per_dim=4

TEST(RGG3DTest, n32_r0125_cells2) {
    // cells_per_dim=2 boundary
    test_configuration(32, 0.125);
}

TEST(RGG3DTest, n32_r008) {
    // cells_per_dim=3
    test_configuration(32, 0.08);
}

TEST(RGG3DTest, n32_r006) {
    // cells_per_dim=4
    test_configuration(32, 0.06);
}
