#ifndef KAGEN_NOMPI
    #include "kagen/comm/mpi_comm.h"
#else
    #include "kagen/comm/seq_comm.h"
#endif
#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/hyperbolic/hyperbolic.h"

#include <gtest/gtest.h>

#include "tests/gather.h"
#include "tests/geometric/utils.h"

using namespace kagen;

namespace {
void validate_graph(const Graph& local_graph, const PGeneratorConfig& config, kagen::Comm& comm) {
    Graph global_graph = kagen::testing::GatherEdgeLists(local_graph, comm);

    PEID rank = comm.Rank();

    if (rank == 0) {
        EXPECT_EQ(config.n, global_graph.coordinates.first.size());

        // Creating the correct edge list as a test instance
        std::vector<std::pair<SInt, SInt>> expected_edges =
            kagen::testing::CreateExpectedHyperbolicEdges<HPFloat>(config, global_graph);

        // Sorting both lists before comparing them
        std::sort(global_graph.edges.begin(), global_graph.edges.end());
        std::sort(expected_edges.begin(), expected_edges.end());

        EXPECT_EQ(global_graph.edges, expected_edges);
    }
}

void test_configuration(const SInt n, const double avg_degree, const double plexp, const int seed = 1) {
    PGeneratorConfig config;
    config.n           = n;
    config.avg_degree  = avg_degree;
    config.plexp       = plexp;
    config.seed        = seed;
    config.hp_floats   = true;
    config.coordinates = true;

    HyperbolicFactory factory;
#ifndef KAGEN_NOMPI
    kagen::MPIComm mpi_world(MPI_COMM_WORLD);
#else
    kagen::SeqComm mpi_world;
#endif
    PEID size = mpi_world.Size();
    PEID rank = mpi_world.Rank();
    config         = factory.NormalizeParameters(config, rank, size, false);
    auto generator = factory.Create(config, rank, size);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    generator->Finalize(mpi_world);
    const auto local_graph = generator->Take();

    validate_graph(local_graph, config, mpi_world);
}
} // namespace

// All tests use n >= 8 to avoid MPI buffer aliasing in gather infrastructure.
// Parameters: n, avg_degree, plexp (gamma)

// === Original tests ===

TEST(HyperbolicTest, n16_d4_g3) {
    test_configuration(16, 4.0, 3.0);
}

TEST(HyperbolicTest, n512_d4_g3) {
    test_configuration(512, 4.0, 3.0);
}

// === Known problematic case (vertex_range inconsistencies) ===

TEST(HyperbolicTest, n1024_d16_g3) {
    test_configuration(1024, 16.0, 3.0);
}

// === Small graphs ===

TEST(HyperbolicTest, n8_d2_g3) {
    test_configuration(8, 2.0, 3.0);
}

TEST(HyperbolicTest, n8_d4_g3) {
    test_configuration(8, 4.0, 3.0);
}

TEST(HyperbolicTest, n16_d2_g3) {
    test_configuration(16, 2.0, 3.0);
}

TEST(HyperbolicTest, n32_d4_g3) {
    test_configuration(32, 4.0, 3.0);
}

// === Varying avg_degree ===

TEST(HyperbolicTest, n64_d2_g3) {
    // Low average degree (sparse)
    test_configuration(64, 2.0, 3.0);
}

TEST(HyperbolicTest, n64_d4_g3) {
    test_configuration(64, 4.0, 3.0);
}

TEST(HyperbolicTest, n64_d8_g3) {
    // Higher average degree
    test_configuration(64, 8.0, 3.0);
}

TEST(HyperbolicTest, n128_d4_g3) {
    test_configuration(128, 4.0, 3.0);
}

TEST(HyperbolicTest, n128_d8_g3) {
    test_configuration(128, 8.0, 3.0);
}

TEST(HyperbolicTest, n128_d16_g3) {
    // High average degree
    test_configuration(128, 16.0, 3.0);
}

TEST(HyperbolicTest, n256_d4_g3) {
    test_configuration(256, 4.0, 3.0);
}

TEST(HyperbolicTest, n256_d16_g3) {
    test_configuration(256, 16.0, 3.0);
}

// === Varying plexp (gamma / power law exponent) ===

TEST(HyperbolicTest, n64_d4_g22) {
    // Low gamma: flatter degree distribution, more hub nodes
    test_configuration(64, 4.0, 2.2);
}

TEST(HyperbolicTest, n64_d4_g25) {
    test_configuration(64, 4.0, 2.5);
}

TEST(HyperbolicTest, n64_d4_g35) {
    // Higher gamma: steeper power law, fewer hubs
    test_configuration(64, 4.0, 3.5);
}

TEST(HyperbolicTest, n64_d4_g5) {
    // Very high gamma: very steep power law
    test_configuration(64, 4.0, 5.0);
}

TEST(HyperbolicTest, n128_d4_g22) {
    test_configuration(128, 4.0, 2.2);
}

TEST(HyperbolicTest, n128_d4_g25) {
    test_configuration(128, 4.0, 2.5);
}

TEST(HyperbolicTest, n256_d4_g25) {
    test_configuration(256, 4.0, 2.5);
}

TEST(HyperbolicTest, n256_d4_g5) {
    test_configuration(256, 4.0, 5.0);
}

// === Combined: varying both avg_degree and plexp ===

TEST(HyperbolicTest, n64_d2_g25) {
    // Sparse + low gamma
    test_configuration(64, 2.0, 2.5);
}

TEST(HyperbolicTest, n64_d8_g25) {
    // Dense + low gamma
    test_configuration(64, 8.0, 2.5);
}

TEST(HyperbolicTest, n64_d2_g5) {
    // Sparse + high gamma
    test_configuration(64, 2.0, 5.0);
}

TEST(HyperbolicTest, n64_d8_g5) {
    // Dense + high gamma
    test_configuration(64, 8.0, 5.0);
}

// === Larger graphs ===

TEST(HyperbolicTest, n1000_d4_g3) {
    test_configuration(1000, 4.0, 3.0);
}

TEST(HyperbolicTest, n1000_d8_g25) {
    test_configuration(1000, 8.0, 2.5);
}

TEST(HyperbolicTest, n2000_d4_g3) {
    test_configuration(2000, 4.0, 3.0);
}

TEST(HyperbolicTest, n2000_d8_g3) {
    test_configuration(2000, 8.0, 3.0);
}

TEST(HyperbolicTest, n2000_d16_g3) {
    test_configuration(2000, 16.0, 3.0);
}

TEST(HyperbolicTest, n3000_d4_g25) {
    test_configuration(3000, 4.0, 2.5);
}

// === Different seeds ===

TEST(HyperbolicTest, n64_d4_g3_seed42) {
    test_configuration(64, 4.0, 3.0, 42);
}

TEST(HyperbolicTest, n64_d4_g3_seed123) {
    test_configuration(64, 4.0, 3.0, 123);
}

TEST(HyperbolicTest, n256_d4_g3_seed7) {
    test_configuration(256, 4.0, 3.0, 7);
}

TEST(HyperbolicTest, n256_d8_g25_seed999) {
    test_configuration(256, 8.0, 2.5, 999);
}

// === Non-power-of-two node counts ===

TEST(HyperbolicTest, n9_d4_g3) {
    test_configuration(9, 4.0, 3.0);
}

TEST(HyperbolicTest, n13_d4_g3) {
    test_configuration(13, 4.0, 3.0);
}

TEST(HyperbolicTest, n50_d4_g3) {
    test_configuration(50, 4.0, 3.0);
}

TEST(HyperbolicTest, n99_d4_g3) {
    test_configuration(99, 4.0, 3.0);
}

TEST(HyperbolicTest, n127_d4_g3) {
    test_configuration(127, 4.0, 3.0);
}

TEST(HyperbolicTest, n500_d4_g3) {
    test_configuration(500, 4.0, 3.0);
}

TEST(HyperbolicTest, n1500_d8_g3) {
    test_configuration(1500, 8.0, 3.0);
}
