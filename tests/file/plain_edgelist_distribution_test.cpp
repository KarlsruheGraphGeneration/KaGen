#include "kagen/context.h"
#include "kagen/generators/file/file_graph.h"

#include <gtest/gtest.h>
#include <mpi.h>

#include "tests/gather.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <string>
#include <unistd.h>
#include <vector>

using namespace kagen;

// Plain edgelist is a REQUIRES_REDISTRIBUTION reader (its Read() treats the given range as an arbitrary byte
// range, not vertex IDs, and FindNodeByEdge() is an unimplemented stub): balance-edges used to throw
// "not implemented" for it, and balance-edges-strict needs the same "read an arbitrary partition, then
// postprocess" treatment. This exercises both now succeeding, covering the one reader category the METIS/ParHIP
// tests in generic_file_generator_test.cpp cannot (they're FindNodeByEdge-capable and never need postprocessing
// for balance-edges).

namespace {
// A small hub-and-ring graph: vertex 0 connects to every other vertex (hub), plus a ring among the rest, so
// there's a real hub to (potentially) split under balance-edges-strict.
std::vector<std::pair<SInt, SInt>> BuildHubAndRingEdges(const SInt n) {
    std::vector<std::pair<SInt, SInt>> edges;
    for (SInt v = 1; v < n; ++v) {
        edges.emplace_back(0, v);
        edges.emplace_back(v, 0);
    }
    for (SInt v = 1; v < n; ++v) {
        const SInt next = 1 + (v % (n - 1));
        edges.emplace_back(v, next);
        edges.emplace_back(next, v);
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
}

void WritePlainEdgelist(const std::string& filename, const std::vector<std::pair<SInt, SInt>>& edges) {
    std::ofstream out(filename);
    for (const auto& [u, v]: edges) {
        out << u << " " << v << "\n";
    }
}

Graph ReadPlainEdgelist(const std::string& filename, const GraphDistribution distribution) {
    PGeneratorConfig config;
    config.input_graph.filename     = filename;
    config.input_graph.format       = FileFormat::PLAIN_EDGE_LIST;
    config.input_graph.distribution = distribution;

    PEID size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    FileGraphGenerator generator(config, rank, size);
    generator.Generate(GraphRepresentation::EDGE_LIST);
    generator.Finalize(MPI_COMM_WORLD);
    return generator.Take();
}
} // namespace

class PlainEdgelistDistributionTest : public ::testing::Test {
protected:
    void SetUp() override {
        const SInt n    = 2000;
        expected_edges_ = BuildHubAndRingEdges(n);

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // Unique per test process: ctest runs several core-count variants of this same binary concurrently, and
        // a shared filename would race (one variant's rank 0 could overwrite the file while another is reading
        // it). getpid() differs per rank (each MPI rank is its own OS process), so only rank 0 computes the
        // name and broadcasts it -- otherwise every other rank would look for a nonexistent file of its own.
        int name_length = 0;
        if (rank == 0) {
            filename_ = "plain_edgelist_distribution_test_input_" + std::to_string(getpid()) + ".txt";
            WritePlainEdgelist(filename_, expected_edges_);
            name_length = static_cast<int>(filename_.size());
        }
        MPI_Bcast(&name_length, 1, MPI_INT, 0, MPI_COMM_WORLD);
        filename_.resize(name_length);
        MPI_Bcast(filename_.data(), name_length, MPI_CHAR, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
    }

    void TearDown() override {
        MPI_Barrier(MPI_COMM_WORLD);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            std::remove(filename_.c_str());
        }
    }

    std::string                        filename_;
    std::vector<std::pair<SInt, SInt>> expected_edges_;
};

TEST_F(PlainEdgelistDistributionTest, BalanceVerticesSucceeds) {
    Graph local_graph  = ReadPlainEdgelist(filename_, GraphDistribution::BALANCE_VERTICES);
    Graph global_graph = kagen::testing::GatherEdgeLists(local_graph);
    std::sort(global_graph.edges.begin(), global_graph.edges.end());
    EXPECT_EQ(global_graph.edges, expected_edges_);
}

TEST_F(PlainEdgelistDistributionTest, BalanceEdgesSucceeds) {
    // Before this feature, this threw "not implemented" for any REQUIRES_REDISTRIBUTION reader.
    Graph local_graph  = ReadPlainEdgelist(filename_, GraphDistribution::BALANCE_EDGES);
    Graph global_graph = kagen::testing::GatherEdgeLists(local_graph);
    std::sort(global_graph.edges.begin(), global_graph.edges.end());
    EXPECT_EQ(global_graph.edges, expected_edges_);
}

TEST_F(PlainEdgelistDistributionTest, BalanceEdgesStrictSucceeds) {
    Graph local_graph  = ReadPlainEdgelist(filename_, GraphDistribution::BALANCE_EDGES_TRUE);
    Graph global_graph = kagen::testing::GatherEdgeLists(local_graph);
    std::sort(global_graph.edges.begin(), global_graph.edges.end());
    EXPECT_EQ(global_graph.edges, expected_edges_);
}

TEST_F(PlainEdgelistDistributionTest, BalanceEdgesStrictBalancesWithinOne) {
    Graph local_graph = ReadPlainEdgelist(filename_, GraphDistribution::BALANCE_EDGES_TRUE);

    int local_count = static_cast<int>(local_graph.edges.size());
    int min_count, max_count;
    MPI_Allreduce(&local_count, &min_count, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_count, &max_count, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    EXPECT_LE(max_count - min_count, 1);
}
