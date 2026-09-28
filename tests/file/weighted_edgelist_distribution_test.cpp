#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/file/file_graph.h"

#include <gtest/gtest.h>
#include <mpi.h>

#include "tests/gather.h"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <random>
#include <string>
#include <unistd.h>
#include <vector>

using namespace kagen;

// Covers the direct edge-balanced read of the (weighted) binary edge list, the one edge-list reader eligible for
// it (fixed-width records => edge count known, random-access slices). BALANCE_EDGES heals the strict slice to
// whole-vertex ownership via a neighbor exchange that must carry edge weights; BALANCE_EDGES_TRUE keeps the exact
// ±1 slice. Unsorted input falls back to redistribution, which rejects weighted graphs.

namespace {
SSInt ExpectedWeight(const SInt u, const SInt v) {
    return static_cast<SSInt>((u * 31 + v) % 997 + 1); // deterministic, positive, distinguishes edges
}

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

// Writes a weighted binary edge list: per edge, int64 from, int64 to, int64 weight (matches the default
// vtx_width/adjwgt_width of 64 bits).
void WriteWeightedBinaryEdgelist(const std::string& filename, const std::vector<std::pair<SInt, SInt>>& edges) {
    std::ofstream out(filename, std::ios::binary);
    for (const auto& [u, v]: edges) {
        const std::int64_t from = u, to = v, weight = ExpectedWeight(u, v);
        out.write(reinterpret_cast<const char*>(&from), sizeof(from));
        out.write(reinterpret_cast<const char*>(&to), sizeof(to));
        out.write(reinterpret_cast<const char*>(&weight), sizeof(weight));
    }
}

Graph Read(const std::string& filename, const GraphDistribution distribution) {
    PGeneratorConfig config;
    config.input_graph.filename     = filename;
    config.input_graph.format       = FileFormat::WEIGHTED_BINARY_EDGE_LIST;
    config.input_graph.distribution = distribution;

    PEID size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    FileGraphGenerator generator(config, rank, size);
    generator.Generate(GraphRepresentation::EDGE_LIST);
    generator.Finalize(MPI_COMM_WORLD);
    return generator.Take();
}

// Verifies locally that every edge still carries its correct weight (catches any weight/edge desync in the heal
// exchange), then returns the local edge count.
SInt CheckLocalWeightsAndCount(const Graph& g) {
    EXPECT_EQ(g.edges.size(), g.edge_weights.size());
    for (std::size_t i = 0; i < g.edges.size(); ++i) {
        EXPECT_EQ(g.edge_weights[i], ExpectedWeight(g.edges[i].first, g.edges[i].second));
    }
    return static_cast<SInt>(g.edges.size());
}
} // namespace

class WeightedEdgelistDistributionTest : public ::testing::Test {
protected:
    void SetUp() override {
        expected_edges_ = BuildHubAndRingEdges(500);

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        int name_length = 0;
        if (rank == 0) {
            base_        = "weighted_edgelist_distribution_test_" + std::to_string(getpid());
            sorted_file_ = base_ + "_sorted.bin";
            WriteWeightedBinaryEdgelist(sorted_file_, expected_edges_);

            auto shuffled = expected_edges_;
            std::mt19937 rng(12345);
            std::shuffle(shuffled.begin(), shuffled.end(), rng);
            unsorted_file_ = base_ + "_unsorted.bin";
            WriteWeightedBinaryEdgelist(unsorted_file_, shuffled);

            name_length = static_cast<int>(base_.size());
        }
        MPI_Bcast(&name_length, 1, MPI_INT, 0, MPI_COMM_WORLD);
        base_.resize(name_length);
        MPI_Bcast(base_.data(), name_length, MPI_CHAR, 0, MPI_COMM_WORLD);
        sorted_file_   = base_ + "_sorted.bin";
        unsorted_file_ = base_ + "_unsorted.bin";
        MPI_Barrier(MPI_COMM_WORLD);
    }

    void TearDown() override {
        MPI_Barrier(MPI_COMM_WORLD);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            std::remove(sorted_file_.c_str());
            std::remove(unsorted_file_.c_str());
        }
    }

    std::string                        base_, sorted_file_, unsorted_file_;
    std::vector<std::pair<SInt, SInt>> expected_edges_;
};

TEST_F(WeightedEdgelistDistributionTest, BalanceEdgesDirectPreservesEdgesAndWeights) {
    Graph local = Read(sorted_file_, GraphDistribution::BALANCE_EDGES);
    CheckLocalWeightsAndCount(local);

    // Vertex-atomic: no vertex is split across PEs.
    int has_split = local.has_split_vertices ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &has_split, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    EXPECT_EQ(has_split, 0);

    Graph global = kagen::testing::GatherEdgeLists(local);
    std::sort(global.edges.begin(), global.edges.end());
    EXPECT_EQ(global.edges, expected_edges_);
}

TEST_F(WeightedEdgelistDistributionTest, BalanceEdgesStrictPreservesEdgesAndWeights) {
    Graph local = Read(sorted_file_, GraphDistribution::BALANCE_EDGES_TRUE);
    CheckLocalWeightsAndCount(local);

    int local_count = static_cast<int>(local.edges.size());
    int min_count, max_count;
    MPI_Allreduce(&local_count, &min_count, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_count, &max_count, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    EXPECT_LE(max_count - min_count, 1); // strict ±1

    Graph global = kagen::testing::GatherEdgeLists(local);
    std::sort(global.edges.begin(), global.edges.end());
    EXPECT_EQ(global.edges, expected_edges_);
}

TEST_F(WeightedEdgelistDistributionTest, UnsortedBalanceEdgesFallsBackAndRejectsWeights) {
    // Unsorted input cannot use the direct read, so it falls back to RedistributeEdgesBalanced, which does not
    // carry weights and therefore rejects weighted input (unchanged behavior).
    EXPECT_THROW(Read(unsorted_file_, GraphDistribution::BALANCE_EDGES), ConfigurationError);
}

TEST_F(WeightedEdgelistDistributionTest, UnsortedBalanceEdgesStrictFallsBackAndRejectsWeights) {
    EXPECT_THROW(Read(unsorted_file_, GraphDistribution::BALANCE_EDGES_TRUE), ConfigurationError);
}
