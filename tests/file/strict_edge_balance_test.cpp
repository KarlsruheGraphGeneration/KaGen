#include "kagen/context.h"
#include "kagen/generators/file/file_graph.h"

#include <gtest/gtest.h>
#include <mpi.h>
#include <unistd.h>

#include "tests/gather.h"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using namespace kagen;

// Exercises the direct strict edge-balanced read path (GraphReader::ReadStrictEdgeRange), which serves each PE's
// exact [from_edge, to_edge) slice straight from disk with no redistribution. Runs over both ParHIP (native CSR
// on disk) and METIS (CSR built locally from the tail-sorted slice) so the split-aware CSR output is covered --
// it cannot go through the generic GatherCSR helper, whose xadj.size()-1 == vertex_range span assumption is
// deliberately broken by a split boundary vertex (which owns a partial row on two PEs).

namespace {
// A hub-and-ring graph: vertex 0 is adjacent to every other vertex (a high-degree hub that will be split across
// PEs under strict edge balance once its degree exceeds a PE's fair share), plus a ring among the rest. Returned
// sorted and deduplicated, i.e. a clean simple graph.
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

// Writes an unweighted, 64-bit ParHIP file (version 3) from a sorted edge list.
void WriteParhip(const std::string& filename, const SInt n, const std::vector<std::pair<SInt, SInt>>& edges) {
    std::vector<std::uint64_t> xadj(n + 1, 0);
    for (const auto& [u, v]: edges) {
        ++xadj[u + 1];
    }
    for (SInt v = 0; v < n; ++v) {
        xadj[v + 1] += xadj[v];
    }
    std::vector<std::uint64_t> adjncy;
    adjncy.reserve(edges.size());
    for (const auto& [u, v]: edges) {
        adjncy.push_back(static_cast<std::uint64_t>(v));
    }

    const std::uint64_t        m             = edges.size();
    const std::uint64_t        header_bytes  = 3 * sizeof(std::uint64_t);
    const std::uint64_t        offsets_bytes = (n + 1) * sizeof(std::uint64_t);
    std::vector<std::uint64_t> offsets(n + 1);
    for (SInt v = 0; v <= n; ++v) {
        offsets[v] = header_bytes + offsets_bytes + xadj[v] * sizeof(std::uint64_t);
    }

    std::ofstream       out(filename, std::ios::binary);
    const std::uint64_t version = 3;
    out.write(reinterpret_cast<const char*>(&version), sizeof(version));
    out.write(reinterpret_cast<const char*>(&n), sizeof(std::uint64_t));
    out.write(reinterpret_cast<const char*>(&m), sizeof(std::uint64_t));
    out.write(reinterpret_cast<const char*>(offsets.data()), offsets.size() * sizeof(std::uint64_t));
    out.write(reinterpret_cast<const char*>(adjncy.data()), adjncy.size() * sizeof(std::uint64_t));
}

// Writes an unweighted METIS text file from a sorted edge list.
void WriteMetis(const std::string& filename, const SInt n, const std::vector<std::pair<SInt, SInt>>& edges) {
    std::vector<std::vector<SInt>> adj(n);
    for (const auto& [u, v]: edges) {
        adj[u].push_back(v);
    }
    std::ofstream out(filename);
    out << n << " " << edges.size() / 2 << "\n";
    for (SInt u = 0; u < n; ++u) {
        for (std::size_t i = 0; i < adj[u].size(); ++i) {
            out << (adj[u][i] + 1) << (i + 1 < adj[u].size() ? " " : "");
        }
        out << "\n";
    }
}

Graph ReadFile(
    const std::string& filename, const FileFormat format, const GraphDistribution distribution,
    const GraphRepresentation rep) {
    PGeneratorConfig config;
    config.input_graph.filename     = filename;
    config.input_graph.format       = format;
    config.input_graph.distribution = distribution;

    PEID size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    FileGraphGenerator generator(config, rank, size);
    generator.Generate(rep);
    generator.Finalize(MPI_COMM_WORLD);
    return generator.Take();
}

// Turns a (possibly split-aware) local CSR fragment into its local edge list, using the row space
// [PhysicalVertexRange().first, PhysicalVertexRange().second) that xadj actually indexes.
Graph LocalCsrToEdgeList(const Graph& csr) {
    Graph out;
    out.representation         = GraphRepresentation::EDGE_LIST;
    const VertexRange rows     = csr.PhysicalVertexRange();
    const SInt        num_rows = rows.second - rows.first;
    for (SInt i = 0; i < num_rows; ++i) {
        const SInt tail = rows.first + i;
        for (SInt e = csr.xadj[i]; e < csr.xadj[i + 1]; ++e) {
            out.edges.emplace_back(tail, csr.adjncy[e]);
        }
    }
    return out;
}

void ExpectBalancedWithinOne(const SInt local_edges) {
    int local = static_cast<int>(local_edges);
    int min_count, max_count;
    MPI_Allreduce(&local, &min_count, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local, &max_count, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    EXPECT_LE(max_count - min_count, 1);
}
} // namespace

class StrictEdgeBalanceTest : public ::testing::TestWithParam<FileFormat> {
protected:
    void SetUp() override {
        n_              = 500;
        expected_edges_ = BuildHubAndRingEdges(n_);
        format_         = GetParam();

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int name_length = 0;
        if (rank == 0) {
            const char* ext = format_ == FileFormat::PARHIP ? ".parhip" : ".metis";
            filename_       = "strict_edge_balance_test_input_" + std::to_string(getpid()) + ext;
            if (format_ == FileFormat::PARHIP) {
                WriteParhip(filename_, n_, expected_edges_);
            } else {
                WriteMetis(filename_, n_, expected_edges_);
            }
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

    Graph Read(const GraphRepresentation rep) {
        return ReadFile(filename_, format_, GraphDistribution::BALANCE_EDGES_TRUE, rep);
    }

    SInt                               n_ = 0;
    FileFormat                         format_;
    std::string                        filename_;
    std::vector<std::pair<SInt, SInt>> expected_edges_;
};

INSTANTIATE_TEST_SUITE_P(Formats, StrictEdgeBalanceTest, ::testing::Values(FileFormat::PARHIP, FileFormat::METIS));

TEST_P(StrictEdgeBalanceTest, EdgeListPreservesAllEdges) {
    Graph global = kagen::testing::GatherEdgeLists(Read(GraphRepresentation::EDGE_LIST));
    std::sort(global.edges.begin(), global.edges.end());
    EXPECT_EQ(global.edges, expected_edges_);
}

TEST_P(StrictEdgeBalanceTest, EdgeListBalancesWithinOne) {
    ExpectBalancedWithinOne(Read(GraphRepresentation::EDGE_LIST).edges.size());
}

TEST_P(StrictEdgeBalanceTest, EdgeListGapFreeVertexRange) {
    Graph      local   = Read(GraphRepresentation::EDGE_LIST);
    const SInt local_n = local.vertex_range.second - local.vertex_range.first;
    SInt       total_n = local_n;
    MPI_Allreduce(MPI_IN_PLACE, &total_n, 1, KAGEN_MPI_SINT, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(total_n, n_); // gap-free partition of [0, n)
}

TEST_P(StrictEdgeBalanceTest, CsrIsInternallyConsistent) {
    Graph local = Read(GraphRepresentation::CSR);
    ASSERT_EQ(local.representation, GraphRepresentation::CSR);
    ASSERT_FALSE(local.xadj.empty());
    // xadj indexes PhysicalVertexRange() (the physically-present, possibly-overlapping row space), not
    // vertex_range (the gap-free ownership range) -- they differ by one row on a left-partial split boundary.
    const VertexRange rows = local.PhysicalVertexRange();
    EXPECT_EQ(local.xadj.size() - 1, rows.second - rows.first);
    EXPECT_EQ(local.xadj.back(), local.adjncy.size());
}

TEST_P(StrictEdgeBalanceTest, CsrBalancesWithinOne) {
    ExpectBalancedWithinOne(Read(GraphRepresentation::CSR).adjncy.size());
}

TEST_P(StrictEdgeBalanceTest, CsrPreservesAllEdges) {
    Graph global = kagen::testing::GatherEdgeLists(LocalCsrToEdgeList(Read(GraphRepresentation::CSR)));
    std::sort(global.edges.begin(), global.edges.end());
    // Each directed edge is stored on exactly one PE (in the row of its tail), so gathering every PE's rows --
    // even the partial rows of a split hub -- reconstructs the original edge set exactly.
    EXPECT_EQ(global.edges, expected_edges_);
}

TEST_P(StrictEdgeBalanceTest, CsrVertexRangeIsGapFree) {
    Graph      local   = Read(GraphRepresentation::CSR);
    const SInt local_n = local.NumberOfLocalVertices();

    SInt total_n = local_n;
    MPI_Allreduce(MPI_IN_PLACE, &total_n, 1, KAGEN_MPI_SINT, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(total_n, n_); // gap-free partition of [0, n), even though PhysicalVertexRange() overlaps at splits
}

TEST_P(StrictEdgeBalanceTest, SplitMetadataIsCollectivelyConsistent) {
    Graph local = Read(GraphRepresentation::CSR);

    int has_split = local.has_split_vertices ? 1 : 0;
    int min_split, max_split;
    MPI_Allreduce(&has_split, &min_split, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&has_split, &max_split, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    EXPECT_EQ(min_split, max_split);

    // right_partial_vertex is credited into vertex_range (its vertex is the last one in it); left_partial_vertex
    // is ceded to the lower-rank neighbor and so sits one below vertex_range.first (= PhysicalVertexRange().first).
    if (local.right_partial_vertex) {
        EXPECT_EQ(local.right_partial_vertex->vertex, local.vertex_range.second - 1);
    }
    if (local.left_partial_vertex) {
        EXPECT_EQ(local.left_partial_vertex->vertex, local.PhysicalVertexRange().first);
    }
    for (const auto& pv: {local.left_partial_vertex, local.right_partial_vertex}) {
        if (!pv) {
            continue;
        }
        EXPECT_LE(pv->local_offset + pv->local_count, static_cast<SInt>(local.adjncy.size()));
    }
}

TEST_P(StrictEdgeBalanceTest, DirectMatchesVertexBalancedEdgeSet) {
    Graph strict = kagen::testing::GatherEdgeLists(Read(GraphRepresentation::EDGE_LIST));
    Graph vbal   = kagen::testing::GatherEdgeLists(
        ReadFile(filename_, format_, GraphDistribution::BALANCE_VERTICES, GraphRepresentation::EDGE_LIST));
    std::sort(strict.edges.begin(), strict.edges.end());
    std::sort(vbal.edges.begin(), vbal.edges.end());
    EXPECT_EQ(strict.edges, vbal.edges);
}
