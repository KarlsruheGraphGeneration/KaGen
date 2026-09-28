#include "kagen/context.h"
#include "kagen/in_memory_facade.h"
#include "kagen/kagen.h"
#include "kagen/tools/utils.h"

#include <gtest/gtest.h>
#include <kagen/tools/random_permutation.h>

#include "../gather.h"
#include "../utils.h"
#include "factories.h"
#include <unordered_map>

using namespace kagen;

template <typename Transformer>
void transform_vertices(std::vector<kagen::testing::SrcDstEdgeWeight>& weighted_edges, Transformer&& op) {
    for (auto& [src, dst, _]: weighted_edges) {
        src = op(src);
        dst = op(dst);
    }
}

void test_distribution_of_edges(SSInt n_from_config, const Graph& permuted_graph) {
    int      size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &size);

    // vertices should be split equally
    const SInt n_local = permuted_graph.NumberOfLocalVertices();
    EXPECT_LE(n_local, n_from_config / size + 2);

    auto       permuted_graph_edges = kagen::testing::ConvertToWeightedEdgelist(permuted_graph);
    const auto vertex_range         = permuted_graph.vertex_range;
    for (const auto& [src, _, __]: permuted_graph_edges) {
        EXPECT_GE(src, vertex_range.first);
        EXPECT_LT(src, vertex_range.second);
    }
}

void test_equality_of_permuted_graph(SSInt n_from_config, const Graph& permuted_graph, const Graph& graph) {
    auto gathered_graph                = kagen::testing::GatherGraph(graph);
    auto gathered_graph_edges          = kagen::testing::ConvertToWeightedEdgelist(gathered_graph);
    auto gathered_permuted_graph       = kagen::testing::GatherGraph(permuted_graph);
    auto gathered_permuted_graph_edges = kagen::testing::ConvertToWeightedEdgelist(gathered_permuted_graph);

    // we do not want the 'random' permutation to be the ID
    EXPECT_NE(gathered_graph_edges, gathered_permuted_graph_edges);

    // apply inverse permutation to permuted graph
    auto permutator = random_permutation::FeistelPseudoRandomPermutation::buildPermutation(n_from_config - 1, 0);
    auto permute    = [&permutator](SInt v) {
        return permutator.finv(v);
    };
    transform_vertices(gathered_permuted_graph_edges, permute);
    std::sort(gathered_graph_edges.begin(), gathered_graph_edges.end());
    std::sort(gathered_permuted_graph_edges.begin(), gathered_permuted_graph_edges.end());

    EXPECT_EQ(gathered_graph_edges, gathered_permuted_graph_edges);

    // test whether vertex weights were permuted correctly
    std::vector<SSInt> vertex_weights(gathered_permuted_graph.vertex_weights.size(), 0);
    for (std::size_t i = 0; i < vertex_weights.size(); ++i) {
        vertex_weights[permute(i)] = gathered_permuted_graph.vertex_weights[i];
    }
    EXPECT_EQ(gathered_graph.vertex_weights, vertex_weights);
}

TEST(GraphPermutation, check_applied_permutation_rgg2d_edgelist_with_vtx_weights) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 64'000;

    kagen::KaGen generator(comm);
    generator.UseEdgeListRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    generator.ConfigureVertexWeightGeneration(VertexWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateRGG2D_NM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateRGG2D_NM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}

TEST(GraphPermutation, check_applied_permutation_rgg2d_csr_with_vtx_weights) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 64'000;

    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    generator.ConfigureVertexWeightGeneration(VertexWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateRGG2D_NM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateRGG2D_NM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}

TEST(GraphPermutation, check_applied_permutation_rgg2d_edgelist_sparse_with_vtx_weights) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 4'000;

    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    generator.ConfigureVertexWeightGeneration(VertexWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateRGG2D_NM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateRGG2D_NM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}

TEST(GraphPermutation, check_applied_permutation_rgg2d_csr_sparse_with_vtx_weights) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 4'000;

    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    generator.ConfigureVertexWeightGeneration(VertexWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateRGG2D_NM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateRGG2D_NM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}

TEST(GraphPermutation, check_applied_permutation_gnm_edgelist_sparse_with_vtx_weights) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 4'000;

    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    generator.ConfigureVertexWeightGeneration(VertexWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateUndirectedGNM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateUndirectedGNM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}

TEST(GraphPermutation, check_applied_permutation_gnm_csr_sparse_with_vtx_weights) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 4'000;

    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateUndirectedGNM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateUndirectedGNM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}

TEST(GraphPermutation, check_applied_permutation_rgg2d_edgelist) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 64'000;

    kagen::KaGen generator(comm);
    generator.UseEdgeListRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateRGG2D_NM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateRGG2D_NM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}

TEST(GraphPermutation, check_applied_permutation_rgg2d_csr) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 64'000;

    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateRGG2D_NM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateRGG2D_NM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}

TEST(GraphPermutation, check_applied_permutation_rgg2d_edgelist_sparse) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 4'000;

    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateRGG2D_NM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateRGG2D_NM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}

TEST(GraphPermutation, check_applied_permutation_rgg2d_csr_sparse) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 4'000;

    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateRGG2D_NM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateRGG2D_NM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}

TEST(GraphPermutation, check_applied_permutation_gnm_edgelist_sparse) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 4'000;

    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateUndirectedGNM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateUndirectedGNM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}

TEST(GraphPermutation, check_applied_permutation_gnm_csr_sparse) {
    MPI_Comm comm = MPI_COMM_WORLD;

    const SInt n = 16'000;
    const SInt m = 4'000;

    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(EdgeWeightGeneratorType::UNIFORM_RANDOM, 1, 100'000'000);
    auto graph = generator.GenerateUndirectedGNM(n, m);
    generator.EnableVertexPermutation();
    auto permuted_graph = generator.GenerateUndirectedGNM(n, m);

    test_distribution_of_edges(n, permuted_graph);
    test_equality_of_permuted_graph(n, permuted_graph, graph);
}

// ---------------------------------------------------------------------------
// Permuting an edge-balanced graph.
//
// A permutation only relabels vertices -- it does not preserve how many *edges* a PE holds, since each PE's new
// vertex set is an essentially random subset of the old one. So whatever balance --redistribution established
// has to be re-established in the permuted ID space, otherwise --permute silently downgrades it to a plain
// vertex-balanced distribution. These tests pin down both halves of that: the permuted graph is still the same
// graph (its edge set is preserved under the inverse permutation), and it is still balanced.
//
// balance-edges-strict additionally splits a vertex's own edges across PEs, which means the permutation has to
// cope with a graph whose CSR row space is Graph::PhysicalVertexRange() rather than vertex_range, and whose
// vertices can arrive at their new owner in several partial pieces.
// ---------------------------------------------------------------------------
namespace {
using EdgeVec = std::vector<std::pair<SInt, SInt>>;

// Split-aware, unlike tests/utils.h's ConvertToWeightedEdgelist: CSR rows are indexed by PhysicalVertexRange(),
// which includes a left-partial boundary vertex that vertex_range excludes.
EdgeVec LocalEdges(const Graph& graph) {
    if (graph.representation == GraphRepresentation::EDGE_LIST) {
        return EdgeVec(graph.edges.begin(), graph.edges.end());
    }

    const VertexRange row_range = graph.PhysicalVertexRange();
    EdgeVec           edges;
    edges.reserve(graph.adjncy.size());
    for (std::size_t i = 0; i + 1 < graph.xadj.size(); ++i) {
        for (SInt e = graph.xadj[i]; e < graph.xadj[i + 1]; ++e) {
            edges.emplace_back(row_range.first + i, graph.adjncy[e]);
        }
    }
    return edges;
}

// Collects all edges on the root PE, sorted, so two distributions of the same graph can be compared directly.
EdgeVec GatherEdges(const EdgeVec& local, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::vector<SInt> flat;
    flat.reserve(local.size() * 2);
    for (const auto& [src, dst]: local) {
        flat.push_back(src);
        flat.push_back(dst);
    }

    const int        local_count = static_cast<int>(flat.size());
    std::vector<int> counts(static_cast<std::size_t>(size));
    MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);

    std::vector<int> displs(static_cast<std::size_t>(size), 0);
    int              total = 0;
    for (PEID pe = 0; pe < size; ++pe) {
        displs[static_cast<std::size_t>(pe)] = total;
        total += counts[static_cast<std::size_t>(pe)];
    }

    std::vector<SInt> received(rank == 0 ? static_cast<std::size_t>(total) : 0);
    MPI_Gatherv(
        flat.data(), local_count, KAGEN_MPI_SINT, received.data(), counts.data(), displs.data(), KAGEN_MPI_SINT, 0,
        comm);

    EdgeVec edges;
    edges.reserve(received.size() / 2);
    for (std::size_t i = 0; i + 1 < received.size(); i += 2) {
        edges.emplace_back(received[i], received[i + 1]);
    }
    std::sort(edges.begin(), edges.end());
    return edges;
}

SInt MaxDegree(const Graph& graph, MPI_Comm comm) {
    std::unordered_map<SInt, SInt> degrees;
    for (const auto& [src, dst]: LocalEdges(graph)) {
        ++degrees[src];
    }
    SInt max_degree = 0;
    for (const auto& [vertex, degree]: degrees) {
        max_degree = std::max(max_degree, degree);
    }
    MPI_Allreduce(MPI_IN_PLACE, &max_degree, 1, KAGEN_MPI_SINT, MPI_MAX, comm);
    return max_degree;
}

PGeneratorConfig EdgeBalancedConfig(const std::string& redistribution) {
    PGeneratorConfig config =
        CreateConfigFromString("gnm_undirected;N=10;M=14;redistribution=" + redistribution, PGeneratorConfig{});
    config.quiet            = true;
    config.statistics_level = StatisticsLevel::NONE;
    return config;
}

// Every vertex of [0, n) is owned by exactly one PE, in ascending rank order.
void ExpectGapFreeVertexRanges(const Graph& graph, const SInt n, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const std::vector<VertexRange> ranges = AllgatherVertexRange(graph.vertex_range, comm);
    if (rank != 0) {
        return;
    }
    EXPECT_EQ(ranges.front().first, 0u);
    EXPECT_EQ(ranges.back().second, n);
    for (PEID pe = 0; pe + 1 < size; ++pe) {
        const auto& range = ranges[static_cast<std::size_t>(pe)];
        EXPECT_LE(range.first, range.second) << "PE " << pe << " has an inverted vertex range";
        EXPECT_EQ(range.second, ranges[static_cast<std::size_t>(pe + 1)].first)
            << "vertex ownership has a gap (or an overlap) between PE " << pe << " and PE " << pe + 1;
    }
}

void TestPermutationOfEdgeBalancedGraph(
    const GraphRepresentation representation, const std::string& redistribution, const bool strict) {
    MPI_Comm comm = MPI_COMM_WORLD;
    PEID     rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const SInt n = 1024; // N=10

    PGeneratorConfig plain_config    = EdgeBalancedConfig(redistribution);
    PGeneratorConfig permuted_config = EdgeBalancedConfig(redistribution);
    permuted_config.permute          = true;

    const Graph plain    = GenerateInMemory(plain_config, representation, comm);
    const Graph permuted = GenerateInMemory(permuted_config, representation, comm);

    // The permuted graph must still be the same graph: undo the permutation and compare the global edge sets.
    // This is what catches edges being dropped (or duplicated, or reattached to the wrong vertex) while a
    // vertex that arrives in several partial pieces is reassembled.
    auto    permutator     = random_permutation::FeistelPseudoRandomPermutation::buildPermutation(n - 1, 0);
    EdgeVec local_permuted = LocalEdges(permuted);
    for (auto& [src, dst]: local_permuted) {
        src = permutator.finv(src);
        dst = permutator.finv(dst);
    }
    const EdgeVec gathered_plain    = GatherEdges(LocalEdges(plain), comm);
    const EdgeVec gathered_permuted = GatherEdges(local_permuted, comm);
    if (rank == 0) {
        EXPECT_EQ(gathered_plain, gathered_permuted);
    }

    // ... and it must still be balanced. Without re-establishing the balance after permuting, a PE's edge count
    // is whatever the permutation happens to hand it.
    const SInt local_edges = permuted.NumberOfLocalEdges();
    SInt       max_edges = local_edges, min_edges = local_edges, total_edges = local_edges;
    MPI_Allreduce(MPI_IN_PLACE, &max_edges, 1, KAGEN_MPI_SINT, MPI_MAX, comm);
    MPI_Allreduce(MPI_IN_PLACE, &min_edges, 1, KAGEN_MPI_SINT, MPI_MIN, comm);
    MPI_Allreduce(MPI_IN_PLACE, &total_edges, 1, KAGEN_MPI_SINT, MPI_SUM, comm);

    if (strict) {
        EXPECT_LE(max_edges - min_edges, 1u) << "permuted graph is not strictly edge-balanced";
    } else {
        // balance-edges keeps a vertex's whole adjacency on one PE, so a PE absorbs at most one vertex's
        // overflow past its fair share -- see ComputeBalancedEdgeDistribution().
        const SInt bucket_size = (total_edges + size - 1) / size;
        EXPECT_LE(max_edges, bucket_size + MaxDegree(plain, comm)) << "permuted graph is not edge-balanced";
    }

    ExpectGapFreeVertexRanges(permuted, n, comm);

    if (representation == GraphRepresentation::CSR) {
        const VertexRange row_range = permuted.PhysicalVertexRange();
        EXPECT_EQ(permuted.xadj.size(), row_range.second - row_range.first + 1)
            << "CSR row space does not match PhysicalVertexRange()";
        EXPECT_EQ(permuted.xadj.back(), permuted.adjncy.size());
    }

    // Split metadata must describe the permuted graph, not the one it was built from.
    bool any_split = permuted.has_split_vertices;
    MPI_Allreduce(MPI_IN_PLACE, &any_split, 1, MPI_C_BOOL, MPI_LOR, comm);
    if (strict && size > 1) {
        // Keeps the split-specific assertions below from passing vacuously: exact edge balance across more than
        // one PE puts a boundary inside some vertex's adjacency for this instance.
        EXPECT_TRUE(any_split) << "expected the strictly edge-balanced permuted graph to have a split vertex";
    }
    if (!any_split) {
        EXPECT_FALSE(permuted.left_partial_vertex.has_value());
        EXPECT_FALSE(permuted.right_partial_vertex.has_value());
        EXPECT_EQ(permuted.PhysicalVertexRange(), permuted.vertex_range);
    }
    if (permuted.left_partial_vertex) {
        EXPECT_EQ(permuted.left_partial_vertex->vertex, permuted.vertex_range.first - 1);
        EXPECT_TRUE(permuted.has_split_vertices);
    }
    if (permuted.right_partial_vertex) {
        EXPECT_EQ(permuted.right_partial_vertex->vertex, permuted.vertex_range.second - 1);
        EXPECT_TRUE(permuted.has_split_vertices);
    }
}
} // namespace

TEST(GraphPermutation, permuted_graph_stays_strictly_edge_balanced_edgelist) {
    TestPermutationOfEdgeBalancedGraph(GraphRepresentation::EDGE_LIST, "balance-edges-strict", true);
}

TEST(GraphPermutation, permuted_graph_stays_strictly_edge_balanced_csr) {
    TestPermutationOfEdgeBalancedGraph(GraphRepresentation::CSR, "balance-edges-strict", true);
}

TEST(GraphPermutation, permuted_graph_stays_edge_balanced_edgelist) {
    TestPermutationOfEdgeBalancedGraph(GraphRepresentation::EDGE_LIST, "balance-edges", false);
}

TEST(GraphPermutation, permuted_graph_stays_edge_balanced_csr) {
    TestPermutationOfEdgeBalancedGraph(GraphRepresentation::CSR, "balance-edges", false);
}
