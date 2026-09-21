#include "kagen/kagen.h"
#include "kagen/tools/postprocessor.h"

#include <gtest/gtest.h>
#include <mpi.h>

#include "../gather.h"
#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

using namespace kagen;

using GeneratorFunc    = std::function<Graph(KaGen&, SInt, SInt)>;
using RedistributeFunc = std::function<VertexRange(Edgelist&, Edgelist&, SInt, bool, MPI_Comm)>;

// Wraps RedistributeEdgesTrueBalance's richer EdgeBalancedDistribution return type down to a plain VertexRange
// (its vertex_range, the complete gap-free partition -- see EdgeBalancedDistribution) so it fits
// RedistributeFunc and can be used interchangeably with the other redistribution functions in tests that only
// care about the final edge multiset, not about exact ownership of boundary/split vertices.
static VertexRange RedistributeEdgesTrueBalanceWrapped(
    Edgelist& source, Edgelist& destination, SInt n, bool remap_round_robin, MPI_Comm comm) {
    return RedistributeEdgesTrueBalance(source, destination, n, remap_round_robin, comm).vertex_range;
}

static Graph MakeEdgeListGraph(const Edgelist& edges) {
    Graph g;
    g.representation = GraphRepresentation::EDGE_LIST;
    g.edges          = edges;
    return g;
}

static Edgelist GatherAllEdges(const Edgelist& local_edges) {
    return kagen::testing::GatherEdgeLists(MakeEdgeListGraph(local_edges)).edges;
}

// Build a star graph (vertex 0 connected to all others, both directions) with all edges on PE 0.
static Edgelist BuildStarOnPE0(SInt n) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Edgelist edges;
    if (rank == 0) {
        for (SInt v = 1; v < n; ++v) {
            edges.emplace_back(0, v);
            edges.emplace_back(v, 0);
        }
    }
    return edges;
}

// Build a graph with one dominant hub vertex (degree = hub_degree, several times any single PE's fair
// share) plus a uniform overlay of edges among the remaining vertices. The hub's edges are scattered
// evenly across all PEs (as they would be for a real generator like RMAT before redistribution, where
// a hub vertex's edges are not pre-located on a single source PE).
static Edgelist BuildHubWithUniformOverlay(SInt n, SInt hub_degree, SInt overlay_edges_per_pe) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Edgelist edges;

    const SInt per_pe          = hub_degree / static_cast<SInt>(size);
    const SInt remainder       = hub_degree % static_cast<SInt>(size);
    const SInt local_hub_edges = per_pe + (static_cast<SInt>(rank) < remainder ? 1 : 0);
    const SInt local_hub_start = static_cast<SInt>(rank) * per_pe + std::min<SInt>(rank, remainder);
    for (SInt i = 0; i < local_hub_edges; ++i) {
        const SInt target = 1 + (local_hub_start + i) % (n - 1);
        edges.emplace_back(0, target);
    }

    for (SInt i = 0; i < overlay_edges_per_pe; ++i) {
        const SInt u = 1 + (static_cast<SInt>(rank) * 37 + i * 7) % (n - 1);
        const SInt v = 1 + (static_cast<SInt>(rank) * 53 + i * 11 + 1) % (n - 1);
        if (u != v) {
            edges.emplace_back(u, v);
        }
    }

    return edges;
}

// Build a graph where *every* PE independently contributes an identical copy of the same hub_degree hub edges
// (0, 1), ..., (0, hub_degree) -- modeling how a generator like RMAT can independently sample the same edge on
// different PEs (e.g. via symmetrization of two independently-sampled directed pairs (u, v) and (v, u); see
// discussion in DedupsHubEdgesDuplicatedAcrossSourcePEs) -- plus a uniform overlay so non-hub vertices have
// plenty of legitimate edges too.
static Edgelist BuildHubWithDuplicatesAcrossPEs(SInt n, SInt hub_degree, SInt overlay_edges_per_pe) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Edgelist edges;
    for (SInt target = 1; target <= hub_degree; ++target) {
        edges.emplace_back(0, target);
    }

    const SInt overlay_base = hub_degree + 2;
    for (SInt i = 0; i < overlay_edges_per_pe; ++i) {
        const SInt u = overlay_base + (static_cast<SInt>(rank) * 37 + i * 7) % (n - overlay_base);
        const SInt v = overlay_base + (static_cast<SInt>(rank) * 53 + i * 11 + 1) % (n - overlay_base);
        if (u != v) {
            edges.emplace_back(u, v);
        }
    }

    return edges;
}

// ---- Parameterized tests over generator, redistribution function, and round-robin remapping ---

using GeneratorParam    = std::tuple<std::string, GeneratorFunc>;
using RedistributeParam = std::tuple<std::string, RedistributeFunc>;
using FixtureParam      = std::tuple<GeneratorParam, RedistributeParam, bool>;

struct RedistributeEdgesFixture : public ::testing::TestWithParam<FixtureParam> {};

INSTANTIATE_TEST_SUITE_P(
    RedistributeEdgesTests, RedistributeEdgesFixture,
    ::testing::Combine(
        ::testing::Values(
            std::make_tuple("GNM", GeneratorFunc([](KaGen& gen, SInt n, SInt m) {
                                return gen.GenerateUndirectedGNM(n, m);
                            })),
            std::make_tuple("RMAT", GeneratorFunc([](KaGen& gen, SInt n, SInt m) {
                                return gen.GenerateRMAT(n, m, 0.56, 0.19, 0.19);
                            })),
            std::make_tuple("RGG2D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) {
                                return gen.GenerateRGG2D_NM(n, m);
                            })),
            std::make_tuple("RHG", GeneratorFunc([](KaGen& gen, SInt n, SInt m) {
                                return gen.GenerateRHG_NM(2.8, n, m);
                            })),
            std::make_tuple("Grid2D", GeneratorFunc([](KaGen& gen, SInt n, SInt m) {
                                return gen.GenerateGrid2D_N(n, 1.0);
                            }))),
        ::testing::Values(
            std::make_tuple("EdgeBalanced", RedistributeFunc(RedistributeEdgesBalanced)),
            std::make_tuple("VertexBalanced", RedistributeFunc(RedistributeEdges)),
            std::make_tuple("EdgeBalancedTrue", RedistributeFunc(RedistributeEdgesTrueBalanceWrapped))),
        ::testing::Values(true, false)),
    [](const ::testing::TestParamInfo<FixtureParam>& info) {
        return std::get<0>(std::get<0>(info.param)) + "_" + std::get<0>(std::get<1>(info.param))
               + (std::get<2>(info.param) ? "_remap" : "_noremap");
    });

TEST_P(RedistributeEdgesFixture, PreservesEdgeSet) {
    auto [gen_pair, redist_pair, remap_round_robin] = GetParam();
    auto generate                                   = std::get<1>(gen_pair);
    auto redistribute                               = std::get<1>(redist_pair);

    const SInt n = 1000;
    const SInt m = 4 * n;

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    Graph graph = generate(generator, n, m);

    Edgelist input        = graph.edges;
    SInt     num_vertices = graph.vertex_range.second;
    MPI_Allreduce(MPI_IN_PLACE, &num_vertices, 1, KAGEN_MPI_SINT, MPI_MAX, MPI_COMM_WORLD);

    Edgelist expected = [&](auto edges) {
        if (remap_round_robin) {
            RoundRobinRemapping(edges, num_vertices, MPI_COMM_WORLD);
        }
        edges = GatherAllEdges(edges);
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
        return edges;
    }(input);

    Edgelist redistributed_edges;
    redistribute(input, redistributed_edges, num_vertices, remap_round_robin, MPI_COMM_WORLD);

    Edgelist result = GatherAllEdges(redistributed_edges);
    std::sort(result.begin(), result.end());

    EXPECT_EQ(expected, result);
}

TEST_P(RedistributeEdgesFixture, OwnershipInvariant) {
    auto [gen_pair, redist_pair, remap_round_robin] = GetParam();
    auto generate                                   = std::get<1>(gen_pair);
    auto [redist_name, redistribute]                = redist_pair;

    if (redist_name == "EdgeBalancedTrue") {
        GTEST_SKIP() << "true edge-balancing intentionally splits a vertex's edges across PEs, so a single "
                        "PE's vertex range does not bound the tails of its edges -- see the dedicated "
                        "true-balance tests instead";
    }

    const SInt n = 1000;
    const SInt m = 4 * n;

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    Graph graph = generate(generator, n, m);

    Edgelist input        = graph.edges;
    SInt     num_vertices = graph.vertex_range.second;
    MPI_Allreduce(MPI_IN_PLACE, &num_vertices, 1, KAGEN_MPI_SINT, MPI_MAX, MPI_COMM_WORLD);

    Edgelist    redistributed_edges;
    VertexRange vr = redistribute(input, redistributed_edges, num_vertices, remap_round_robin, MPI_COMM_WORLD);

    for (const auto& edge: redistributed_edges) {
        EXPECT_GE(edge.first, vr.first);
        EXPECT_LT(edge.first, vr.second);
    }
}

TEST_P(RedistributeEdgesFixture, NoDuplicatesInOutput) {
    auto [gen_pair, redist_pair, remap_round_robin] = GetParam();
    auto generate                                   = std::get<1>(gen_pair);
    auto [redist_name, redistribute]                = redist_pair;

    const SInt n = 1000;
    const SInt m = 4 * n;

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    Graph graph = generate(generator, n, m);

    Edgelist input        = graph.edges;
    SInt     num_vertices = graph.vertex_range.second;
    MPI_Allreduce(MPI_IN_PLACE, &num_vertices, 1, KAGEN_MPI_SINT, MPI_MAX, MPI_COMM_WORLD);

    Edgelist redistributed_edges;
    redistribute(input, redistributed_edges, num_vertices, remap_round_robin, MPI_COMM_WORLD);

    EXPECT_TRUE(std::is_sorted(redistributed_edges.begin(), redistributed_edges.end()));
    auto dup = std::adjacent_find(redistributed_edges.begin(), redistributed_edges.end());
    EXPECT_EQ(dup, redistributed_edges.end());
}

// ---- Star graph and edge case tests -----------------------------------------

using SimpleFixtureParam = std::tuple<RedistributeParam, bool>;

struct RedistributeEdgesSimpleFixture : public ::testing::TestWithParam<SimpleFixtureParam> {};

INSTANTIATE_TEST_SUITE_P(
    RedistributeEdgesSimpleTests, RedistributeEdgesSimpleFixture,
    ::testing::Combine(
        ::testing::Values(
            std::make_tuple("EdgeBalanced", RedistributeFunc(RedistributeEdgesBalanced)),
            std::make_tuple("VertexBalanced", RedistributeFunc(RedistributeEdges)),
            std::make_tuple("EdgeBalancedTrue", RedistributeFunc(RedistributeEdgesTrueBalanceWrapped))),
        ::testing::Values(true, false)),
    [](const ::testing::TestParamInfo<SimpleFixtureParam>& info) {
        return std::get<0>(std::get<0>(info.param)) + (std::get<1>(info.param) ? "_remap" : "_noremap");
    });

TEST_P(RedistributeEdgesSimpleFixture, PreservesEdgeSet_Star) {
    auto [redist_pair, remap_round_robin] = GetParam();
    auto redistribute                     = std::get<1>(redist_pair);

    const SInt n     = 100;
    Edgelist   input = BuildStarOnPE0(n);

    Edgelist reference = input;
    if (remap_round_robin) {
        RoundRobinRemapping(reference, n, MPI_COMM_WORLD);
    }
    Edgelist expected = GatherAllEdges(reference);
    std::sort(expected.begin(), expected.end());
    expected.erase(std::unique(expected.begin(), expected.end()), expected.end());

    Edgelist redistributed_edges;
    redistribute(input, redistributed_edges, n, remap_round_robin, MPI_COMM_WORLD);

    Edgelist result = GatherAllEdges(redistributed_edges);
    std::sort(result.begin(), result.end());

    EXPECT_EQ(expected, result);
}

TEST_P(RedistributeEdgesSimpleFixture, OwnershipInvariant_Star) {
    auto [redist_pair, remap_round_robin] = GetParam();
    auto [redist_name, redistribute]      = redist_pair;

    if (redist_name == "EdgeBalancedTrue") {
        GTEST_SKIP() << "the star's hub vertex may legitimately be split across PEs under true edge-balancing";
    }

    const SInt  n     = 100;
    Edgelist    input = BuildStarOnPE0(n);
    Edgelist    redistributed_edges;
    VertexRange vr = redistribute(input, redistributed_edges, n, remap_round_robin, MPI_COMM_WORLD);

    for (const auto& edge: redistributed_edges) {
        EXPECT_GE(edge.first, vr.first);
        EXPECT_LT(edge.first, vr.second);
    }
}

TEST_P(RedistributeEdgesSimpleFixture, EmptyInput) {
    auto [redist_pair, remap_round_robin] = GetParam();
    auto redistribute                     = std::get<1>(redist_pair);

    const SInt n = 100;
    Edgelist   input;
    Edgelist   redistributed_edges;

    redistribute(input, redistributed_edges, n, remap_round_robin, MPI_COMM_WORLD);

    EXPECT_TRUE(redistributed_edges.empty());
}

TEST_P(RedistributeEdgesSimpleFixture, SingleEdge) {
    auto [redist_pair, remap_round_robin] = GetParam();
    auto redistribute                     = std::get<1>(redist_pair);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const SInt n = 10;
    Edgelist   input;
    if (rank == 0) {
        input.emplace_back(2, 7);
    }

    Edgelist reference = input;
    if (remap_round_robin) {
        RoundRobinRemapping(reference, n, MPI_COMM_WORLD);
    }
    Edgelist expected = GatherAllEdges(reference);
    std::sort(expected.begin(), expected.end());

    Edgelist    redistributed_edges;
    VertexRange vr = redistribute(input, redistributed_edges, n, remap_round_robin, MPI_COMM_WORLD);

    Edgelist result = GatherAllEdges(redistributed_edges);
    std::sort(result.begin(), result.end());

    EXPECT_EQ(expected, result);

    for (const auto& edge: redistributed_edges) {
        EXPECT_GE(edge.first, vr.first);
        EXPECT_LT(edge.first, vr.second);
    }
}

// ---- Adversarial hub vertex test (regression test for the vertex-atomic breakpoint bug) --------

// Before the fix, a single vertex whose degree spans multiple bucket boundaries would be emitted as a
// breakpoint once per boundary crossed, making Distribution hand several *other*, unrelated PEs an
// empty (v, v) vertex range -- i.e. those PEs get spuriously starved to zero vertices/edges even though
// there are plenty of edges (and other vertices) to go around.
TEST(RedistributeEdgesHubStarvationTest, NoPEIsStarvedWithDominantHub) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const SInt n = 4000;
    // The hub degree comfortably overflows a single PE's fair bucket share (to exercise the "absorb the
    // overflow" path), but the uniform overlay is dense enough that every *other* PE still has plenty of
    // edges to legitimately own -- unlike an extreme hub that swamps the entire graph, this does not
    // require more than one PE to go empty, so a real regression here is distinguishable from the
    // inherent, documented limit of vertex-atomic balancing (a hub that alone exceeds several PEs' worth
    // of edges necessarily forces neighboring PEs to be under-filled; that is what
    // RedistributeEdgesTrueBalance is for).
    const SInt hub_degree           = 1000;
    const SInt overlay_edges_per_pe = 400;

    Edgelist input = BuildHubWithUniformOverlay(n, hub_degree, overlay_edges_per_pe);

    Edgelist    redistributed_edges;
    VertexRange vr =
        RedistributeEdgesBalanced(input, redistributed_edges, n, /*remap_round_robin=*/false, MPI_COMM_WORLD);

    EXPECT_LT(vr.first, vr.second) << "PE " << rank << " was starved to an empty vertex range";
}

// True edge-balancing exists precisely to handle a hub so dominant that vertex-atomic balancing cannot give
// every PE a fair share (see NoPEIsStarvedWithDominantHub's doc comment): here the hub comfortably exceeds a
// single PE's fair bucket share even after accounting for the whole rest of the graph, so splitting the hub's
// own edges across multiple PEs is unavoidable for exact balance.
TEST(RedistributeEdgesTrueBalanceHubTest, SplitsDominantHubForExactBalance) {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // n must comfortably exceed hub_degree: BuildHubWithUniformOverlay's hub targets are only distinct modulo
    // (n - 1), so hub_degree > n - 1 would wrap around and create spurious duplicate (and thus deduplicated)
    // hub edges, undermining the exact-balance assertion below.
    const SInt n                    = 10000;
    const SInt hub_degree           = 5000;
    const SInt overlay_edges_per_pe = 50;

    Edgelist input = BuildHubWithUniformOverlay(n, hub_degree, overlay_edges_per_pe);

    Edgelist                 destination;
    EdgeBalancedDistribution dist =
        RedistributeEdgesTrueBalance(input, destination, n, /*remap_round_robin=*/false, MPI_COMM_WORLD);

    // (a) even the hub-absorbing PE is not left with an outsized share: every PE's local edge count is within
    // 1 of every other PE's.
    int local_count = static_cast<int>(destination.size());
    int min_count, max_count;
    MPI_Allreduce(&local_count, &min_count, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_count, &max_count, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    EXPECT_LE(max_count - min_count, 1) << "local edge counts are not balanced to within 1 across PEs";

    // (b) the hub vertex's tail appears on more than one PE once there is more than one PE to split it across.
    int has_hub_edge =
        std::any_of(destination.begin(), destination.end(), [](const auto& edge) { return edge.first == 0; }) ? 1 : 0;
    int num_pes_with_hub_edge = 0;
    MPI_Allreduce(&has_hub_edge, &num_pes_with_hub_edge, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (size > 1) {
        EXPECT_GT(num_pes_with_hub_edge, 1) << "hub vertex's edges ended up entirely on a single PE";
        EXPECT_TRUE(dist.has_split_vertices);

        // (c) fully_owned_vertex_range is a strict subset of vertex_range on every PE that actually shares a
        // boundary vertex (left_partial_vertex/right_partial_vertex set), and every vertex it reports is
        // excluded from fully_owned_vertex_range -- i.e. the trim is neither missing nor over-eager. For
        // vertex_range itself (the complete partition), a shared vertex spanning more than 2 PEs is resolved to
        // exactly the lowest-rank PE holding any of its edges: on that PE it's credited into vertex_range (as
        // that PE's right_partial_vertex, shared with a higher-rank neighbor); on every other PE holding a share
        // of it, it was already ceded away and so must fall strictly below vertex_range.first (as that PE's
        // left_partial_vertex, shared with a lower-rank neighbor).
        for (const auto& split_vertex: {dist.left_partial_vertex, dist.right_partial_vertex}) {
            if (!split_vertex) {
                continue;
            }
            EXPECT_FALSE(
                split_vertex->vertex >= dist.fully_owned_vertex_range.first
                && split_vertex->vertex < dist.fully_owned_vertex_range.second)
                << "shared boundary vertex " << split_vertex->vertex << " must be excluded from "
                << "fully_owned_vertex_range";
            const bool credited_here =
                split_vertex->vertex >= dist.vertex_range.first && split_vertex->vertex < dist.vertex_range.second;
            const bool ceded_to_lower_rank = split_vertex->vertex < dist.vertex_range.first;
            EXPECT_TRUE(credited_here || ceded_to_lower_rank)
                << "shared boundary vertex " << split_vertex->vertex << " neither credited to this PE's "
                << "vertex_range nor ceded to a lower rank";
        }
        if (!dist.left_partial_vertex && !dist.right_partial_vertex) {
            EXPECT_EQ(dist.fully_owned_vertex_range, dist.vertex_range)
                << "with no shared boundary vertex on this PE, fully_owned_vertex_range should equal vertex_range";
        }
    }
}

// Regression test: a generator like RMAT can independently sample the same undirected edge twice (once as
// (u, v), once as (v, u)) on two *different* PEs, since it symmetrizes every directed sample unconditionally
// and samples independently per chunk. If a vertex involved in such a duplicate gets split across PEs, the two
// copies are not guaranteed to be routed to the same target PE by position alone, so the final per-PE
// SortAndRemoveDuplicates -- which only catches duplicates that happen to co-locate -- misses them. This models
// the extreme case where the duplication is total (every PE contributes an identical copy of the whole hub).
//
// This also regression-tests the balance side of the same bug: the bucket boundaries used to route edges are
// computed from *pre-deduplication* degree counts, so if duplicates aren't accounted for before dedup happens,
// the post-dedup per-PE edge counts drift away from their +/-1 target (holes left behind by removed
// duplicates). See Stage 2 in RedistributeEdgesTrueBalance.
TEST(RedistributeEdgesTrueBalanceHubTest, DedupsHubEdgesDuplicatedAcrossSourcePEs) {
    const SInt n                    = 4000;
    const SInt hub_degree           = 1000;
    const SInt overlay_edges_per_pe = 50;

    Edgelist input = BuildHubWithDuplicatesAcrossPEs(n, hub_degree, overlay_edges_per_pe);

    Edgelist destination;
    RedistributeEdgesTrueBalance(input, destination, n, /*remap_round_robin=*/false, MPI_COMM_WORLD);

    Edgelist result = GatherAllEdges(destination);
    std::sort(result.begin(), result.end());
    auto dup = std::adjacent_find(result.begin(), result.end());
    EXPECT_EQ(dup, result.end()) << "duplicate edge survived redistribution";

    const auto hub_edges_in_result =
        std::count_if(result.begin(), result.end(), [](const auto& edge) { return edge.first == 0; });
    EXPECT_EQ(static_cast<SInt>(hub_edges_in_result), hub_degree)
        << "expected exactly hub_degree unique hub edges; duplicates contributed by different source PEs should "
           "have been removed";

    int local_count = static_cast<int>(destination.size());
    int min_count, max_count;
    MPI_Allreduce(&local_count, &min_count, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_count, &max_count, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    EXPECT_LE(max_count - min_count, 1)
        << "local edge counts are not balanced to within 1 across PEs once cross-PE duplicates are removed";
}

// Regression test for the "isolated vertices silently dropped" bug: vertex_range must be a complete, gap-free
// partition of [0, n) -- covering isolated (degree-0) vertices too, and resolving a shared boundary vertex to
// exactly one PE -- not just the vertices that happen to have local edges. Before the fix, a graph with
// mostly-isolated vertices (e.g. 16 vertices, 1 edge) would have its global vertex count silently collapse to
// just the handful of vertices touched by an edge.
TEST(RedistributeEdgesTrueBalanceRangeTest, VertexRangeCoversAllVerticesIncludingIsolated) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const SInt n = 16;
    Edgelist   input;
    if (rank == 0) {
        input.emplace_back(2, 7);
    }

    Edgelist                 destination;
    EdgeBalancedDistribution dist =
        RedistributeEdgesTrueBalance(input, destination, n, /*remap_round_robin=*/false, MPI_COMM_WORLD);

    // fully_owned_vertex_range is always a subset of vertex_range, and this graph is far too sparse to trigger
    // an actual split, so the two should coincide exactly here.
    EXPECT_EQ(dist.fully_owned_vertex_range, dist.vertex_range)
        << "no split should have occurred for a single edge among mostly-isolated vertices";

    std::vector<SInt> starts(static_cast<std::size_t>(size));
    std::vector<SInt> ends(static_cast<std::size_t>(size));
    MPI_Allgather(&dist.vertex_range.first, 1, KAGEN_MPI_SINT, starts.data(), 1, KAGEN_MPI_SINT, MPI_COMM_WORLD);
    MPI_Allgather(&dist.vertex_range.second, 1, KAGEN_MPI_SINT, ends.data(), 1, KAGEN_MPI_SINT, MPI_COMM_WORLD);

    EXPECT_EQ(starts.front(), 0) << "the partition of [0, n) must start at 0";
    EXPECT_EQ(ends.back(), n) << "the partition of [0, n) must end at n, including any trailing isolated vertices";
    for (int pe = 0; pe < size; ++pe) {
        EXPECT_LE(starts[static_cast<std::size_t>(pe)], ends[static_cast<std::size_t>(pe)])
            << "PE " << pe << " has an invalid (start > end) range";
    }
    for (int pe = 0; pe + 1 < size; ++pe) {
        EXPECT_EQ(ends[static_cast<std::size_t>(pe)], starts[static_cast<std::size_t>(pe + 1)])
            << "gap or overlap between PE " << pe << " and PE " << pe + 1;
    }
}

// General regression net: for every redistribution mode (including true edge-balancing), the returned vertex
// range's sizes must sum to exactly the graph's vertex count across all PEs -- i.e. no PE's range overlaps
// another's and no vertex (isolated or not) is left uncovered.
TEST_P(RedistributeEdgesFixture, VertexRangeSizesSumToN) {
    auto [gen_pair, redist_pair, remap_round_robin] = GetParam();
    auto generate                                   = std::get<1>(gen_pair);
    auto redistribute                               = std::get<1>(redist_pair);

    const SInt n = 1000;
    const SInt m = 4 * n;

    kagen::KaGen generator(MPI_COMM_WORLD);
    generator.UseEdgeListRepresentation();
    Graph graph = generate(generator, n, m);

    Edgelist input        = graph.edges;
    SInt     num_vertices = graph.vertex_range.second;
    MPI_Allreduce(MPI_IN_PLACE, &num_vertices, 1, KAGEN_MPI_SINT, MPI_MAX, MPI_COMM_WORLD);

    Edgelist    redistributed_edges;
    VertexRange vr = redistribute(input, redistributed_edges, num_vertices, remap_round_robin, MPI_COMM_WORLD);

    SInt local_n = vr.second - vr.first;
    SInt total_n = 0;
    MPI_Allreduce(&local_n, &total_n, 1, KAGEN_MPI_SINT, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(total_n, num_vertices);
}
