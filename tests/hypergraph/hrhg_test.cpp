#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/hyper/h_hyperbolic/hyper_hyperbolic.h"

#include <gtest/gtest.h>
#include <mpi.h>

#include "hypergraph/utils.h"
#include <algorithm>
#include <cstddef>
#include <unordered_set>
#include <utility>

using namespace kagen;
using namespace kagen::testing;

namespace {

PGeneratorConfig DefaultConfig() {
    PGeneratorConfig config;

    config.generator = GeneratorType::H_RHG;

    config.n = 5000;
    config.m = 1000;

    config.avg_degree = 16.0;
    config.plexp      = 2.8;

    config.seed  = 42;
    config.quiet = true;
    config.debug = false;

    config.random_radius = false;
    config.r             = 4.0;

    config.partial_cell_mode = PartialCellMode::GenerateAndCheck;

    return config;
}

bool IsRoot() {
    PEID rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    return rank == 0;
}

Graph GenerateHyperbolic(PGeneratorConfig config) {
    Hyper_HyperbolicFactory factory;

    PEID rank;
    PEID size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    config = factory.NormalizeParameters(config, rank, size, false);

    auto generator = factory.Create(config, rank, size);

    generator->Generate(GraphRepresentation::CSR);

    generator->Finalize(MPI_COMM_WORLD);

    return generator->Take();
}

ExpandedHypergraph GenerateAndGather(PGeneratorConfig config) {
    const Graph local_graph = GenerateHyperbolic(std::move(config));

    return GatherHypergraph(local_graph);
}

ExpandedHypergraph Canonicalize(ExpandedHypergraph hypergraph) {
    for (auto& edge: hypergraph) {
        std::sort(edge.begin(), edge.end());
    }

    std::sort(hypergraph.begin(), hypergraph.end());

    return hypergraph;
}

void ExpectValidHypergraph(const ExpandedHypergraph& hypergraph, const SInt num_vertices) {
    for (std::size_t e = 0; e < hypergraph.size(); ++e) {
        const auto& edge = hypergraph[e];

        std::unordered_set<SInt> seen;

        for (const SInt pin: edge) {
            EXPECT_GE(pin, SInt{0}) << "negative pin in hyperedge " << e;

            EXPECT_LT(pin, num_vertices) << "invalid pin in hyperedge " << e;

            EXPECT_TRUE(seen.insert(pin).second) << "duplicate pin " << pin << " in hyperedge " << e;
        }
    }
}

std::size_t CountPins(const ExpandedHypergraph& hypergraph) {
    std::size_t pins = 0;

    for (const auto& edge: hypergraph) {
        pins += edge.size();
    }

    return pins;
}

} // namespace

TEST(HHyperbolic, GeneratesValidHypergraphExact) {
    auto config = DefaultConfig();

    config.random_radius = false;
    config.r             = 4.0;

    config.partial_cell_mode = PartialCellMode::GenerateAndCheck;

    const auto hypergraph = GenerateAndGather(config);

    if (!IsRoot()) {
        return;
    }

    ASSERT_EQ(hypergraph.size(), static_cast<std::size_t>(config.m));

    ExpectValidHypergraph(hypergraph, config.n);

    EXPECT_GT(CountPins(hypergraph), std::size_t{0});
}

TEST(HHyperbolic, GeneratesValidHypergraphCoverageRange) {
    auto config = DefaultConfig();

    config.partial_cell_mode = PartialCellMode::EstimateByCoverageRange;

    const auto hypergraph = GenerateAndGather(config);

    if (!IsRoot()) {
        return;
    }

    ASSERT_EQ(hypergraph.size(), static_cast<std::size_t>(config.m));

    ExpectValidHypergraph(hypergraph, config.n);

    EXPECT_GT(CountPins(hypergraph), std::size_t{0});
}

TEST(HHyperbolic, GeneratesValidHypergraphCoverageFloyd) {
    auto config = DefaultConfig();

    config.partial_cell_mode = PartialCellMode::EstimateByCoverageFloyd;

    const auto hypergraph = GenerateAndGather(config);

    if (!IsRoot()) {
        return;
    }

    ASSERT_EQ(hypergraph.size(), static_cast<std::size_t>(config.m));

    ExpectValidHypergraph(hypergraph, config.n);

    EXPECT_GT(CountPins(hypergraph), std::size_t{0});
}

TEST(HHyperbolic, SupportsRandomRadius) {
    auto config = DefaultConfig();

    config.random_radius = true;

    config.min_hyperedge_radius = 0.25;
    config.max_hyperedge_radius = 10.0;

    config.hyperedge_radius_exponent = 30.0;

    const auto hypergraph = GenerateAndGather(config);

    if (!IsRoot()) {
        return;
    }
    ASSERT_EQ(hypergraph.size(), static_cast<std::size_t>(config.m));

    ExpectValidHypergraph(hypergraph, config.n);
}

TEST(HHyperbolic, SupportsPinBudgetMode) {
    auto config = DefaultConfig();

    config.random_radius = true;

    config.size_dist_lower_bound = 2;
    config.size_dist_upper_bound = 16;

    //
    // Mean expected hyperedge size:
    //
    //     8000 / 1000 = 8
    //
    config.size_dist_pin_budget = 8000;

    const auto hypergraph = GenerateAndGather(config);

    if (!IsRoot()) {
        return;
    }

    ASSERT_EQ(hypergraph.size(), static_cast<std::size_t>(config.m));

    ExpectValidHypergraph(hypergraph, config.n);
}

TEST(HHyperbolic, SameSeedProducesSameHypergraph) {
    auto config = DefaultConfig();

    config.seed = 42;

    const auto first = GenerateAndGather(config);

    const auto second = GenerateAndGather(config);

    if (!IsRoot()) {
        return;
    }
    ASSERT_GT(CountPins(first), std::size_t{0});
    ASSERT_GT(CountPins(second), std::size_t{0});
    EXPECT_EQ(Canonicalize(first), Canonicalize(second));
}

TEST(HHyperbolic, DifferentSeedsProduceDifferentHypergraphs) {
    auto first_config = DefaultConfig();

    auto second_config = DefaultConfig();

    first_config.seed  = 1;
    second_config.seed = 42;

    first_config.r  = 4.0;
    second_config.r = 4.0;

    const auto first = GenerateAndGather(first_config);

    const auto second = GenerateAndGather(second_config);

    if (!IsRoot()) {
        return;
    }
    ASSERT_GT(CountPins(first), std::size_t{0});
    ASSERT_GT(CountPins(second), std::size_t{0});
    EXPECT_NE(Canonicalize(first), Canonicalize(second));
}