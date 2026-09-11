#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/hyper/h_erdos/cigam.h"

#include <gtest/gtest.h>
#include <mpi.h>

#include "hypergraph/utils.h"
#include <algorithm>
#include <cstddef>
#include <set>
#include <utility>
#include <vector>

using namespace kagen;

namespace {

Graph GenerateCIGAM(PGeneratorConfig config) {
    HyperCIGAMFactory factory;

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

bool IsRoot() {
    PEID rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank == 0;
}

PGeneratorConfig BaseConfig() {
    PGeneratorConfig config;

    config.n    = 128;
    config.k    = 8;
    config.seed = 1;

    config.cigam_mode        = CIGAMMode::EXACT;
    config.cigam_lambda      = 1.0;
    config.cigam_c           = {2.0};
    config.cigam_breakpoints = {1.0};
    config.cigam_sizes       = {3};

    // Edge-budget scaling keeps these unit tests small while exercising all
    // dominant-position blocks. The realized count remains random.
    config.edge_budget = 64;

    config.allow_duplicates = false;
    return config;
}

kagen::testing::ExpandedHypergraph Canonicalize(kagen::testing::ExpandedHypergraph hypergraph) {
    for (auto& edge: hypergraph) {
        std::sort(edge.begin(), edge.end());
    }
    std::sort(hypergraph.begin(), hypergraph.end());
    return hypergraph;
}

void ExpectValidPins(const kagen::testing::ExpandedHypergraph& hypergraph, const SInt n) {
    for (std::size_t edge_index = 0; edge_index < hypergraph.size(); ++edge_index) {
        const auto&    edge = hypergraph[edge_index];
        std::set<SInt> pins;

        for (const SInt pin: edge) {
            EXPECT_GE(pin, 0) << "in hyperedge " << edge_index;
            EXPECT_LT(pin, n) << "in hyperedge " << edge_index;
            EXPECT_TRUE(pins.insert(pin).second) << "duplicate pin " << pin << " in hyperedge " << edge_index;
        }
    }
}

void ExpectNoDuplicateHyperedges(kagen::testing::ExpandedHypergraph hypergraph) {
    hypergraph = Canonicalize(std::move(hypergraph));

    for (std::size_t i = 1; i < hypergraph.size(); ++i) {
        EXPECT_NE(hypergraph[i - 1], hypergraph[i]) << "duplicate hyperedge at canonical index " << i;
    }
}

} // namespace

TEST(CIGAMTest, ExactModeProducesValidFixedSizeHyperedges) {
    const auto  config            = BaseConfig();
    const Graph local_graph       = GenerateCIGAM(config);
    const auto  global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    ASSERT_FALSE(global_hypergraph.empty());
    for (const auto& edge: global_hypergraph) {
        EXPECT_EQ(edge.size(), std::size_t{3});
    }
    ExpectValidPins(global_hypergraph, config.n);
}

TEST(CIGAMTest, ExactModeSupportsConfiguredSizeRange) {
    auto config = BaseConfig();
    config.cigam_sizes.clear();
    config.size_dist_lower_bound = 2;
    config.size_dist_upper_bound = 5;
    config.edge_budget           = 32;

    const Graph local_graph       = GenerateCIGAM(config);
    const auto  global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    ASSERT_FALSE(global_hypergraph.empty());
    for (const auto& edge: global_hypergraph) {
        EXPECT_GE(edge.size(), std::size_t{2});
        EXPECT_LE(edge.size(), std::size_t{5});
    }
    ExpectValidPins(global_hypergraph, config.n);
}

TEST(CIGAMTest, ExactModeProducesNoDuplicateHyperedgesWhenDisabled) {
    auto config             = BaseConfig();
    config.allow_duplicates = false;

    const Graph local_graph       = GenerateCIGAM(config);
    auto        global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    ASSERT_FALSE(global_hypergraph.empty());
    ExpectNoDuplicateHyperedges(std::move(global_hypergraph));
}

TEST(CIGAMTest, SameSeedProducesSameExactHypergraph) {
    auto config = BaseConfig();
    config.seed = 42;

    const Graph first_local  = GenerateCIGAM(config);
    const auto  first_global = kagen::testing::GatherHypergraph(first_local);

    const Graph second_local  = GenerateCIGAM(config);
    const auto  second_global = kagen::testing::GatherHypergraph(second_local);

    if (!IsRoot()) {
        return;
    }

    EXPECT_EQ(Canonicalize(first_global), Canonicalize(second_global));
}

TEST(CIGAMTest, MultipleLayersPreserveBasicHypergraphConstraints) {
    auto config              = BaseConfig();
    config.cigam_c           = {1.5, 2.5};
    config.cigam_breakpoints = {0.25, 1.0};

    const Graph local_graph       = GenerateCIGAM(config);
    const auto  global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    ASSERT_FALSE(global_hypergraph.empty());
    for (const auto& edge: global_hypergraph) {
        EXPECT_EQ(edge.size(), std::size_t{3});
    }
    ExpectValidPins(global_hypergraph, config.n);
}

TEST(CIGAMTest, ApproximateModePreservesBasicHypergraphConstraints) {
    auto config        = BaseConfig();
    config.cigam_mode  = CIGAMMode::APPROX;
    config.cigam_sizes = {3, 5};

    const Graph local_graph       = GenerateCIGAM(config);
    const auto  global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    ASSERT_FALSE(global_hypergraph.empty());
    for (const auto& edge: global_hypergraph) {
        EXPECT_TRUE(edge.size() == std::size_t{3} || edge.size() == std::size_t{5});
    }
    ExpectValidPins(global_hypergraph, config.n);
}

TEST(CIGAMTest, PaperModePreservesBasicHypergraphConstraints) {
    auto config             = BaseConfig();
    config.n                = 32;
    config.cigam_mode       = CIGAMMode::PAPER;
    config.cigam_sizes      = {2};
    config.edge_budget      = 0;
    config.allow_duplicates = false;

    const Graph local_graph       = GenerateCIGAM(config);
    const auto  global_hypergraph = kagen::testing::GatherHypergraph(local_graph);

    if (!IsRoot()) {
        return;
    }

    ASSERT_FALSE(global_hypergraph.empty());
    for (const auto& edge: global_hypergraph) {
        EXPECT_EQ(edge.size(), std::size_t{2});
    }
    ExpectValidPins(global_hypergraph, config.n);
    ExpectNoDuplicateHyperedges(global_hypergraph);
}

TEST(CIGAMTest, RejectsNonPositiveLambda) {
    auto config         = BaseConfig();
    config.cigam_lambda = 0.0;
    EXPECT_THROW(GenerateCIGAM(config), ConfigurationError);
}

TEST(CIGAMTest, RejectsBreakpointsThatDoNotCoverAllLayers) {
    auto config              = BaseConfig();
    config.cigam_c           = {1.5, 2.5};
    config.cigam_breakpoints = {1.0};
    EXPECT_THROW(GenerateCIGAM(config), ConfigurationError);
}

TEST(CIGAMTest, RejectsInvalidExplicitHyperedgeSize) {
    auto config        = BaseConfig();
    config.cigam_sizes = {1};
    EXPECT_THROW(GenerateCIGAM(config), ConfigurationError);
}
