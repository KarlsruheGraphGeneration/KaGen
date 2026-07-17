#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/in_memory_facade.h"
#include "kagen/kagen.h"

#include <gtest/gtest.h>

using namespace kagen;

// CheckSplitVertexCompatibility is pure decision logic over an already-reduced any_split flag and a
// PGeneratorConfig, so it's tested directly here without driving any MPI collectives or actual generation.

TEST(CompatibilityGuardTest, NoSplitAlwaysSucceeds) {
    PGeneratorConfig config;
    config.output_graph.formats = {FileFormat::METIS};
    config.validate_simple_graph = true;
    config.edge_weights.generator_type = EdgeWeightGeneratorType::HASHING_BASED;
    config.statistics_level             = StatisticsLevel::ADVANCED;

    EXPECT_NO_THROW(CheckSplitVertexCompatibility(false, GraphRepresentation::CSR, config));
    EXPECT_NO_THROW(CheckSplitVertexCompatibility(false, GraphRepresentation::EDGE_LIST, config));
}

TEST(CompatibilityGuardTest, SplitWithEdgeListFormatSucceeds) {
    PGeneratorConfig config;
    config.output_graph.formats = {FileFormat::EDGE_LIST, FileFormat::BINARY_EDGE_LIST};
    config.statistics_level     = StatisticsLevel::NONE;

    EXPECT_NO_THROW(CheckSplitVertexCompatibility(true, GraphRepresentation::EDGE_LIST, config));
}

TEST(CompatibilityGuardTest, SplitWithCSRRepresentationSucceeds) {
    // A split-vertex CSR graph is only ever produced by the direct strict edge-balanced file read
    // (ParhipReader::ReadStrictEdgeRange), which builds valid xadj/adjncy with partial leading/trailing rows
    // (described by partial_vertices) -- it never goes through BuildCSRFromEdgeList. The paths that WOULD build
    // a corrupt CSR from a split edge list throw earlier (EdgeListOnlyGenerator::FinalizeCSR,
    // FileGraphGenerator::FinalizeCSR), so this guard no longer rejects CSR here.
    PGeneratorConfig config;
    config.output_graph.formats = {FileFormat::EDGE_LIST}; // otherwise-compatible output format
    config.statistics_level     = StatisticsLevel::NONE;

    EXPECT_NO_THROW(CheckSplitVertexCompatibility(true, GraphRepresentation::CSR, config));
}

struct AdjacencyGroupedFormatParam : public ::testing::TestWithParam<FileFormat> {};

INSTANTIATE_TEST_SUITE_P(
    CompatibilityGuardTests, AdjacencyGroupedFormatParam,
    ::testing::Values(
        FileFormat::METIS, FileFormat::HMETIS, FileFormat::HMETIS_DIRECTED, FileFormat::HMETIS_EP, FileFormat::DOT,
        FileFormat::DOT_DIRECTED, FileFormat::PARHIP, FileFormat::XTRAPULP, FileFormat::FREIGHT_NETL,
        FileFormat::FREIGHT_NETL_EP, FileFormat::NETD_ARE));

TEST_P(AdjacencyGroupedFormatParam, SplitWithAdjacencyGroupedFormatThrows) {
    PGeneratorConfig config;
    config.output_graph.formats = {GetParam()};
    config.statistics_level     = StatisticsLevel::NONE;

    EXPECT_THROW(CheckSplitVertexCompatibility(true, GraphRepresentation::EDGE_LIST, config), ConfigurationError);
}

TEST(CompatibilityGuardTest, SplitWithValidateSimpleGraphThrows) {
    PGeneratorConfig config;
    config.output_graph.formats  = {FileFormat::EDGE_LIST};
    config.validate_simple_graph = true;
    config.statistics_level      = StatisticsLevel::NONE;

    EXPECT_THROW(CheckSplitVertexCompatibility(true, GraphRepresentation::EDGE_LIST, config), ConfigurationError);
}

TEST(CompatibilityGuardTest, SplitWithHashingBasedEdgeWeightsThrows) {
    PGeneratorConfig config;
    config.output_graph.formats        = {FileFormat::EDGE_LIST};
    config.edge_weights.generator_type = EdgeWeightGeneratorType::HASHING_BASED;
    config.statistics_level            = StatisticsLevel::NONE;

    EXPECT_THROW(CheckSplitVertexCompatibility(true, GraphRepresentation::EDGE_LIST, config), ConfigurationError);
}

TEST(CompatibilityGuardTest, SplitWithUniformRandomEdgeWeightsThrows) {
    PGeneratorConfig config;
    config.output_graph.formats        = {FileFormat::EDGE_LIST};
    config.edge_weights.generator_type = EdgeWeightGeneratorType::UNIFORM_RANDOM;
    config.statistics_level            = StatisticsLevel::NONE;

    EXPECT_THROW(CheckSplitVertexCompatibility(true, GraphRepresentation::EDGE_LIST, config), ConfigurationError);
}

TEST(CompatibilityGuardTest, SplitWithDefaultOrVoidingEdgeWeightsSucceeds) {
    PGeneratorConfig config;
    config.output_graph.formats = {FileFormat::EDGE_LIST};
    config.statistics_level     = StatisticsLevel::NONE;

    config.edge_weights.generator_type = EdgeWeightGeneratorType::DEFAULT;
    EXPECT_NO_THROW(CheckSplitVertexCompatibility(true, GraphRepresentation::EDGE_LIST, config));

    config.edge_weights.generator_type = EdgeWeightGeneratorType::VOIDING;
    EXPECT_NO_THROW(CheckSplitVertexCompatibility(true, GraphRepresentation::EDGE_LIST, config));

    config.edge_weights.generator_type = EdgeWeightGeneratorType::EUCLIDEAN_DISTANCE;
    EXPECT_NO_THROW(CheckSplitVertexCompatibility(true, GraphRepresentation::EDGE_LIST, config));
}

TEST(CompatibilityGuardTest, SplitWithStatisticsThrows) {
    PGeneratorConfig config;
    config.output_graph.formats = {FileFormat::EDGE_LIST};
    config.quiet                = false;

    config.statistics_level = StatisticsLevel::BASIC;
    EXPECT_THROW(CheckSplitVertexCompatibility(true, GraphRepresentation::EDGE_LIST, config), ConfigurationError);

    config.statistics_level = StatisticsLevel::ADVANCED;
    EXPECT_THROW(CheckSplitVertexCompatibility(true, GraphRepresentation::EDGE_LIST, config), ConfigurationError);
}

TEST(CompatibilityGuardTest, SplitWithStatisticsButQuietSucceeds) {
    // --quiet suppresses statistics computation entirely, regardless of statistics_level (matches the actual
    // gating in GenerateInMemory: "if (!config.quiet) { ... if (statistics_level >= BASIC) ... }").
    PGeneratorConfig config;
    config.output_graph.formats = {FileFormat::EDGE_LIST};
    config.quiet                = true;
    config.statistics_level     = StatisticsLevel::ADVANCED;

    EXPECT_NO_THROW(CheckSplitVertexCompatibility(true, GraphRepresentation::EDGE_LIST, config));
}

TEST(CompatibilityGuardTest, SplitWithMultipleFormatsThrowsIfAnyIsIncompatible) {
    PGeneratorConfig config;
    config.output_graph.formats = {FileFormat::EDGE_LIST, FileFormat::METIS, FileFormat::BINARY_EDGE_LIST};
    config.statistics_level     = StatisticsLevel::NONE;

    EXPECT_THROW(CheckSplitVertexCompatibility(true, GraphRepresentation::EDGE_LIST, config), ConfigurationError);
}
