#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/hyperbolic/hyperbolic.h"

#include <gtest/gtest.h>

#include "tests/gather.h"
#include "tests/geometric/utils.h"

using namespace kagen;

namespace {
void validate_graph(const Graph& local_graph, const PGeneratorConfig& config) {
    std::cout << "Range: " << local_graph.vertex_range.first << ".." << local_graph.vertex_range.second << std::endl;
    Graph global_graph = kagen::testing::GatherEdgeLists(local_graph);

    PEID rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        EXPECT_EQ(global_graph.NumberOfLocalVertices(), global_graph.coordinates.first.size());

        // Creating the correct edge list as a test instance
        std::vector<std::pair<SInt, SInt>> expected_edges =
            kagen::testing::CreateExpectedHyperbolicEdges<HPFloat>(config, global_graph);

        // Sorting both lists before comparing them
        std::sort(global_graph.edges.begin(), global_graph.edges.end());
        std::sort(expected_edges.begin(), expected_edges.end());

        EXPECT_EQ(global_graph.edges.size(), expected_edges.size());
        EXPECT_EQ(global_graph.edges, expected_edges);

        SInt missing_edges = 0;
        SInt extra_edges   = 0;

        std::size_t i = 0;
        std::size_t j = 0;
        while (i < global_graph.edges.size() && j < expected_edges.size()) {
            if (global_graph.edges[i] == expected_edges[j]) {
                ++i;
                ++j;
            } else if (global_graph.edges[i] < expected_edges[j]) {
                extra_edges++;
                ++i;
            } else if (global_graph.edges[i] > expected_edges[j]) {
                missing_edges++;
                ++j;
            }
        }
        missing_edges += expected_edges.size() - j;
        extra_edges += global_graph.edges.size() - i;

        std::cout << "No of vertices: " << global_graph.NumberOfLocalVertices() << std::endl;
        std::cout << "No of expected edges: " << expected_edges.size() << std::endl;
        std::cout << "No of actual edges: " << global_graph.edges.size() << std::endl;
        std::cout << "MISSING: " << missing_edges << ", EXTRA: " << extra_edges << std::endl;

        bool warned = false;
        for (std::size_t i = 0; !warned && i < global_graph.edges.size(); ++i) {
            auto [u, v] = global_graph.edges[i];

            if (i >= expected_edges.size()) {
                std::cout << "Warning: edge missing from expected: " << u << " -> " << v << std::endl;
                warned = true;
            } else if (expected_edges[i] != global_graph.edges[i]) {
                std::cout << "Warning: expected at pos " << i << ": " << expected_edges[i].first << " -> "
                          << expected_edges[i].second << ", got " << u << " -> " << v << std::endl;

                const auto    dist  = kagen::testing::ComputeHyperbolicDistance<HPFloat>(global_graph, u, v);
                const HPFloat alpha = (config.plexp - 1.0) / 2.0;
                const HPFloat r =
                    PGGeometry<HPFloat>::GetTargetRadius(config.n, config.n * config.avg_degree / 2, alpha);
                std::cout << "Distance: " << dist << ", vs radius: " << r << std::endl;

                warned = true;
            }
        }
    }
}

void test_configuration(const SInt n, const double avg_degree, const double plexp) {
    PGeneratorConfig config;
    config.n           = n;
    config.avg_degree  = avg_degree;
    config.plexp       = plexp;
    config.hp_floats   = true;
    config.coordinates = true;
    config.quiet       = false;

    HyperbolicFactory factory;
    PEID              size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    config         = factory.NormalizeParameters(config, rank, size, true);
    auto generator = factory.Create(config, rank, size);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    generator->Finalize(MPI_COMM_WORLD);
    const auto local_graph = generator->Take();

    validate_graph(local_graph, config);
}
} // namespace

TEST(RhgTest, generates_graph_on_np_PE_n16_d4_g3) {
    test_configuration(16, 4.0, 3.0);
}

TEST(RhgTest, generates_graph_on_np_PE_n512_d4_g3) {
    test_configuration(512, 4.0, 3.0);
}

TEST(RhgTest, generates_graph_on_np_PE_n1024_d16_g3) {
    test_configuration(1024, 16.0, 3.0);
}
