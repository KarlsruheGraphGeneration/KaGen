#pragma once

#include "kagen/context.h"
#include "kagen/definitions.h"

#include <gtest/gtest.h>
#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <unordered_set>
#include <vector>

namespace kagen::testing {

using ExpandedHypergraph = std::vector<std::vector<SInt>>;

inline ExpandedHypergraph ExpandLocalHypergraph(const Graph& graph) {
    const SInt num_hyperedges = graph.NumberOfLocalHyperedges();

    ExpandedHypergraph hyperedges(static_cast<std::size_t>(num_hyperedges));

    for (SInt e = 0; e < num_hyperedges; ++e) {
        auto& pins = hyperedges[static_cast<std::size_t>(e)];

        //
        // Explicit pins.
        //
        const SInt pin_begin = graph.hyperedge_offsets[e];

        const SInt pin_end = graph.hyperedge_offsets[e + 1];

        pins.insert(pins.end(), graph.hyperedge_pins.begin() + pin_begin, graph.hyperedge_pins.begin() + pin_end);

        //
        // Compressed pin ranges.
        //
        if (!graph.hyperedge_range_offsets.empty()) {
            const SInt range_begin = graph.hyperedge_range_offsets[e];

            const SInt range_end = graph.hyperedge_range_offsets[e + 1];

            for (SInt r = range_begin; r < range_end; ++r) {
                const auto& range = graph.hyperedge_ranges[r];

                for (SInt vertex = range.begin; vertex < range.end; ++vertex) {
                    pins.push_back(vertex);
                }
            }
        }

        std::sort(pins.begin(), pins.end());
    }

    return hyperedges;
}

struct HypergraphStats {
    SInt num_vertices   = 0;
    SInt num_hyperedges = 0;
    SInt num_pins       = 0;

    SInt min_size = 0;
    SInt max_size = 0;

    double mean_size     = 0.0;
    double variance_size = 0.0;
};

inline std::vector<SInt> HyperedgeSizes(const Graph& graph) {
    const auto hyperedges = ExpandLocalHypergraph(graph);

    std::vector<SInt> sizes;

    sizes.reserve(hyperedges.size());

    for (const auto& edge: hyperedges) {
        sizes.push_back(static_cast<SInt>(edge.size()));
    }

    return sizes;
}

inline HypergraphStats ComputeHypergraphStats(const Graph& graph, const SInt num_vertices) {
    HypergraphStats stats;

    stats.num_vertices = num_vertices;

    if (graph.hyperedge_offsets.empty()) {
        return stats;
    }

    stats.num_hyperedges = static_cast<SInt>(graph.hyperedge_offsets.size() - 1);

    stats.num_pins = graph.NumberOfLocalPins();

    const auto sizes = HyperedgeSizes(graph);

    if (sizes.empty()) {
        return stats;
    }

    const auto [min_it, max_it] = std::minmax_element(sizes.begin(), sizes.end());

    stats.min_size = *min_it;
    stats.max_size = *max_it;

    const double sum = std::accumulate(sizes.begin(), sizes.end(), 0.0);

    stats.mean_size = sum / static_cast<double>(sizes.size());

    double squared_error = 0.0;

    for (const SInt size: sizes) {
        const double diff = static_cast<double>(size) - stats.mean_size;

        squared_error += diff * diff;
    }

    stats.variance_size = squared_error / static_cast<double>(sizes.size());

    return stats;
}

inline void ValidateHypergraphStructure(const Graph& graph, const SInt num_vertices, const SInt expected_hyperedges) {
    ASSERT_FALSE(graph.hyperedge_offsets.empty());

    ASSERT_EQ(graph.hyperedge_offsets.front(), 0);

    ASSERT_EQ(graph.hyperedge_offsets.size(), static_cast<std::size_t>(expected_hyperedges + 1));

    EXPECT_EQ(graph.hyperedge_offsets.back(), graph.hyperedge_pins.size());

    for (std::size_t i = 1; i < graph.hyperedge_offsets.size(); ++i) {
        EXPECT_LE(graph.hyperedge_offsets[i - 1], graph.hyperedge_offsets[i]);
    }

    if (!graph.hyperedge_range_offsets.empty()) {
        ASSERT_EQ(graph.hyperedge_range_offsets.size(), static_cast<std::size_t>(expected_hyperedges + 1));

        ASSERT_EQ(graph.hyperedge_range_offsets.front(), 0);

        EXPECT_EQ(graph.hyperedge_range_offsets.back(), graph.hyperedge_ranges.size());

        for (std::size_t i = 1; i < graph.hyperedge_range_offsets.size(); ++i) {
            EXPECT_LE(graph.hyperedge_range_offsets[i - 1], graph.hyperedge_range_offsets[i]);
        }
    }

    const auto hyperedges = ExpandLocalHypergraph(graph);

    for (std::size_t e = 0; e < hyperedges.size(); ++e) {
        for (const SInt pin: hyperedges[e]) {
            EXPECT_LT(pin, num_vertices) << "invalid pin in hyperedge " << e;
        }
    }
}

inline void ValidateNoDuplicatePins(const Graph& graph) {
    const auto hyperedges = ExpandLocalHypergraph(graph);

    for (std::size_t e = 0; e < hyperedges.size(); ++e) {
        const auto& edge = hyperedges[e];

        for (std::size_t i = 1; i < edge.size(); ++i) {
            EXPECT_NE(edge[i - 1], edge[i]) << "duplicate vertex " << edge[i] << " in hyperedge " << e;
        }
    }
}

inline void ValidateSizeBounds(const Graph& graph, const SInt lower, const SInt upper) {
    const auto sizes = HyperedgeSizes(graph);

    for (const SInt size: sizes) {
        EXPECT_GE(size, lower);

        if (upper > 0) {
            EXPECT_LE(size, upper);
        }
    }
}

inline void ExpectNearRelative(const double actual, const double expected, const double relative_tolerance) {
    const double tolerance = relative_tolerance * std::max(1.0, std::abs(expected));

    EXPECT_NEAR(actual, expected, tolerance);
}

inline ExpandedHypergraph GatherHypergraph(const Graph& local_graph) {
    const ExpandedHypergraph local_hyperedges = ExpandLocalHypergraph(local_graph);

    //
    // Convert the local hypergraph into:
    //
    //   local_sizes = [|e_0|, |e_1|, ...]
    //   local_pins  = [pins(e_0), pins(e_1), ...]
    //
    std::vector<SInt> local_sizes;
    std::vector<SInt> local_pins;

    local_sizes.reserve(local_hyperedges.size());

    std::size_t local_pin_count = 0;

    for (const auto& edge: local_hyperedges) {
        local_pin_count += edge.size();
    }

    local_pins.reserve(local_pin_count);

    for (const auto& edge: local_hyperedges) {
        local_sizes.push_back(static_cast<SInt>(edge.size()));

        local_pins.insert(local_pins.end(), edge.begin(), edge.end());
    }

    //
    // MPI_Allgatherv uses int for element counts and displacements.
    // This is fine for the small instances used by the test suite,
    // but guard the conversion explicitly.
    //
    if (local_sizes.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        ADD_FAILURE() << "too many local hyperedges "
                         "for test gathering";

        return {};
    }

    if (local_pins.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        ADD_FAILURE() << "too many local pins "
                         "for test gathering";

        return {};
    }

    const int num_local_hyperedges = static_cast<int>(local_sizes.size());

    const int num_local_pins = static_cast<int>(local_pins.size());

    PEID num_pes;

    MPI_Comm_size(MPI_COMM_WORLD, &num_pes);

    //
    // First gather how many hyperedges and pins each PE owns.
    //
    std::vector<int> hyperedge_counts(num_pes);

    std::vector<int> pin_counts(num_pes);

    MPI_Allgather(&num_local_hyperedges, 1, MPI_INT, hyperedge_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    MPI_Allgather(&num_local_pins, 1, MPI_INT, pin_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    //
    // Calculate receive offsets for MPI_Allgatherv.
    //
    std::vector<int> hyperedge_displacements(num_pes);

    std::vector<int> pin_displacements(num_pes);

    std::exclusive_scan(hyperedge_counts.begin(), hyperedge_counts.end(), hyperedge_displacements.begin(), 0);

    std::exclusive_scan(pin_counts.begin(), pin_counts.end(), pin_displacements.begin(), 0);

    const int num_global_hyperedges = hyperedge_displacements.back() + hyperedge_counts.back();

    const int num_global_pins = pin_displacements.back() + pin_counts.back();

    //
    // Gather hyperedge sizes.
    //
    std::vector<SInt> global_sizes(static_cast<std::size_t>(num_global_hyperedges));

    MPI_Allgatherv(
        local_sizes.data(), num_local_hyperedges, KAGEN_MPI_SINT, global_sizes.data(), hyperedge_counts.data(),
        hyperedge_displacements.data(), KAGEN_MPI_SINT, MPI_COMM_WORLD);

    //
    // Gather all pins.
    //
    std::vector<SInt> global_pins(static_cast<std::size_t>(num_global_pins));

    MPI_Allgatherv(
        local_pins.data(), num_local_pins, KAGEN_MPI_SINT, global_pins.data(), pin_counts.data(),
        pin_displacements.data(), KAGEN_MPI_SINT, MPI_COMM_WORLD);

    //
    // Sanity-check that the gathered edge sizes account for exactly
    // the number of gathered pins.
    //
    const SInt expected_global_pins = std::accumulate(global_sizes.begin(), global_sizes.end(), SInt{0});

    if (expected_global_pins != static_cast<SInt>(global_pins.size())) {
        ADD_FAILURE() << "gathered hyperedge sizes describe " << expected_global_pins << " pins, but "
                      << global_pins.size() << " pins were gathered";

        return {};
    }

    //
    // Reconstruct the canonical expanded hypergraph.
    //
    ExpandedHypergraph global_hyperedges;

    global_hyperedges.reserve(global_sizes.size());

    std::size_t pin_offset = 0;

    for (const SInt edge_size: global_sizes) {
        const std::size_t size = static_cast<std::size_t>(edge_size);

        global_hyperedges.emplace_back(global_pins.begin() + pin_offset, global_pins.begin() + pin_offset + size);

        pin_offset += size;
    }

    return global_hyperedges;
}

template <typename Factory>
Graph GenerateHypergraph(PGeneratorConfig config) {
    PEID rank;
    PEID size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Factory factory;

    config = factory.NormalizeParameters(config, rank, size, false);

    auto generator = factory.Create(config, rank, size);

    generator->Generate(GraphRepresentation::CSR);
    generator->Finalize(MPI_COMM_WORLD);

    return generator->Take();
}

} // namespace kagen::testing