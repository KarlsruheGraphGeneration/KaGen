#pragma once

#include "kagen/comm/comm.h"
#include "kagen/comm/comm_types.h"
#include "kagen/context.h"

#include <gtest/gtest.h>

#include <cmath>
#include <numeric>
#include <utility>

namespace kagen::testing {
inline void GatherCoordinates2D(const Graph& local_graph, Graph& global_graph, Comm& comm) {
    PEID size = comm.Size();

    // Count the amount of coordinates per PE
    std::vector<HPFloat> x_values;
    std::vector<HPFloat> y_values;
    for (const auto& coord: local_graph.coordinates.first) {
        auto [x, y] = coord;
        x_values.push_back(x);
        y_values.push_back(y);
    }

    // Exchange amount of coordinates between PEs
    const int        local_size = x_values.size();
    std::vector<int> sizes(size);
    comm.Allgather(&local_size, 1, CommDatatype::INT, sizes.data(), 1, CommDatatype::INT);

    // Gather amount of vertices
    const int total_size = std::accumulate(sizes.begin(), sizes.end(), 0);

    // Exchange coordinates
    std::vector<int> displacements(size);
    std::exclusive_scan(sizes.begin(), sizes.end(), displacements.begin(), 0);

    global_graph.coordinates.first.resize(total_size);

    std::vector<HPFloat> recv_x_values(total_size);
    std::vector<HPFloat> recv_y_values(total_size);

    comm.Allgatherv(
        x_values.data(), local_size, CommDatatype::LONG_DOUBLE, recv_x_values.data(), sizes.data(),
        displacements.data(), CommDatatype::LONG_DOUBLE);
    comm.Allgatherv(
        y_values.data(), local_size, CommDatatype::LONG_DOUBLE, recv_y_values.data(), sizes.data(),
        displacements.data(), CommDatatype::LONG_DOUBLE);

    for (int i = 0; i < total_size; i++) {
        global_graph.coordinates.first[i] = {recv_x_values[i], recv_y_values[i]};
    }
}

inline void GatherCoordinates3D(const Graph& local_graph, Graph& global_graph, Comm& comm) {
    PEID size = comm.Size();

    // Count the amount of coordinates per PE
    std::vector<HPFloat> x_values;
    std::vector<HPFloat> y_values;
    std::vector<HPFloat> z_values;
    for (const auto& coord: local_graph.coordinates.second) {
        auto [x, y, z] = coord;
        x_values.push_back(x);
        y_values.push_back(y);
        z_values.push_back(z);
    }

    // Exchange amount of coordinates between PEs
    const int        local_size = x_values.size();
    std::vector<int> sizes(size);
    comm.Allgather(&local_size, 1, CommDatatype::INT, sizes.data(), 1, CommDatatype::INT);

    // Gather amount of vertices
    const int total_size = std::accumulate(sizes.begin(), sizes.end(), 0);

    // Exchange coordinates
    std::vector<int> displacements(size);
    std::exclusive_scan(sizes.begin(), sizes.end(), displacements.begin(), 0);

    global_graph.coordinates.second.resize(total_size);

    std::vector<HPFloat> recv_x_values(total_size);
    std::vector<HPFloat> recv_y_values(total_size);
    std::vector<HPFloat> recv_z_values(total_size);

    comm.Allgatherv(
        x_values.data(), local_size, CommDatatype::LONG_DOUBLE, recv_x_values.data(), sizes.data(),
        displacements.data(), CommDatatype::LONG_DOUBLE);
    comm.Allgatherv(
        y_values.data(), local_size, CommDatatype::LONG_DOUBLE, recv_y_values.data(), sizes.data(),
        displacements.data(), CommDatatype::LONG_DOUBLE);
    comm.Allgatherv(
        z_values.data(), local_size, CommDatatype::LONG_DOUBLE, recv_z_values.data(), sizes.data(),
        displacements.data(), CommDatatype::LONG_DOUBLE);

    for (int i = 0; i < total_size; i++) {
        global_graph.coordinates.second[i] = {recv_x_values[i], recv_y_values[i], recv_z_values[i]};
    }
}

inline void GatherWeights(const Graph& local_graph, Graph& global_graph, Comm& comm) {
    PEID size = comm.Size();

    int has_vertex_weights = !local_graph.vertex_weights.empty();
    comm.Allreduce(COMM_IN_PLACE, &has_vertex_weights, 1, CommDatatype::INT, CommOp::MAX);
    if (has_vertex_weights) {
        const int        num_local_vertices = local_graph.vertex_range.second - local_graph.vertex_range.first;
        std::vector<int> recvcounts(size);
        comm.Allgather(&num_local_vertices, 1, CommDatatype::INT, recvcounts.data(), 1, CommDatatype::INT);
        std::vector<int> displs(size);
        std::exclusive_scan(recvcounts.begin(), recvcounts.end(), displs.begin(), 0);

        global_graph.vertex_weights.resize(recvcounts.back() + displs.back());
        comm.Allgatherv(
            local_graph.vertex_weights.data(), num_local_vertices, CommDatatype::LONG_LONG,
            global_graph.vertex_weights.data(), recvcounts.data(), displs.data(), CommDatatype::LONG_LONG);
    }

    int has_edge_weights = !local_graph.edge_weights.empty();
    comm.Allreduce(COMM_IN_PLACE, &has_edge_weights, 1, CommDatatype::INT, CommOp::MAX);
    if (has_edge_weights) {
        const int        num_local_edges = std::max<int>(local_graph.edges.size(), local_graph.adjncy.size());
        std::vector<int> recvcounts(size);
        comm.Allgather(&num_local_edges, 1, CommDatatype::INT, recvcounts.data(), 1, CommDatatype::INT);
        std::vector<int> displs(size);
        std::exclusive_scan(recvcounts.begin(), recvcounts.end(), displs.begin(), 0);

        global_graph.edge_weights.resize(recvcounts.back() + displs.back());
        comm.Allgatherv(
            local_graph.edge_weights.data(), num_local_edges, CommDatatype::LONG_LONG,
            global_graph.edge_weights.data(), recvcounts.data(), displs.data(), CommDatatype::LONG_LONG);
    }
}

inline Graph GatherEdgeLists(const Graph& local_graph, Comm& comm) {
    PEID size = comm.Size();

    const int        num_local_edges = local_graph.edges.size() * 2;
    std::vector<int> recvcounts(size);
    comm.Allgather(&num_local_edges, 1, CommDatatype::INT, recvcounts.data(), 1, CommDatatype::INT);

    std::vector<int> displs(size);
    std::exclusive_scan(recvcounts.begin(), recvcounts.end(), displs.begin(), 0u);
    const SInt num_global_edges = displs.back() + recvcounts.back();

    Graph global_graph;
    global_graph.representation = local_graph.representation;
    global_graph.edges.resize(num_global_edges / 2);
    comm.Allgatherv(
        local_graph.edges.data(), num_local_edges, CommDatatype::UNSIGNED_LONG_LONG, global_graph.edges.data(),
        recvcounts.data(), displs.data(), CommDatatype::UNSIGNED_LONG_LONG);

    GatherWeights(local_graph, global_graph, comm);
    GatherCoordinates2D(local_graph, global_graph, comm);
    GatherCoordinates3D(local_graph, global_graph, comm);

    return global_graph;
}

inline Graph GatherCSR(const Graph& local_graph, Comm& comm) {
    EXPECT_EQ(local_graph.xadj.size() - 1, local_graph.vertex_range.second - local_graph.vertex_range.first);
    EXPECT_EQ(local_graph.adjncy.size(), local_graph.xadj.back());

    PEID size = comm.Size();

    const int        num_local_vertices = local_graph.xadj.size() - 1;
    std::vector<int> degree_recvcounts(size);
    comm.Allgather(&num_local_vertices, 1, CommDatatype::INT, degree_recvcounts.data(), 1, CommDatatype::INT);
    std::vector<int> degree_displs(size);
    std::exclusive_scan(degree_recvcounts.begin(), degree_recvcounts.end(), degree_displs.begin(), 0);

    const int        num_local_edges = local_graph.adjncy.size();
    std::vector<int> edges_recvcounts(size);
    comm.Allgather(&num_local_edges, 1, CommDatatype::INT, edges_recvcounts.data(), 1, CommDatatype::INT);
    std::vector<int> edges_displs(size);
    std::exclusive_scan(edges_recvcounts.begin(), edges_recvcounts.end(), edges_displs.begin(), 0);

    std::vector<SInt> local_degrees(num_local_vertices);
    for (SInt u = 0; u < num_local_vertices; ++u) {
        local_degrees[u] = local_graph.xadj[u + 1] - local_graph.xadj[u];
    }

    Graph global_graph;
    global_graph.representation = local_graph.representation;
    global_graph.xadj.resize(degree_recvcounts.back() + degree_displs.back() + 1);
    global_graph.adjncy.resize(edges_recvcounts.back() + edges_displs.back());

    comm.Allgatherv(
        local_degrees.data(), num_local_vertices, CommDatatype::UNSIGNED_LONG_LONG, global_graph.xadj.data(),
        degree_recvcounts.data(), degree_displs.data(), CommDatatype::UNSIGNED_LONG_LONG);
    comm.Allgatherv(
        local_graph.adjncy.data(), num_local_edges, CommDatatype::UNSIGNED_LONG_LONG, global_graph.adjncy.data(),
        edges_recvcounts.data(), edges_displs.data(), CommDatatype::UNSIGNED_LONG_LONG);

    std::exclusive_scan(global_graph.xadj.begin(), global_graph.xadj.end(), global_graph.xadj.begin(), 0);

    EXPECT_EQ(global_graph.xadj.back(), global_graph.adjncy.size());

    GatherWeights(local_graph, global_graph, comm);
    GatherCoordinates2D(local_graph, global_graph, comm);
    GatherCoordinates3D(local_graph, global_graph, comm);

    return global_graph;
}

inline Graph GatherGraph(const Graph& local_graph, Comm& comm) {
    switch (local_graph.representation) {
        case GraphRepresentation::CSR:
            return GatherCSR(local_graph, comm);

        case GraphRepresentation::EDGE_LIST:
            return GatherEdgeLists(local_graph, comm);
    }

    __builtin_unreachable();
}
} // namespace kagen::testing
