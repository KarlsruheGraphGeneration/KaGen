#pragma once

#include "kagen/kagen.h"

#include <mpi.h>

namespace kagen {
/**
 * @brief Adds any missing reverse edges that goes between different PEs.
 * More precisely, for any edge `(u, v)` where `u` is assigned to another PE than `v`, this function adds the reverse
 * edge `(v, u)` if it is not already present.
 * **Note:** in order to do so, the edge list is sorted and duplicate edges are removed.
 *
 * @param edge_list The edge list to add reverse edges to.
 * @param vertex_range The vertex range assigned to this PE.
 * @param comm The MPI communicator.
 */
void AddNonlocalReverseEdges(Edgelist& edge_list, VertexRange vertex_range, MPI_Comm comm);

/**
 * @brief Redistributes the edges such that each PE owns the edges that it should owned due to the given vertex range.
 * @param edge_list The edge list to redistribute.
 * @param vertex_range The vertex range assigned to this PE.
 * @param comm The MPI communicator.
 */
void RedistributeEdgesByVertexRange(Edgelist& edge_list, VertexRange vertex_range, MPI_Comm comm);

/**
 * @brief Assigns vertices round-robin to PEs and redistributes the edges accordingly.
 * More precisely, this function assigns vertex `v` to PE `v mod <nproc>`, relabels the vertices, and
 * migrates the edges accordingly.
 * **Note:** In order to do so, both edge lists are sorted and duplicate edges are removed.
 *
 * @param source The edge list to redistribute.
 * @param destination The edge list to store the redistributed edges in.
 * @param n The number of vertices in the graph.
 * @param comm The MPI communicator.
 * @return The vertex range assigned to this PE.
 */
VertexRange RedistributeEdgesRoundRobin(Edgelist32& source, Edgelist& destination, SInt n, MPI_Comm comm);
} // namespace kagen
