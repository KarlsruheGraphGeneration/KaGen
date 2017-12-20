/******************************************************************************
 * generator_io.h
 *
 * Source of the graph generator
 ******************************************************************************
 * Copyright (C) 2016 Sebastian Lamm <lamm@ira.uka.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef _GENERATOR_IO_H_
#define _GENERATOR_IO_H_

#include <mpi.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <vector>

#include "generator_config.h"

template <typename T>
struct identity {
  typedef T type;
};

template <typename Edge = std::tuple<SInt, SInt>>
class GeneratorIO {
 public:
  GeneratorIO(const PGeneratorConfig& config) : config_(config), local_num_edges_(0) {
    dist_.resize(config_.dist_size);
  }

  inline void UpdateDist(SInt node_id) {
    // if ((CRCHash::hash(node_id) % config_.n) < dist_.size()) dist_[node_id]++;
    if (node_id < dist_.size()) dist_[node_id]++;
    local_num_edges_++;
  }

  void OutputDist() const {
    // Exchange local dist
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<SInt> global_dist(dist_.size(), 0);
    MPI_Reduce(&dist_[0], &global_dist[0], dist_.size(), MPI_LONG, MPI_SUM,
               ROOT, MPI_COMM_WORLD);
    if (rank == ROOT) {
      FILE* fout = fopen(config_.output_file.c_str(), "w+");
      for (SInt i = 0; i < global_dist.size(); ++i) {
        fprintf(fout, "%llu\n", global_dist[i]);
      }
      fclose(fout);
    }
  }

  void ReserveEdges(SInt num_edges) { edges_.reserve(num_edges); }

  template <typename... Args>
  inline void PushEdge(Args... args) {
    edges_.emplace_back(std::make_tuple(args...));
    local_num_edges_++;
  }

  void OutputEdges() const { 
#ifdef SINGLE_LIST
    GatherPrint(identity<Edge>());
#else
    Print(identity<Edge>()); 
#endif
  }

  SInt NumEdges() const { 
    return edges_.size() > 0 ? edges_.size() : local_num_edges_/2; 
  }

 private:
  PGeneratorConfig config_;

  std::vector<SInt> dist_;
  std::vector<Edge> edges_;

  SInt local_num_edges_;

  void GatherPrint(identity<std::tuple<SInt, SInt>>) const {
    // Exchange local dist
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Gather number of edges for each PE
    std::vector<int> displ(size);
    std::vector<int> num_edges(size);
    MPI_Gather(&local_num_edges_, 1, MPI_INT, &num_edges[0], 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    int current_displ = 0;
    int total_num_edges = 0;
    if (rank == ROOT) {
      for (SInt i = 0; i < num_edges.size(); ++i) {
        displ[i] = current_displ;
        total_num_edges += num_edges[i];
        current_displ = total_num_edges;
      }
    }

    // Gather actual edges
    MPI_Datatype MPI_EDGE;
    MPI_Type_vector(1, 2, 0, MPI_LONG, &MPI_EDGE);
    MPI_Type_commit(&MPI_EDGE);
    std::vector<Edge> edges(total_num_edges);
    MPI_Gatherv(&edges_[0], local_num_edges_, MPI_EDGE, &edges[0], &num_edges[0], &displ[0], MPI_EDGE, ROOT, MPI_COMM_WORLD);


    if (rank == ROOT) {
      // Sort edges and remove duplicates
      std::sort(std::begin(edges), std::end(edges));
      edges.erase(unique(edges.begin(), edges.end()), edges.end());
      
      // Output edges
      FILE* fout = fopen(config_.output_file.c_str(), "w+");
      fprintf(fout, "%llu %llu\n", config_.n, edges.size());
      for (auto edge : edges) fprintf(fout, "%llu %llu\n", std::get<0>(edge), std::get<1>(edge));
      fclose(fout);
    }
  }

  // node id output
  void Print(identity<std::tuple<SInt, SInt>>) const {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    FILE* fout =
        fopen((config_.output_file + std::to_string(rank)).c_str(), "w+");
    for (auto edge : edges_) {
      fprintf(fout, "%llu %llu\n", std::get<0>(edge), std::get<1>(edge));
    }
    fclose(fout);
  };

  // ABUSE: adjacency list output
  void Print(identity<std::tuple<SInt, std::vector<SInt>>>) const {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // fugly but saves memory as this is a const method
    auto& nodes = const_cast<std::vector<Edge>&>(edges_);

    // sort edges by node id
    std::sort(std::begin(nodes), std::end(nodes),
              [](const Edge& t1, const Edge& t2) {
                return std::get<0>(t1) <
                       std::get<0>(t2);  // or use a custom compare function
              });

    // compute edge count
    SInt edgeCount = std::accumulate(
        std::begin(nodes), std::end(nodes), SInt(0),
        [](SInt a, const Edge& b) { return a + std::get<1>(b).size(); });

    FILE* fout =
        fopen((config_.output_file + std::to_string(rank)).c_str(), "w+");
    fprintf(fout, "%lu %llu\n", edges_.size(), edgeCount);

    for (auto& node : nodes) {
      auto& edges = std::get<1>(node);
      edgeCount += edges.size();
      std::sort(std::begin(edges), std::end(edges));

      std::string sep = "";
      for (const auto& edge : edges) {
        fprintf(fout, "%s%llu", sep.c_str(), edge);
        sep = " ";
      }
      fprintf(fout, "\n");
    }

    fclose(fout);
  };
};

#endif
