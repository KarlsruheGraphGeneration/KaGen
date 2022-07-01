#pragma once

#include "kagen/generators/generator.h"

#include <mpi.h>

namespace kagen {
class Graph500Generator : public Generator {
public:
    Graph500Generator(const PGeneratorConfig& config_, MPI_Comm comm);

protected:
    inline void PushLocalEdge(const int from, const int to) {
        if (config_.self_loops || from != to) {
            local_edges_.emplace_back(from, to);
        }
        if (!config_.directed && from != to) {
            local_edges_.emplace_back(to, from);
        }
    }

    void DistributeRoundRobin(SInt n);

private:
    const PGeneratorConfig&           config_;
    MPI_Comm                          comm_;
    std::vector<std::tuple<int, int>> local_edges_;
};
} // namespace kagen
