#pragma once

#include "kagen/generators/generator.h"

#include <mpi.h>

namespace kagen {
class Graph500Generator : public virtual Generator, private EdgeListOnlyGenerator {
public:
    Graph500Generator(const PGeneratorConfig& config);

protected:
    void FinalizeEdgeList(MPI_Comm comm) final;

    inline void PushLocalEdge(const SInt from, const SInt to) {
        if (config_.self_loops || from != to) {
            local_edges_.emplace_back(from, to);
        }
        if (!config_.directed && from != to) {
            local_edges_.emplace_back(to, from);
        }
    }

private:
    const PGeneratorConfig& config_;
    Edgelist32              local_edges_;
};
} // namespace kagen
