#pragma once

#include "kagen/context.h"
#include "kagen/kagen.h"
#include "kagen/vertexweight_generators/vertex_weight_generator.h"

namespace kagen {
/*!
 * Locally draws a weight for each vertex from a pseudorandom number generator.
 */
class UniformRandomVertexWeightGenerator : public VertexWeightGenerator {
public:
    UniformRandomVertexWeightGenerator(VertexWeightConfig config, MPI_Comm comm);
    void GenerateVertexWeights(const VertexRange& vertex_range, const Edgelist& edgelist, VertexWeights& weights) final;
    void GenerateVertexWeights(
        const VertexRange& vertex_range, const XadjArray& xadj, const AdjncyArray& adjncy,
        VertexWeights& weights) final;

private:
    const VertexWeightConfig config_;
    MPI_Comm                 comm_;
};
} // namespace kagen
