#pragma once

#include "kagen/context.h"
#include "kagen/edge_range.h"
#include "kagen/edgeweight_generators/edge_weight_generator.h"
#include "kagen/kagen.h"

namespace kagen {
template <typename Derived>
class PerEdgeWeightGenerator : public EdgeWeightGenerator {
public:
    PerEdgeWeightGenerator(VertexRange vertex_range) :
          vertex_range_{vertex_range} {}

    void GenerateEdgeWeights(const Edgelist& edgelist, EdgeWeights& weights) final {
        EdgeRange edge_range(edgelist);
        GenerateEdgeWeightsImpl(edge_range, weights);
    }

    void GenerateEdgeWeights(const XadjArray& xadj, const AdjncyArray& adjncy, EdgeWeights& weights) final {
        EdgeRange edge_range(xadj, adjncy, vertex_range_);
        GenerateEdgeWeightsImpl(edge_range, weights);
    }

private:
    VertexRange vertex_range_;

    void GenerateEdgeWeightsImpl(EdgeRange range, EdgeWeights& weights) {
        weights.reserve(range.size());
        for (auto it = range.begin(); it != range.end(); ++it) {
            auto [u, v]  = *it;
            SSInt weight = static_cast<Derived*>(this)->GenerateEdgeWeight(u, v);
            weights.push_back(weight);
        }
    }
};
} // namespace kagen
