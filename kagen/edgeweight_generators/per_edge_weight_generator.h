#pragma once

#include "kagen/edgeweight_generators/edge_weight_generator.h"
#include "kagen/kagen.h"

namespace kagen {
template <typename Derived>
class PerEdgeWeightGenerator : public EdgeWeightGenerator {
public:
    EdgeWeights GenerateEdgeWeights(const Edgelist& edgelist) final {
        EdgeWeights weights;
        weights.reserve(edgelist.size());

        for (const auto& [u, v]: edgelist) {
            SSInt weight = static_cast<Derived*>(this)->GenerateEdgeWeight(u, v);
            weights.push_back(weight);
        }
        return weights;
    }

    EdgeWeights
    GenerateEdgeWeights(const XadjArray& xadj, const AdjncyArray& adjncy) final {
        EdgeWeights weights;
        weights.reserve(adjncy.size());

        for (SInt u = 0; u + 1 < xadj.size(); ++u) {
            for (SInt e = xadj[u]; e < xadj[u + 1]; ++e) {
                SInt   v      = adjncy[e];
                SSInt weight = static_cast<Derived*>(this)->GenerateEdgeWeight(u, v);
                weights.push_back(weight);
                weights.push_back(weight);
            }
        }

        return weights;
    }
};
} // namespace kagen
