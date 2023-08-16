#ifndef PER_EDGE_WEIGHT_GENERATOR_H
#define PER_EDGE_WEIGHT_GENERATOR_H

#include "kagen/kagen.h"
#include "kagen/generators/edgeweights/edge_weight_generator.h"

namespace kagen {
    template<typename Derived>
    class PerEdgeWeightGenerator : public EdgeWeightGenerator {
    public:
        EdgeWeights GenerateEdgeWeights(const Edgelist &edgelist, const Coordinates &coordinates) final {
            EdgeWeights weights;
            weights.reserve(edgelist.size());

            for (const auto &[u, v]: edgelist) {
                SSInt weight;
                if(!coordinates.first.empty()) {
                    weight = _GenerateEdgeWeight(u, coordinates.first[u], v, coordinates.first[v]);
                } else {
                    weight = _GenerateEdgeWeight3D(u, coordinates.second[u], v, coordinates.second[v]);
                }
                weights.push_back(weight);
            }

            return weights;
        }

        EdgeWeights GenerateEdgeWeights(const XadjArray &xadj, const AdjncyArray &adjncy,
                                        const Coordinates &coordinates) final {
            EdgeWeights weights;
            weights.reserve(adjncy.size());

            for (int u = 0; u + 1 < xadj.size(); ++u) {
                for (int e = xadj[u]; e < xadj[u + 1]; ++e) {
                    int v = adjncy[e];
                    SSInt weight;
                    if(!coordinates.first.empty()) {
                        weight = _GenerateEdgeWeight(u, coordinates.first[u], v, coordinates.first[v]);
                    } else {
                        weight = _GenerateEdgeWeight3D(u, coordinates.second[u], v, coordinates.second[v]);
                    }
                    weights.push_back(weight);
                }
            }

            return weights;
        }

    private:
        SSInt _GenerateEdgeWeight(SInt u, std::tuple<HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat> cv) {
            return static_cast<Derived *>(this)->GenerateEdgeWeight(u, cu, v, cv);
        }

        SSInt _GenerateEdgeWeight3D(SInt u, std::tuple<HPFloat, HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat, HPFloat> cv) {
            return static_cast<Derived *>(this)->GenerateEdgeWeight3D(u, cu, v, cv);
        }
    };

}

#endif  // PER_EDGE_WEIGHT_GENERATOR_H