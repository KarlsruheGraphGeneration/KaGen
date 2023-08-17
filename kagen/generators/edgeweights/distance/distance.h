#ifndef DISTANCE_H
#define DISTANCE_H

#include "kagen/kagen.h"
#include "kagen/generators/edgeweights/per_edge_weight_generator.h"
#include "kagen/generators/edgeweights/edge_weight_generator.h"

namespace kagen {
    class DistanceEdgeWeightGenerator : public PerEdgeWeightGenerator<DistanceEdgeWeightGenerator> {
    public:
        SSInt GenerateEdgeWeight(SInt u, SInt v);

        SSInt GenerateEdgeWeight2D(SInt u, std::tuple<HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat> cv);

        SSInt GenerateEdgeWeight3D(SInt u, std::tuple<HPFloat, HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat, HPFloat> cv);
    };

    class DistanceEdgeWeightGeneratorFactory : public EdgeWeightGeneratorFactory {
    public:
        std::unique_ptr<EdgeWeightGenerator> Create() const final;
    };
}

#endif