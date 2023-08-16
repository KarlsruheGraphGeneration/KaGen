#ifndef CONSTANT_H
#define CONSTANT_H

#include "kagen/kagen.h"
#include "kagen/generators/edgeweights/per_edge_weight_generator.h"
#include "kagen/generators/edgeweights/edge_weight_generator.h"

namespace kagen {
    class ConstantEdgeWeightGenerator : public PerEdgeWeightGenerator<ConstantEdgeWeightGenerator> {
    public:
        SSInt GenerateEdgeWeight(SInt u, std::tuple<HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat> cv);
        SSInt GenerateEdgeWeight3D(SInt u, std::tuple<HPFloat, HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat, HPFloat> cv);
    };

    class ConstantEdgeWeightGeneratorFactory : public EdgeWeightGeneratorFactory {
    public:
        std::unique_ptr<EdgeWeightGenerator> Create() const final;
    };
}

#endif