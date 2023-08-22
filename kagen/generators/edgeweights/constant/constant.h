#ifndef CONSTANT_H
#define CONSTANT_H

#include "kagen/kagen.h"
#include "kagen/generators/edgeweights/per_edge_weight_generator.h"
#include "kagen/generators/edgeweights/edge_weight_generator.h"
#include "kagen/context.h"

namespace kagen {
    class ConstantEdgeWeightGenerator : public PerEdgeWeightGenerator<ConstantEdgeWeightGenerator> {
    public:
        ConstantEdgeWeightGenerator(EdgeWeightConfig config);

        SSInt GenerateEdgeWeight(SInt u, SInt v);

        SSInt GenerateEdgeWeight2D(SInt u, std::tuple<HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat> cv);

        SSInt GenerateEdgeWeight3D(SInt u, std::tuple<HPFloat, HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat, HPFloat> cv);

    private:
        const EdgeWeightConfig config_;
    };

    class ConstantEdgeWeightGeneratorFactory : public EdgeWeightGeneratorFactory {
    public:
        std::unique_ptr<EdgeWeightGenerator> Create(EdgeWeightConfig config) const final;
    };
}

#endif