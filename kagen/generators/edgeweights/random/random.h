#ifndef RANDOM_H
#define RANDOM_H


#include "kagen/kagen.h"
#include "kagen/generators/edgeweights/per_edge_weight_generator.h"
#include "kagen/generators/edgeweights/edge_weight_generator.h"

namespace kagen {
    class RandomEdgeWeightGenerator : public PerEdgeWeightGenerator<RandomEdgeWeightGenerator> {
    public:
        SSInt GenerateEdgeWeight(SInt u, std::tuple<HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat> cv);
        SSInt GenerateEdgeWeight3D(SInt u, std::tuple<HPFloat, HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat, HPFloat> cv);
    };

    class RandomEdgeWeightGeneratorFactory : public EdgeWeightGeneratorFactory {
    public:
        std::unique_ptr<EdgeWeightGenerator> Create() const final;
    };
}

#endif // RANDOM_H