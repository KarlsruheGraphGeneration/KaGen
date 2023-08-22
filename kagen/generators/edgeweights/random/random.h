#ifndef RANDOM_H
#define RANDOM_H


#include "kagen/kagen.h"
#include "kagen/generators/edgeweights/per_edge_weight_generator.h"
#include "kagen/generators/edgeweights/edge_weight_generator.h"
#include "kagen/context.h"

namespace kagen {
    class RandomEdgeWeightGenerator : public PerEdgeWeightGenerator<RandomEdgeWeightGenerator> {
    public:
        RandomEdgeWeightGenerator(EdgeWeightConfig config);

        SSInt GenerateEdgeWeight(SInt u, SInt v);

        SSInt GenerateEdgeWeight2D(SInt u, std::tuple<HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat> cv);

        SSInt GenerateEdgeWeight3D(SInt u, std::tuple<HPFloat, HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat, HPFloat> cv);

    private:
        const EdgeWeightConfig config_;
    };

    class RandomEdgeWeightGeneratorFactory : public EdgeWeightGeneratorFactory {
    public:
        std::unique_ptr<EdgeWeightGenerator> Create(EdgeWeightConfig config) const final;
    };
}

#endif // RANDOM_H