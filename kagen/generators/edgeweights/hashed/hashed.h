#ifndef HASHED_H
#define HASHED_H

#include "kagen/kagen.h"
#include "kagen/generators/edgeweights/per_edge_weight_generator.h"
#include "kagen/generators/edgeweights/edge_weight_generator.h"

namespace kagen {
    class HashedEdgeWeightGenerator : public PerEdgeWeightGenerator<HashedEdgeWeightGenerator> {
    public:
        SSInt GenerateEdgeWeight(SInt u, std::tuple<HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat> cv);
        SSInt GenerateEdgeWeight3D(SInt u, std::tuple<HPFloat, HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat, HPFloat> cv);
    };

    class HashedEdgeWeightGeneratorFactory : public EdgeWeightGeneratorFactory {
    public:
        std::unique_ptr<EdgeWeightGenerator> Create() const final;
    };
}

#endif


//#pragma once
//
//#include "kagen/kagen.h"
//#include "kagen/generators/edgeweights/edgeweightgenerator.h"
//
//namespace kagen {
//    class HashedEdgeWeightGenerator : public EdgeWeightGenerator {
//    public:
//        EdgeWeights GenerateEdgeWeights(Edgelist edges, Coordinates coor);
//    };
//
//    class HashedEdgeWeightGeneratorFactory : public EdgeWeightGeneratorFactory {
//    public:
//        std::unique_ptr<EdgeWeightGenerator> Create(const PGeneratorConfig &config) const;
//    };
//
//} // namespace kagen
