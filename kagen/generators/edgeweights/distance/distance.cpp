#include "distance.h"
#include "kagen/kagen.h"
#include <exception>
#include <cmath>

namespace kagen {
    SSInt DistanceEdgeWeightGenerator::GenerateEdgeWeight(SInt u, std::tuple<HPFloat, HPFloat> cu, SInt v,
                                                          std::tuple<HPFloat, HPFloat> cv) {
        auto [x1, y1] = cu;
        auto [x2, y2] = cv;
        return std::hypot(x1 - x2, y1 - y2) * 100;
    }

    SSInt DistanceEdgeWeightGenerator::GenerateEdgeWeight3D(SInt u, std::tuple<HPFloat, HPFloat, HPFloat> cu, SInt v,
                                                            std::tuple<HPFloat, HPFloat, HPFloat> cv) {
        auto [x1, y1, z1] = cu;
        auto [x2, y2, z2] = cv;
        return std::hypot(x1 - x2, y1 - y2, z1 - z2) * 100;
    }


    std::unique_ptr<EdgeWeightGenerator> DistanceEdgeWeightGeneratorFactory::Create() const {
        return std::make_unique<DistanceEdgeWeightGenerator>();
    }
}
