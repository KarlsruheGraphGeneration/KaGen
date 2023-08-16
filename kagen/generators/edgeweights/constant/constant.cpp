#include "constant.h"
#include "kagen/kagen.h"

namespace kagen {
    SSInt ConstantEdgeWeightGenerator::GenerateEdgeWeight(SInt u, std::tuple<HPFloat, HPFloat> cu, SInt v,
                                                          std::tuple<HPFloat, HPFloat> cv) {
        return 1;
    }

    SSInt ConstantEdgeWeightGenerator::GenerateEdgeWeight3D(SInt u, std::tuple<HPFloat, HPFloat, HPFloat> cu, SInt v,
                                                            std::tuple<HPFloat, HPFloat, HPFloat> cv) {
        return 1;
    }


    std::unique_ptr<EdgeWeightGenerator> ConstantEdgeWeightGeneratorFactory::Create() const {
        return std::make_unique<ConstantEdgeWeightGenerator>();
    }
}