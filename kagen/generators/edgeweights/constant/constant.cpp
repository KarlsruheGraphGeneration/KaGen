#include "constant.h"
#include "kagen/kagen.h"
#include "kagen/context.h"

namespace kagen {
    ConstantEdgeWeightGenerator::ConstantEdgeWeightGenerator(EdgeWeightConfig config) : config_(config) {
    }

    SSInt ConstantEdgeWeightGenerator::GenerateEdgeWeight(SInt u, SInt v) {
        return 1;
    }

    SSInt ConstantEdgeWeightGenerator::GenerateEdgeWeight2D(SInt u, std::tuple<HPFloat, HPFloat> cu, SInt v,
                                                          std::tuple<HPFloat, HPFloat> cv) {
        return 1;
    }

    SSInt ConstantEdgeWeightGenerator::GenerateEdgeWeight3D(SInt u, std::tuple<HPFloat, HPFloat, HPFloat> cu, SInt v,
                                                            std::tuple<HPFloat, HPFloat, HPFloat> cv) {
        return 1;
    }


    std::unique_ptr<EdgeWeightGenerator>
    ConstantEdgeWeightGeneratorFactory::Create(EdgeWeightConfig config) const {
        return std::make_unique<ConstantEdgeWeightGenerator>(config);
    }
}