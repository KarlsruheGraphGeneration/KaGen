#include <random>
#include "random.h"
#include "kagen/kagen.h"

namespace kagen {
    RandomEdgeWeightGenerator::RandomEdgeWeightGenerator(EdgeWeightConfig config) : config_(config) {
    }

    SSInt RandomEdgeWeightGenerator::GenerateEdgeWeight(SInt u, SInt v) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dist(config_.weight_minimum, config_.weight_maximum);

        return dist(gen);
    }

    SSInt RandomEdgeWeightGenerator::GenerateEdgeWeight2D(SInt u, std::tuple<HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat> cv) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dist(config_.weight_minimum, config_.weight_maximum);

        return dist(gen);
    }

    SSInt RandomEdgeWeightGenerator::GenerateEdgeWeight3D(SInt u, std::tuple<HPFloat, HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat, HPFloat> cv) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dist(config_.weight_minimum, config_.weight_maximum);

        return dist(gen);
    }

    std::unique_ptr<EdgeWeightGenerator> RandomEdgeWeightGeneratorFactory::Create(EdgeWeightConfig config) const {
        return std::make_unique<RandomEdgeWeightGenerator>(config);
    }
}
