#include <random>
#include "random.h"
#include "kagen/kagen.h"

namespace kagen {
    SSInt RandomEdgeWeightGenerator::GenerateEdgeWeight(SInt u, std::tuple<HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat> cv) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dist(0, 100); // ToDo: What range

        return dist(gen);
    }

    SSInt RandomEdgeWeightGenerator::GenerateEdgeWeight3D(SInt u, std::tuple<HPFloat, HPFloat, HPFloat> cu, SInt v, std::tuple<HPFloat, HPFloat, HPFloat> cv) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dist(0, 100); // ToDo: What range

        return dist(gen);
    }

    std::unique_ptr<EdgeWeightGenerator> RandomEdgeWeightGeneratorFactory::Create() const {
        return std::make_unique<RandomEdgeWeightGenerator>();
    }
}
