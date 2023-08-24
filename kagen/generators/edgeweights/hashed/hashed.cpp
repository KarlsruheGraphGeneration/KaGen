#include "hashed.h"
#include "kagen/kagen.h"
#include "xxhash.h"
#include "kagen/context.h"

namespace kagen {
    HashedEdgeWeightGenerator::HashedEdgeWeightGenerator(EdgeWeightConfig config) : config_(config) {
    }

    SSInt HashedEdgeWeightGenerator::GenerateEdgeWeight(SInt u, SInt v) {
        SSInt hash1, hash2;

        hash1 = XXH64(&u, 1, 0);
        hash2 = XXH64(&v, 1, 0);
        return (hash1 ^ hash2) % config_.weight_maximum + config_.weight_minimum; // XOR combination
    }

    SSInt HashedEdgeWeightGenerator::GenerateEdgeWeight2D(SInt u, std::tuple<HPFloat, HPFloat> cu, SInt v,
                                                          std::tuple<HPFloat, HPFloat> cv) {
        SSInt hash1, hash2;

        hash1 = XXH64(&u, 1, 0);
        hash2 = XXH64(&v, 1, 0);
        return (hash1 ^ hash2) % config_.weight_maximum + config_.weight_minimum; // XOR combination
    }

    SSInt HashedEdgeWeightGenerator::GenerateEdgeWeight3D(SInt u, std::tuple<HPFloat, HPFloat, HPFloat> cu, SInt v,
                                                          std::tuple<HPFloat, HPFloat, HPFloat> cv) {
        SSInt hash1, hash2;

        hash1 = XXH64(&u, 1, 0);
        hash2 = XXH64(&v, 1, 0);
        return (hash1 ^ hash2) % config_.weight_maximum + config_.weight_minimum; // XOR combination
    }

    std::unique_ptr<EdgeWeightGenerator>
    HashedEdgeWeightGeneratorFactory::Create(kagen::EdgeWeightConfig config) const {
        return std::make_unique<HashedEdgeWeightGenerator>(config);
    }
}


//#include "hashed.h"
//#include "kagen/kagen.h"
//
//namespace kagen {
//    EdgeWeights HashedEdgeWeightGenerator::GenerateEdgeWeights(Edgelist edges, Coordinates coor) {
//        SSInt combinedHash, hash1, hash2;
//        EdgeWeights edgeWeights;
//        for (const auto &pair: edges) {
//            hash1 = XXH64(&pair.first, sizeof(pair.first), 0);
//            hash2 = XXH64(&pair.second, sizeof(pair.second), 0);
//            std::cout << "Edge (" << &pair.first << ", " << &pair.second << ") " << hash1
//                      << "Hash2" << hash2 << std::endl;
//            combinedHash = hash1 ^ hash2; // XOR combination
//            edgeWeights.emplace_back(combinedHash % edgeWeights.size());
//        }
//        return edgeWeights;
//    }
//
//    std::unique_ptr<EdgeWeightGenerator> HashedEdgeWeightGeneratorFactory::Create(const PGeneratorConfig &config) const {
//        return std::make_unique<HashedEdgeWeightGenerator>();
//    }
//} // namespace kagen
