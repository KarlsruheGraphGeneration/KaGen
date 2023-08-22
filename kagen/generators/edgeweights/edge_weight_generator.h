#pragma once
#ifndef EDGE_WEIGHT_GENERATOR_H
#define EDGE_WEIGHT_GENERATOR_H

#include <vector>
#include <tuple>
#include "kagen/definitions.h"
#include "kagen/context.h"

namespace kagen {
    class EdgeWeightGenerator {
    public:
        virtual ~EdgeWeightGenerator() = default;

        virtual EdgeWeights GenerateEdgeWeights(const Edgelist &edgelist, const Coordinates &coordinates) = 0;

        virtual EdgeWeights
        GenerateEdgeWeights(const XadjArray &xadj, const AdjncyArray &adjncy, const Coordinates &coordinates) = 0;
    };

    class EdgeWeightGeneratorFactory {
    public:
        virtual ~EdgeWeightGeneratorFactory() = default;

        virtual std::unique_ptr<EdgeWeightGenerator> Create(EdgeWeightConfig config) const = 0;
    };
}

#endif  // EDGE_WEIGHT_GENERATOR_H



//#pragma once
//
//#include <cstdint>
//#include <exception>
//#include <memory>
//
//#include <mpi.h>
//
//#include <cmath>
//#include "kagen/context.h"
//#include "kagen/definitions.h"
//#include "xxhash.h"
//
//namespace kagen {
//    class EdgeWeightGenerator {
//    public:
//        virtual ~EdgeWeightGenerator() {}
//
//        virtual EdgeWeights GenerateEdgeWeights(const Edgelist &edgelist, const Coordinates &coordinates) = 0;
//
//        virtual EdgeWeights GenerateEdgeWeights(const std::vector<int> &xadj, const std::vector<int> &adjncy,
//                                                     const Coordinates &coordinates) = 0;
//
//    };
//
//    class EdgeWeightGeneratorFactory {
//    public:
//        virtual ~EdgeWeightGeneratorFactory() {}
//
//        virtual std::unique_ptr<EdgeWeightGenerator> Create() const = 0;
//    };
//} // namespace kagen



//##################################################

//#pragma once
//
//#include <cstdint>
//#include <exception>
//#include <memory>
//
//#include <mpi.h>
//
//#include <cmath>
//#include "kagen/context.h"
//#include "kagen/definitions.h"
//#include "xxhash.h"
//
//namespace kagen {
//    class EdgeWeightGenerator {
//    public:
//        virtual ~EdgeWeightGenerator() {}
//
//        virtual EdgeWeights GenerateEdgeWeights(Edgelist edges, Coordinates coor);
//    };
//
//    class EdgeWeightGeneratorFactory {
//    public:
//        virtual ~EdgeWeightGeneratorFactory() {}
//
//        virtual std::unique_ptr<EdgeWeightGenerator> Create() const = 0;
//    };
//} // namespace kagen
//
//
////#pragma once
////
////#include <cstdint>
////#include <exception>
////#include <memory>
////
////#include <mpi.h>
////
////#include <cmath>
////#include "kagen/context.h"
////#include "kagen/definitions.h"
////#include "xxhash.h"
////
////namespace kagen {
////    class EdgeWeightGenerator {
////    public:
////        static EdgeWeights
////        GenerateEdgeWeights(MPI_Comm comm, const EdgeWeightType type, Edgelist edges, Coordinates coor, GraphRepresentation rep) {
////            switch (type) {
////                case EdgeWeightType::HASHED:
////                    return HashedEdgeWeights(edges);
////                case EdgeWeightType::DISTANCE:
////                    return DistanzEdgeWeights(edges, coor);
////                case EdgeWeightType::CONSTANT:
////                    return ConstantEdgeWeights(edges);
////                case EdgeWeightType::RANDOM:
////                    return RandomEdgeWeights(comm, edges);
////            }
////        }
////
////        static EdgeWeights HashedEdgeWeights(Edgelist edges) {
////            SSInt combinedHash, hash1, hash2;
////            EdgeWeights edgeWeights;
////            for (const auto &pair: edges) {
////                hash1 = XXH64(&pair.first, sizeof(pair.first), 0);
////                hash2 = XXH64(&pair.second, sizeof(pair.second), 0);
////                std::cout << "Edge (" << &pair.first << ", " << &pair.second << ") " << hash1
////                          << "Hash2" << hash2 << std::endl;
////                combinedHash = hash1 ^ hash2; // XOR combination
////                edgeWeights.emplace_back(combinedHash % edgeWeights.size());
////            }
////            return edgeWeights;
////        }
////
////        static EdgeWeights DistanzEdgeWeights(Edgelist edges, Coordinates coor) {
////            EdgeWeights edgeWeights;
////            if (!coor.first.empty()) {
////                for (SInt i = 0; i < coor.first.size(); i++) {
////                    auto [x1, y1] = coor.first[i];
////                    for (SInt j = 0; j < coor.first.size(); j++) {
////                        auto [x2, y2] = coor.first[j];
////                        // Comparing all
////                        if (i != j) {
////                            edgeWeights.emplace_back(std::hypot(x1 - x2, y1 - y2));
////                        }
////                    }
////                }
////            } else if (!coor.second.empty()) {
////                for (SInt i = 0; i < coor.first.size(); i++) {
////                    auto [x1, y1, z1] = coor.second[i];
////                    for (SInt j = 0; j < coor.first.size(); j++) {
////                        auto [x2, y2, z2] = coor.second[j];
////                        // Comparing all
////                        if (i != j) {
////                            edgeWeights.emplace_back(std::hypot(x1 - x2, y1 - y2, z1 - z2));
////                        }
////                    }
////                }
////            }
////
////            return edgeWeights;
////        }
////
////        static EdgeWeights ConstantEdgeWeights(Edgelist edges) {
////            EdgeWeights edgeweights;
////            for (const auto &pair: edges) {
////                edgeweights.emplace_back(1);
////            }
////            return edgeweights;
////        }
////
////        static EdgeWeights RandomEdgeWeights(MPI_Comm pCommunicator, Edgelist edges) {
////            return kagen::EdgeWeights();
////        }
////
////    };
////
////    class EdgeWeightGeneratorFactory {
////    public:
////        virtual ~EdgeWeightGeneratorFactory();
////
////        virtual std::unique_ptr<EdgeWeightGenerator>
////        Create(const PGeneratorConfig &config, PEID rank, PEID size) const = 0;
////    };
////} // namespace kagen
