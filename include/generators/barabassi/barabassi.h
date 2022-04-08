/*******************************************************************************
 * include/generators/barabassi/barabassi.h
 *
 * Copyright (C) 2016-2017 Christian Schulz <christian.schulz@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "definitions.h"
#include "generator_config.h"
#include "generator_io.h"
#include "hash.hpp"

namespace kagen {
class Barabassi {
public:
    Barabassi(PGeneratorConfig& config, const PEID rank)
        : config_(config),
          rank_(rank),
          io_(config),
          min_degree_(config_.min_degree),
          total_degree_(2 * config_.min_degree) {
        PEID size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        // Init variables
        from_ = rank * std::ceil(config_.n / (LPFloat)size);
        to_   = std::min((SInt)((rank + 1) * ceil(config_.n / (LPFloat)size) - 1), config_.n - 1);
        std::cout << "f " << from_ << " t " << to_ << std::endl;
    }

    void Generate() {
        GenerateEdges();

        // Additional stats
        if (rank_ == ROOT) {
            HPFloat terra_bytes = 128 * config_.n * min_degree_;
            terra_bytes /= (8);
            terra_bytes /= (1024);
            terra_bytes /= (1024);
            terra_bytes /= (1024);
            terra_bytes /= (1024);

            std::cout << "memory of graph in tera bytes  " << std::setprecision(40) << terra_bytes << std::endl;
        }
    }

    GeneratorIO& IO() {
        return io_;
    }

    std::pair<SInt, SInt> GetVertexRange() {
        return std::make_pair(from_, to_);
    }

private:
    // Config
    PGeneratorConfig& config_;
    PEID              rank_;

    // I/O
    GeneratorIO io_;

    // Constants and variables
    SInt min_degree_;
    SInt total_degree_;
    SInt from_, to_;

    void GenerateEdges() {
        for (SInt v = from_; v <= to_; v++) {
            for (SInt i = 0; i < min_degree_; i++) {
                SInt r = 2 * (v * min_degree_ + i) + 1;
                do {
                    // compute hash h(r)
                    SInt hash = sampling::Spooky::hash(config_.seed + r);
                    r         = hash % r;
                } while (r % 2 == 1);
                SInt w = r / total_degree_;
                io_.PushEdge(v, w);
            }
        }
    };
};
} // namespace kagen
