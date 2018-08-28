/*******************************************************************************
 * include/generators/barabassi/barabassi.h
 *
 * Copyright (C) 2016-2017 Christian Schulz <christian.schulz@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _BARABASSI_H_
#define _BARABASSI_H_

#include <iostream>
#include <vector>

#include "definitions.h"
#include "generator_config.h"
#include "generator_io.h"
#include "hash.hpp"

namespace kagen {

template <typename EdgeCallback> 
class Barabassi {
 public:
  Barabassi(const PGeneratorConfig &config, const PEID rank,
            const EdgeCallback &cb)
      : config_(config),
        rank_(rank),
        io_(config),
        cb_(cb),
        min_degree_(config_.min_degree),
        total_degree_(2 * config_.min_degree) {
    PEID size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Init variables
    from_ = rank * ceil(config_.n / (LPFloat)size);
    to_ = std::min((SInt)((rank + 1) * ceil(config_.n / (LPFloat)size) - 1),
                   config_.n - 1);
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

      std::cout << "memory of graph in tera bytes  " << std::setprecision(40)
                << terra_bytes << std::endl;
    }
  }

  void Output() const {
#ifdef OUTPUT_EDGES
    io_.OutputEdges();
#else
    io_.OutputDist();
#endif
  }

  std::pair<SInt, SInt> GetVertexRange() {
    return std::make_pair(from_, to_);
  }

  SInt NumberOfEdges() const { return io_.NumEdges(); }

 private:
  // Config
  PGeneratorConfig config_;
  PEID rank_;

  // I/O
  GeneratorIO<> io_;
  EdgeCallback cb_; 

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
          r = hash % r;
        } while (r % 2 == 1);
        SInt w = r / total_degree_;
        cb_(v, w);
#ifdef OUTPUT_EDGES
        io_.PushEdge(v, w);
#else
        io_.UpdateDist(v);
        io_.UpdateDist(w);
#endif
      }
    }
  };
};

}
#endif
