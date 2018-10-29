/*******************************************************************************
 * interface/kagen_interface.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _KAGEN_INTERFACE_H_
#define _KAGEN_INTERFACE_H_

#include <iostream>
#include <mpi.h>

#include "definitions.h"
#include "generator_config.h"
#include "parse_parameters.h"

#include "geometric/delaunay/delaunay_2d.h"
#include "geometric/delaunay/delaunay_3d.h"
#include "geometric/rgg/rgg_2d.h"
#include "geometric/rgg/rgg_3d.h"
#include "gnm/gnm_directed.h"
#include "gnm/gnm_undirected.h"
#include "gnp/gnp_directed.h"
#include "gnp/gnp_undirected.h"
#include "hyperbolic/hyperbolic.h"
#include "barabassi/barabassi.h"

namespace kagen {

typedef std::vector<std::pair<SInt, SInt>> EdgeList;

class KaGen {
 public:
  KaGen(const PEID rank, const PEID size) 
    : rank_(rank), size_(size) { 
      SetDefaults();
    }

  virtual ~KaGen() = default;

  EdgeList GenerateDirectedGNM(SInt n, 
                               SInt m, 
                               SInt k = 0, 
                               SInt seed = 1, 
                               const std::string &output = "out", 
                               bool self_loops = false) {
    EdgeList edges; 

    // Update config
    config_.n = n;
    config_.m = m;
    config_.k = (k == 0 ? config_.k : k);
    config_.seed = seed;
    config_.output_file = output;
    config_.self_loops = self_loops;

    // Edge callback
    auto edge_cb = [&](SInt source, SInt target) {
      edges.emplace_back(source, target);
    };

    // Init and run generator
    GNMDirected<decltype(edge_cb)> gen(config_, rank_, edge_cb);
    gen.Generate();

    edges.insert(begin(edges), gen.GetVertexRange());
    return edges;
  }

  EdgeList GenerateUndirectedGNM(SInt n, 
                                 SInt m, 
                                 SInt k = 0, 
                                 SInt seed = 1, 
                                 const std::string &output = "out", 
                                 bool self_loops = false) {
    EdgeList edges; 

    // Update config
    config_.n = n;
    config_.m = m;
    config_.k = (k == 0 ? config_.k : k);
    config_.seed = seed;
    config_.output_file = output;
    config_.self_loops = self_loops;

    // Edge callback
    auto edge_cb = [&](SInt source, SInt target) {
      edges.emplace_back(source, target);
    };

    // Init and run generator
    GNMUndirected<decltype(edge_cb)> gen(config_, rank_, edge_cb);
    gen.Generate();

    edges.insert(begin(edges), gen.GetVertexRange());
    return edges;
  }

  EdgeList GenerateDirectedGNP(SInt n, 
                               LPFloat p, 
                               SInt k = 0, 
                               SInt seed = 1, 
                               const std::string &output = "out", 
                               bool self_loops = false) {
    EdgeList edges; 

    // Update config
    config_.n = n;
    config_.p = p;
    config_.k = (k == 0 ? config_.k : k);
    config_.seed = seed;
    config_.output_file = output;
    config_.self_loops = self_loops;

    // Edge callback
    auto edge_cb = [&](SInt source, SInt target) {
      edges.emplace_back(source, target);
    };

    // Init and run generator
    GNPDirected<decltype(edge_cb)> gen(config_, rank_, edge_cb);
    gen.Generate();

    edges.insert(begin(edges), gen.GetVertexRange());
    return edges;
  }

  EdgeList GenerateUndirectedGNP(SInt n, 
                                 LPFloat p, 
                                 SInt k = 0, 
                                 SInt seed = 1, 
                                 const std::string &output = "out", 
                                 bool self_loops = false) {
    EdgeList edges; 

    // Update config
    config_.n = n;
    config_.p = p;
    config_.k = (k == 0 ? config_.k : k);
    config_.seed = seed;
    config_.output_file = output;
    config_.self_loops = self_loops;

    // Edge callback
    auto edge_cb = [&](SInt source, SInt target) {
      edges.emplace_back(source, target);
    };

    // Init and run generator
    GNPUndirected<decltype(edge_cb)> gen(config_, rank_, edge_cb);
    gen.Generate();

    edges.insert(begin(edges), gen.GetVertexRange());
    return edges;
  }

  EdgeList Generate2DRGG(SInt n, 
                         LPFloat r, 
                         SInt k = 0, 
                         SInt seed = 1, 
                         const std::string &output = "out") {
    EdgeList edges; 

    // Update config
    config_.n = n;
    config_.r = r;
    config_.k = (k == 0 ? config_.k : k);
    config_.seed = seed;
    config_.output_file = output;

    // Edge callback
    auto edge_cb = [&](SInt source, SInt target) {
      edges.emplace_back(source, target);
    };

    // Init and run generator
    RGG2D<decltype(edge_cb)> gen(config_, rank_, edge_cb);
    gen.Generate();

    edges.insert(begin(edges), gen.GetVertexRange());
    return edges;
  }

  EdgeList Generate3DRGG(SInt n, 
                         LPFloat r, 
                         SInt k = 0, 
                         SInt seed = 1, 
                         const std::string &output = "out") {
    EdgeList edges; 

    // Update config
    config_.n = n;
    config_.r = r;
    config_.k = (k == 0 ? config_.k : k);
    config_.seed = seed;
    config_.output_file = output;

    // Edge callback
    auto edge_cb = [&](SInt source, SInt target) {
      edges.emplace_back(source, target);
    };

    // Init and run generator
    RGG3D<decltype(edge_cb)> gen(config_, rank_, edge_cb);
    gen.Generate();

    edges.insert(begin(edges), gen.GetVertexRange());
    return edges;
  }

  EdgeList Generate2DRDG(SInt n, 
                         SInt k = 0, 
                         SInt seed = 1, 
                         const std::string &output = "out") {
    EdgeList edges; 

    // Update config
    config_.n = n;
    config_.k = (k == 0 ? config_.k : k);
    config_.seed = seed;
    config_.output_file = output;

    // Edge callback
    auto edge_cb = [&](SInt source, SInt target) {
      edges.emplace_back(source, target);
    };

    // Init and run generator
    Delaunay2D<decltype(edge_cb)> gen(config_, rank_, edge_cb);
    gen.Generate();

    edges.insert(begin(edges), gen.GetVertexRange());
    return edges;
  }

  EdgeList Generate3DRDG(SInt n, 
                         SInt k = 0, 
                         SInt seed = 1, 
                         const std::string &output = "out") {
    EdgeList edges; 

    // Update config
    config_.n = n;
    config_.k = (k == 0 ? config_.k : k);
    config_.seed = seed;
    config_.output_file = output;

    // Edge callback
    auto edge_cb = [&](SInt source, SInt target) {
      edges.emplace_back(source, target);
    };

    // Init and run generator
    Delaunay3D<decltype(edge_cb)> gen(config_, rank_, edge_cb);
    gen.Generate();

    edges.insert(begin(edges), gen.GetVertexRange());
    return edges;
  }

  EdgeList GenerateBA(SInt n, 
                      SInt d,
                      SInt k = 0, 
                      SInt seed = 1, 
                      const std::string &output = "out") {
    EdgeList edges; 

    // Update config
    config_.n = n;
    config_.min_degree = d;
    config_.k = (k == 0 ? config_.k : k);
    config_.seed = seed;
    config_.output_file = output;

    // Edge callback
    auto edge_cb = [&](SInt source, SInt target) {
      edges.emplace_back(source, target);
    };

    // Init and run generator
    Barabassi<decltype(edge_cb)> gen(config_, rank_, edge_cb);
    gen.Generate();

    edges.insert(begin(edges), gen.GetVertexRange());
    return edges;
  }

  EdgeList GenerateRHG(SInt n, 
                       LPFloat gamma,
                       SInt d,
                       SInt k = 0, 
                       SInt seed = 1, 
                       const std::string &output = "out") {
    EdgeList edges; 

    // Update config
    config_.n = n;
    config_.plexp = gamma;
    config_.avg_degree = d;
    config_.k = (k == 0 ? config_.k : k);
    config_.seed = seed;
    config_.output_file = output;

    // Edge callback
    auto edge_cb = [&](SInt source, SInt target) {
      edges.emplace_back(source, target);
    };

    // Init and run generator
    Hyperbolic<decltype(edge_cb)> gen(config_, rank_, edge_cb);
    gen.Generate();

    edges.insert(begin(edges), gen.GetVertexRange());
    return edges;
  }

 private:
  // PE status
  PEID rank_, size_;

  PGeneratorConfig config_; 

  void SetDefaults() {
    config_.n = 100;
    config_.m = 0;
    config_.k = size_;
    config_.seed = 1;
    config_.hash_sample = false;
    config_.use_binom = false;
    config_.output_file = "out";
    config_.debug_output = "dbg";
    config_.dist_size = 10;
    config_.p = 0.0;
    config_.self_loops = false;
    config_.r = 0.125;
    config_.avg_degree = 5.0;
    config_.plexp = 2.6;
    config_.thres = 0;
    config_.query_both = false;
    config_.min_degree = 4;
    config_.precision = 32;
    config_.base_size = (ULONG)1 << 8;
    config_.hyp_base = (ULONG)1 << 8;
    config_.iterations = 1;
  }
};

}
#endif
