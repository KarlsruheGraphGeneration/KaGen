#include "kagen_library.h"

#include "generator_config.h"

#include "barabassi/barabassi.h"
#include "geometric/rgg/rgg_2d.h"
#include "geometric/rgg/rgg_3d.h"
#include "gnm/gnm_directed.h"
#include "gnm/gnm_undirected.h"
#include "gnp/gnp_directed.h"
#include "gnp/gnp_undirected.h"
#include "grid/grid_2d.h"
#include "hyperbolic/hyperbolic.h"
#include "kronecker/kronecker.h"

#ifdef KAGEN_CGAL_FOUND
#include "geometric/delaunay/delaunay_2d.h"
#include "geometric/delaunay/delaunay_3d.h"
#endif 

#include <iostream>

namespace kagen {
KaGen::KaGen(const PEID rank, const PEID size)
    : rank_(rank), size_(size), config_(std::make_unique<PGeneratorConfig>()) {
  SetDefaults();
}

KaGen::~KaGen() = default;
  

KaGenResult KaGen::GenerateDirectedGMM(const SInt n, const SInt m, const SInt k, const int seed,
                                       const std::string &output, const bool self_loops) {
  // Update config
  config_->n = n;
  config_->m = m;
  config_->k = (k == 0 ? config_->k : k);
  config_->seed= seed;
  config_->output_file = output;
  config_->self_loops = self_loops;

  // Edge callback
  EdgeList edges;
  auto edge_cb = [&](SInt source, SInt target) { edges.emplace_back(source, target); };

  // Init and run generator
  GNMDirected<decltype(edge_cb)> gen(*config_, rank_, edge_cb);
  gen.Generate();

  return {std::move(edges), gen.GetVertexRange()};
}

KaGenResult KaGen::GenerateUndirectedGNM(const SInt n, const SInt m, const SInt k, const int seed,
                                         const std::string &output, const bool self_loops) {
  // Update config
  config_->n = n;
  config_->m = m;
  config_->k = (k == 0 ? config_->k : k);
  config_->seed = seed;
  config_->output_file = output;
  config_->self_loops = self_loops;

  // Edge callback
  EdgeList edges;
  auto edge_cb = [&](SInt source, SInt target) { edges.emplace_back(source, target); };

  // Init and run generator
  GNMUndirected<decltype(edge_cb)> gen(*config_, rank_, edge_cb);
  gen.Generate();

  return {std::move(edges), gen.GetVertexRange()};
}

KaGenResult KaGen::GenerateDirectedGNP(const SInt n, const LPFloat p, const SInt k, const int seed,
                                       const std::string &output, const bool self_loops) {
  // Update config
  config_->n = n;
  config_->p = p;
  config_->k = (k == 0 ? config_->k : k);
  config_->seed = seed;
  config_->output_file = output;
  config_->self_loops = self_loops;

  // Edge callback
  EdgeList edges;
  auto edge_cb = [&](SInt source, SInt target) { edges.emplace_back(source, target); };

  // Init and run generator
  GNPDirected<decltype(edge_cb)> gen(*config_, rank_, edge_cb);
  gen.Generate();

  return {std::move(edges), gen.GetVertexRange()};
}

KaGenResult KaGen::GenerateUndirectedGNP(const SInt n, const LPFloat p, const SInt k, const int seed,
                                         const std::string &output, const bool self_loops) {
  // Update config
  config_->n = n;
  config_->p = p;
  config_->k = (k == 0 ? config_->k : k);
  config_->seed = seed;
  config_->output_file = output;
  config_->self_loops = self_loops;

  // Edge callback
  EdgeList edges;
  auto edge_cb = [&](SInt source, SInt target) { edges.emplace_back(source, target); };

  // Init and run generator
  GNPUndirected<decltype(edge_cb)> gen(*config_, rank_, edge_cb);
  gen.Generate();

  return {std::move(edges), gen.GetVertexRange()};
}

KaGenResult KaGen::Generate2DRGG(const SInt n, const LPFloat r, const SInt k, const int seed,
                                 const std::string &output) {
  // Update config
  config_->n = n;
  config_->r = r;
  config_->k = (k == 0 ? config_->k : k);
  config_->seed = seed;
  config_->output_file = output;

  // Edge callback
  EdgeList edges;
  auto edge_cb = [&](SInt source, SInt target) { edges.emplace_back(source, target); };

  // Init and run generator
  RGG2D<decltype(edge_cb)> gen(*config_, rank_, edge_cb);
  gen.Generate();

  std::cout << "Range for PE " << rank_ << ": " << gen.GetVertexRange().first << " -- " << gen.GetVertexRange().second << std::endl;
  
  return {std::move(edges), gen.GetVertexRange()};
}

KaGenResult KaGen::Generate3DRGG(const SInt n, const LPFloat r, const SInt k, const int seed,
                                 const std::string &output) {
  // Update config
  config_->n = n;
  config_->r = r;
  config_->k = (k == 0 ? config_->k : k);
  config_->seed = seed;
  config_->output_file = output;

  // Edge callback
  EdgeList edges;
  auto edge_cb = [&](SInt source, SInt target) { edges.emplace_back(source, target); };

  // Init and run generator
  RGG3D<decltype(edge_cb)> gen(*config_, rank_, edge_cb);
  gen.Generate();

  return {std::move(edges), gen.GetVertexRange()};
}

#ifdef KAGEN_CGAL_FOUND
KaGenResult KaGen::Generate2DRDG(const SInt n, const SInt k, const int seed, const std::string &output) {
  // Update config
  config_->n = n;
  config_->k = (k == 0 ? config_->k : k);
  config_->seed = seed;
  config_->output_file = output;

  // Edge callback
  EdgeList edges;
  auto edge_cb = [&](SInt source, SInt target) { edges.emplace_back(source, target); };

  // Init and run generator
  Delaunay2D<decltype(edge_cb)> gen(*config_, rank_, edge_cb);
  gen.Generate();

  return {std::move(edges), gen.GetVertexRange()};
}

KaGenResult KaGen::Generate3DRDG(const SInt n, const SInt k, const int seed, const std::string &output) {
  // Update config
  config_->n = n;
  config_->k = (k == 0 ? config_->k : k);
  config_->seed = seed;
  config_->output_file = output;

  // Edge callback
  EdgeList edges;
  auto edge_cb = [&](SInt source, SInt target) { edges.emplace_back(source, target); };

  // Init and run generator
  Delaunay3D<decltype(edge_cb)> gen(*config_, rank_, edge_cb);
  gen.Generate();

  return {std::move(edges), gen.GetVertexRange()};
}
#else // KAGEN_CGAL_FOUND
KaGenResult KaGen::Generate2DRDG(SInt, SInt, int, const std::string &) {
  throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}

KaGenResult KaGen::Generate3DRDG(SInt, SInt, int, const std::string &) {
  throw std::runtime_error("Library was compiled without CGAL. Thus, delaunay generators are not available.");
}
#endif // KAGEN_CGAL_FOUND

KaGenResult KaGen::GenerateBA(const SInt n, const SInt d, const SInt k, const int seed, const std::string &output) {
  // Update config
  config_->n = n;
  config_->min_degree = d;
  config_->k = (k == 0 ? config_->k : k);
  config_->seed = seed;
  config_->output_file = output;

  // Edge callback
  EdgeList edges;
  auto edge_cb = [&](SInt source, SInt target) { edges.emplace_back(source, target); };

  // Init and run generator
  Barabassi<decltype(edge_cb)> gen(*config_, rank_, edge_cb);
  gen.Generate();

  return {std::move(edges), gen.GetVertexRange()};
}

KaGenResult KaGen::GenerateRHG(const SInt n, const LPFloat gamma, const SInt d, const SInt k, const int seed,
                               const std::string &output) {
  // Update config
  config_->n = n;
  config_->plexp = gamma;
  config_->avg_degree = d;
  config_->k = (k == 0 ? config_->k : k);
  config_->seed = seed;
  config_->output_file = output;

  // Edge callback
  EdgeList edges;
  auto edge_cb = [&](SInt source, SInt target) { edges.emplace_back(source, target); };


  // Init and run generator
  Hyperbolic<decltype(edge_cb)> gen(*config_, rank_, edge_cb);
  gen.Generate();

  return {std::move(edges), gen.GetVertexRange()};
}

KaGenResult KaGen::Generate2DGrid(const SInt n, const SInt m, const LPFloat p, const SInt periodic, const SInt k,
                                  const int seed, const std::string &output) {
  // Update config
  config_->n = n;
  config_->m = m;
  config_->p = p;
  config_->periodic = periodic;
  config_->k = (k == 0 ? config_->k : k);
  config_->seed = seed;
  config_->output_file = output;

  // Edge callback
  EdgeList edges;
  auto edge_cb = [&](SInt source, SInt target) { edges.emplace_back(source, target); };

  // Init and run generator
  Grid2D<decltype(edge_cb)> gen(*config_, rank_, edge_cb);
  gen.Generate();

  return {std::move(edges), gen.GetVertexRange()};
}

KaGenResult KaGen::GenerateKronecker(const SInt n, const SInt m, const SInt k,
                                  const int seed, const std::string &output) {
  // Update config
  config_->n = n;
  config_->m = m;
  config_->k = (k == 0 ? config_->k : k);
  config_->seed = seed;
  config_->output_file = output;

  // Edge callback
  EdgeList edges;
  auto edge_cb = [&](SInt source, SInt target) { edges.emplace_back(source, target); };

  // Init and run generator
  Kronecker<decltype(edge_cb)> gen(*config_, rank_, edge_cb);
  gen.Generate();

  return {std::move(edges), gen.GetVertexRange()};
}

void KaGen::SetDefaults() {
  config_->n = 100;
  config_->m = 0;
  config_->k = size_;
  config_->seed = 1;
  config_->hash_sample = false;
  config_->use_binom = false;
  config_->output_file = "out";
  config_->debug_output = "dbg";
  config_->dist_size = 10;
  config_->p = 0.0;
  config_->self_loops = false;
  config_->r = 0.125;
  config_->avg_degree = 5.0;
  config_->plexp = 2.6;
  config_->thres = 0;
  config_->query_both = true;
  config_->min_degree = 4;
  config_->precision = 32;
  config_->base_size = (SInt)1 << 8;
  config_->hyp_base = (SInt)1 << 8;
  config_->iterations = 1;
}
} // namespace kagen 
