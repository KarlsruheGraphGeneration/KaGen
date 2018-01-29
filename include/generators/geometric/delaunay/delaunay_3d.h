/*******************************************************************************
 * include/generators/geometric/delaunay/delaunay_3d.h
 *
 * Copyright (C) 2017 Daniel Funke <funke@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _DELAUNAY_3D_H_
#define _DELAUNAY_3D_H_

#include "geometric/geometric_3d.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/spatial_sort.h>

#include <sys/stat.h>

using K_3d = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb_3d =
    CGAL::Triangulation_vertex_base_with_info_3<SInt, K_3d>; // attach an ID to
                                                             // each point
using Fb_3d = CGAL::Triangulation_cell_base_3<K_3d>;
using Tds_3d = CGAL::Triangulation_data_structure_3<Vb_3d, Fb_3d>;
using Dt_3d = CGAL::Delaunay_triangulation_3<K_3d, Tds_3d>;
using Point_3d = Dt_3d::Point;
using Fh_3d = Dt_3d::Cell_handle;
using Vh_3d = Dt_3d::Vertex_handle;

using Points_3d = std::vector<std::pair<Point_3d, SInt>>; // needed to add
                                                          // points with ID to
                                                          // CGAL

using Circ_3d = CGAL::Sphere_3<K_3d>;
using Box_3d = CGAL::Bbox_3;

class Delaunay3D : public Geometric3D {
public:
  Delaunay3D(const PGeneratorConfig &config, const PEID rank)
      : Geometric3D(config, rank), point_io_(config), edge_io_(config) {
    // Chunk variables
    total_chunks_ = config_.k;
    chunks_per_dim_ = cbrt(config_.k);
    chunk_size_ = 1.0 / chunks_per_dim_;

    // Cell variables

    // upper bound for expected distance to nearest neighbor, according to paper
    LPFloat edNN = 4.8 / cbrt(config_.n);
    LPFloat lCellSize = std::min(edNN, chunk_size_);

    // use chunks_per_dim_ as upper bound for cells_per_dim_
    cells_per_dim_ = std::min(SInt(48), (SInt)floor(chunk_size_ / lCellSize));
    cells_per_chunk_ = cells_per_dim_ * cells_per_dim_ * cells_per_dim_;
    cell_size_ = chunk_size_ / cells_per_dim_;

    max_radius_ = (cells_per_dim_ * chunks_per_dim_);

    InitDatastructures();
  }

  void Output() const override {
#ifdef DEL_STATS
    outputRadiusStats();
#endif
#ifdef OUTPUT_EDGES
    edge_io_.OutputEdges();
#else
    edge_io_.OutputDist();
#endif
    // point_io_.OutputEdges();
    // edge_io_.OutputEdges();
    // adjacency_io_.OutputEdges();
  }

  SInt NumberOfEdges() const override { return edge_io_.NumEdges(); }

private:
  // I/O
  GeneratorIO<std::tuple<LPFloat, LPFloat, LPFloat>> point_io_;
  GeneratorIO<std::tuple<SInt, SInt>> edge_io_;

#ifdef DEL_STATS

  struct del_stats {
    del_stats(int r, SInt p_cell, SInt p_tria, SInt s)
        : radius(r), n_p_cell(p_cell), n_p_tria(p_tria), n_simplices(s) {}

    int radius;
    SInt n_p_cell;
    SInt n_p_tria;
    SInt n_simplices;
  };

  std::map<SInt, del_stats> radius_stats_;

  void outputRadiusStats() const {

    std::string fileName = config_.output_file + "_stats_" + std::to_string(rank_);

    FILE *fout = fopen(fileName.c_str(), "a");

//    fprintf(fout, "STATS begin\n");
//    fprintf(fout, "%llu %llu %llu %f %llu %llu %f %llu\n", config_.n,
//            total_chunks_, chunks_per_dim_, chunk_size_, cells_per_chunk_,
//            cells_per_dim_, cell_size_, max_radius_);

    for (const auto &s : radius_stats_) {
//      SInt x, y, z;
//      Decode(s.first, x, y, z);
      fprintf(fout, "%u %llu %llu %llu %f %llu %llu %f %llu %llu %i %llu %llu %llu\n", 3, config_.n,
                total_chunks_, chunks_per_dim_, chunk_size_, cells_per_chunk_,
                cells_per_dim_, cell_size_, max_radius_, s.first,
                s.second.radius, s.second.n_p_cell, s.second.n_p_tria,
                s.second.n_simplices);
    }

//    fprintf(fout, "STATS end\n");

    fclose(fout);
  };
#endif

  SInt max_radius_;
  static constexpr SInt COPY_FLAG = SInt(1) << (sizeof(SInt) * CHAR_BIT - 1);

  void GenerateEdges(const SInt chunk_row, const SInt chunk_column,
                     const SInt chunk_depth) override {
    SInt chunk_id = Encode(chunk_column, chunk_row, chunk_depth);
    SInt id_low = std::get<5>(chunks_[chunk_id]);
    SInt id_high = id_low + std::get<0>(chunks_[chunk_id]);
    //            CGALCairoPainter basePainter;

    Dt_3d tria;
    Fh_3d hint;

    // bounding box of own chunk
    Box_3d bbChunk(chunk_row * chunk_size_, chunk_column * chunk_size_,
                   chunk_depth * chunk_size_, (chunk_row + 1) * chunk_size_,
                   (chunk_column + 1) * chunk_size_,
                   (chunk_depth + 1) * chunk_size_);

    //            printf("[%llu] chunk (%llu, %llu, %llu): (%f, %f, %f) - (%f,
    //            %f, %f)\n",
    //                   chunk_id, chunk_row, chunk_column, chunk_depth,
    //                   bbChunk.xmin(), bbChunk.ymin(), bbChunk.zmin(),
    //                   bbChunk.xmax(), bbChunk.ymax(), bbChunk.zmax());

    // Iterate chunk cells
    for (SInt cell_row = 0; cell_row < cells_per_dim_; ++cell_row) {
      for (SInt cell_column = 0; cell_column < cells_per_dim_; ++cell_column) {
        for (SInt cell_depth = 0; cell_depth < cells_per_dim_; ++cell_depth) {
          SInt cell_id = cell_row * cells_per_dim_ + cell_column +
                         (cells_per_dim_ * cells_per_dim_) * cell_depth;
          SInt cellIdx = ComputeGlobalCellId(chunk_id, cell_id);

          //                    printf("[%llu] adding %lu points from own cell
          //                    %llu (%llu, %llu)\n",
          //                           chunk_id, vertices_[cellIdx].size(),
          //                           cell_id, cell_row, cell_column);

        SortCellVertices(vertices_[cellIdx]);
          for (const auto &v : vertices_[cellIdx]) {
            Point_3d p(std::get<0>(v), std::get<1>(v), std::get<2>(v));
            assert(bbChunk.xmin() <= p.x() && p.x() <= bbChunk.xmax());
            assert(bbChunk.ymin() <= p.y() && p.y() <= bbChunk.ymax());
            assert(bbChunk.zmin() <= p.z() && p.z() <= bbChunk.zmax());

            auto vh = tria.insert(p, hint);
            vh->info() = std::get<3>(v);
            hint = vh->cell();
          }
        }
      }
    }

    //            printf("[%llu] %lu own points\n",
    //                   chunk_id, points.size());

    bool conflictFree = false;
    SInt cell_radius = 1;
    for (; !conflictFree && cell_radius <= max_radius_; ++cell_radius) {
      // bounding box of neighborhood
      Box_3d bbNH(chunk_row * chunk_size_ - cell_radius * cell_size_,
                  chunk_column * chunk_size_ - cell_radius * cell_size_,
                  chunk_depth * chunk_size_ - cell_radius * cell_size_,
                  (chunk_row + 1) * chunk_size_ + cell_radius * cell_size_,
                  (chunk_column + 1) * chunk_size_ + cell_radius * cell_size_,
                  (chunk_depth + 1) * chunk_size_ + cell_radius * cell_size_);

      //        printf("[%llu] chunk (%llu, %llu, %llu), radius %llu,
      //        neighborhood: (%f, %f, %f) - (%f, %f, %f)\n",
      //               chunk_id, chunk_row, chunk_column, chunk_depth,
      //               cell_radius, bbNH.xmin(), bbNH.ymin(), bbNH.zmin(),
      //               bbNH.xmax(), bbNH.ymax(), bbNH.zmax());

      SInt chunk_radius =
          (SInt)std::ceil(cell_radius / static_cast<LPFloat>((cells_per_dim_)));
      for (SSInt chunk_row_diff = -chunk_radius;
           chunk_row_diff <= (SSInt)chunk_radius; ++chunk_row_diff) {
        for (SSInt chunk_col_diff = -chunk_radius;
             chunk_col_diff <= (SSInt)chunk_radius; ++chunk_col_diff) {
          for (SSInt chunk_dep_diff = -chunk_radius;
               chunk_dep_diff <= (SSInt)chunk_radius; ++chunk_dep_diff) {
            if ((SInt)std::max({std::abs(chunk_row_diff),
                                std::abs(chunk_col_diff),
                                std::abs(chunk_dep_diff)}) < chunk_radius)
              continue; // skip chunks in inner region (already added)

            SSInt neighbor_chunk_col = chunk_column + chunk_col_diff;
            SSInt neighbor_chunk_row = chunk_row + chunk_row_diff;
            SSInt neighbor_chunk_dep = chunk_depth + chunk_dep_diff;

            LPFloat x_offset = 0;
            if (neighbor_chunk_row < 0) {
              x_offset =
                  -1.0 * (std::ceil(std::abs(neighbor_chunk_row) /
                                    static_cast<LPFloat>(chunks_per_dim_)));
            } else {
              x_offset =
                  1.0 * (std::floor(std::abs(neighbor_chunk_row) /
                                    static_cast<LPFloat>(chunks_per_dim_)));
            }

            LPFloat y_offset = 0;
            if (neighbor_chunk_col < 0) {
              y_offset =
                  -1.0 * (std::ceil(std::abs(neighbor_chunk_col) /
                                    static_cast<LPFloat>(chunks_per_dim_)));
            } else {
              y_offset =
                  1.0 * (std::floor(std::abs(neighbor_chunk_col) /
                                    static_cast<LPFloat>(chunks_per_dim_)));
            }

            LPFloat z_offset = 0;
            if (neighbor_chunk_dep < 0) {
              z_offset =
                  -1.0 * (std::ceil(std::abs(neighbor_chunk_dep) /
                                    static_cast<LPFloat>(chunks_per_dim_)));
            } else {
              z_offset =
                  1.0 * (std::floor(std::abs(neighbor_chunk_dep) /
                                    static_cast<LPFloat>(chunks_per_dim_)));
            }

            // Get correct grid cells
            neighbor_chunk_row = (neighbor_chunk_row % (SSInt)chunks_per_dim_ +
                                  chunks_per_dim_) %
                                 chunks_per_dim_;
            neighbor_chunk_col = (neighbor_chunk_col % (SSInt)chunks_per_dim_ +
                                  chunks_per_dim_) %
                                 chunks_per_dim_;
            neighbor_chunk_dep = (neighbor_chunk_dep % (SSInt)chunks_per_dim_ +
                                  chunks_per_dim_) %
                                 chunks_per_dim_;
            SInt neighbor_chunk_id =
                Encode((SInt)neighbor_chunk_col, (SInt)neighbor_chunk_row,
                       (SInt)neighbor_chunk_dep);

            // now we have our neighbor chunk -> dive into its cells
            SInt row_cell_diff = std::min(
                (cell_radius - (std::abs(chunk_row_diff) - 1) * cells_per_dim_),
                cells_per_dim_);
            SInt col_cell_diff = std::min(
                (cell_radius - (std::abs(chunk_col_diff) - 1) * cells_per_dim_),
                cells_per_dim_);
            SInt dep_cell_diff = std::min(
                (cell_radius - (std::abs(chunk_dep_diff) - 1) * cells_per_dim_),
                cells_per_dim_);

            for (SInt dcr = 0; dcr < row_cell_diff; ++dcr) {
              for (SInt dcc = 0; dcc < col_cell_diff; ++dcc) {
                for (SInt dcd = 0; dcd < dep_cell_diff; ++dcd) {
                  SInt neighbor_cell_row =
                      (chunk_row_diff < 0 ? cells_per_dim_ - row_cell_diff
                                          : 0) +
                      dcr;
                  SInt neighbor_cell_col =
                      (chunk_col_diff < 0 ? cells_per_dim_ - col_cell_diff
                                          : 0) +
                      dcc;
                  SInt neighbor_cell_dep =
                      (chunk_dep_diff < 0 ? cells_per_dim_ - dep_cell_diff
                                          : 0) +
                      dcd;

                  SInt total_row_diff =
                      chunk_row_diff != 0
                          ? (chunk_row_diff < 0
                                 ? cells_per_dim_ - neighbor_cell_row
                                 : neighbor_cell_row + 1) +
                                (std::abs(chunk_row_diff) - 1) * cells_per_dim_
                          : 0;
                  SInt total_col_diff =
                      chunk_col_diff != 0
                          ? (chunk_col_diff < 0
                                 ? cells_per_dim_ - neighbor_cell_col
                                 : neighbor_cell_col + 1) +
                                (std::abs(chunk_col_diff) - 1) * cells_per_dim_
                          : 0;
                  SInt total_dep_diff =
                      chunk_dep_diff != 0
                          ? (chunk_dep_diff < 0
                                 ? cells_per_dim_ - neighbor_cell_dep
                                 : neighbor_cell_dep + 1) +
                                (std::abs(chunk_dep_diff) - 1) * cells_per_dim_
                          : 0;

                  if (std::max({total_row_diff, total_col_diff,
                                total_dep_diff}) < cell_radius)
                    continue; // skip cell in inner region (already added)

                  SInt neighbor_cell_id =
                      neighbor_cell_row * cells_per_dim_ + neighbor_cell_col +
                      (cells_per_dim_ * cells_per_dim_) * neighbor_cell_dep;

                  // Check if vertices not generated
                  SInt neighbor_cell_offset =
                      ComputeGlobalCellId(neighbor_chunk_id, neighbor_cell_id);
                  // lazily generate vertices
                  GenerateVertices(neighbor_chunk_id, neighbor_cell_id);

                  // Gather vertices
                  //                            printf("[%llu] adding %lu points
                  //                            from chunk %llu (%llu,%llu,%llu)
                  //                            cell %llu (%llu, %llu, %llu)
                  //                            with offset (%f,%f,%f)\n",
                  //                                   chunk_id,
                  //                                   vertices_[neighbor_cell_offset].size(),
                  //                                   neighbor_chunk_id,
                  //                                   neighbor_chunk_row,
                  //                                   neighbor_chunk_col,
                  //                                   neighbor_chunk_dep,
                  //                                   neighbor_cell_id,
                  //                                   neighbor_cell_row,
                  //                                   neighbor_cell_col,
                  //                                   neighbor_cell_dep,
                  //                                   x_offset, y_offset,
                  //                                   z_offset);

                  SortCellVertices(vertices_[neighbor_cell_offset]);
                  for (const auto &v : vertices_[neighbor_cell_offset]) {
                    Point_3d p(std::get<0>(v) + x_offset,
                               std::get<1>(v) + y_offset,
                               std::get<2>(v) + z_offset);
                    assert(bbNH.xmin() <= p.x() && p.x() <= bbNH.xmax());
                    assert(bbNH.ymin() <= p.y() && p.y() <= bbNH.ymax());
                    assert(bbNH.zmin() <= p.z() && p.z() <= bbNH.zmax());

                    auto vh = tria.insert(p, hint);
                    vh->info() = std::get<3>(v) + COPY_FLAG;
                    hint = vh->cell();
                  }
                }
              }
            }
          }
        }
      }
      // printf("[%llu] level %llu %lu total points\n",
      //       chunk_id, cell_radius, points.size());

      // printf("[%llu] Triangulation of %lu points with %lu faces\n",
      //       chunk_id, tria.number_of_vertices(), tria.number_of_cells());

      // painter.drawTriangles(tria);

      conflictFree = true;
      for (auto f = tria.finite_cells_begin(); f != tria.finite_cells_end(); ++f) {
          bool touchesChunk =    (!(f->vertex(0)->info() & COPY_FLAG))
                              || (!(f->vertex(1)->info() & COPY_FLAG))
                              || (!(f->vertex(2)->info() & COPY_FLAG))
                              || (!(f->vertex(3)->info() & COPY_FLAG));

        if (touchesChunk) {
          Circ_3d c(f->vertex(0)->point(), f->vertex(1)->point(),
                    f->vertex(2)->point(), f->vertex(3)->point());

          if (!boxContains(bbNH, c)) {
            //                printf("simplex (%f, %f, %f) (%f, %f, %f) (%f, %f,
            //                %f) (%f, %f, %f), circumsphere (%f, %f, %f) -
            //                %f\n",
            //                       f->vertex(0)->point().x(),
            //                       f->vertex(0)->point().y(),
            //                       f->vertex(0)->point().z(),
            //                       f->vertex(1)->point().x(),
            //                       f->vertex(1)->point().y(),
            //                       f->vertex(1)->point().z(),
            //                       f->vertex(2)->point().x(),
            //                       f->vertex(2)->point().y(),
            //                       f->vertex(2)->point().z(),
            //                       f->vertex(3)->point().x(),
            //                       f->vertex(3)->point().y(),
            //                       f->vertex(3)->point().z(), c.center().x(),
            //                       c.center().y(), c.center().z(),
            //                       std::sqrt(c.squared_radius()));

            conflictFree = false;
            break;
          }
        }
      }

      //                printf("[%llu] %lu conflicting faces with radius
      //                %llu\n",
      //                       chunk_id, conflicts.size(), cell_radius);
    }

#ifdef DEL_STATS
    radius_stats_.emplace(chunk_id,
                          del_stats(conflictFree ? cell_radius - 1 : -1,
                                    id_high - id_low, tria.number_of_vertices(),
                                    tria.number_of_finite_cells()));
#endif

    if (!conflictFree) {
      fprintf(stderr, "[%llu] EXCEPTION: triangulation did not converge\n",
              chunk_id);
      return;
    }

    // we have a conflict free triangulation, output edges
#ifdef OUTPUT_EDGES
    edge_io_.ReserveEdges(id_high - id_low);
#endif
    for (auto e = tria.finite_edges_begin(); e != tria.finite_edges_end();
         ++e) {
      // we only save outgoing edges from a vertices within our chunk

      if(tria.is_infinite(e->first))
          continue;

      auto v1 = e->first->vertex(e->second);
      auto v2 = e->first->vertex(e->third);

      if(v1 == nullptr || v2 == nullptr)
          continue;

        bool touches = false;
        if (!tria.is_infinite(v1) &&  !(v1->info() & COPY_FLAG)) {
            // v1 is in chunk we save the edge
            edge_io_.UpdateDist(v1->info());
            touches = true;
        }
        if (!tria.is_infinite(v2) &&  !(v2->info() & COPY_FLAG)) {
            // v1 is not in chunk but v2 is in chunk we save the edge
            edge_io_.UpdateDist(v2->info());
            touches = true;
        }

#ifdef OUTPUT_EDGES
        if(touches){
            edge_io_.PushEdge((v1->info() & COPY_FLAG) ? v1->info() - COPY_FLAG : v1->info(), (v2->info() & COPY_FLAG) ? v2->info() - COPY_FLAG : v2->info());
        }
#endif
    }

    /*adjacency_io_.ReserveEdges(id_high - id_low);
    for (auto n = tria.finite_vertices_begin(); n != tria.finite_vertices_end();
    ++n) {
        // we only save outgoing edges from a vertices within our chunk

        if(!(id_low <= n->info() && n->info() < id_high))
            continue;

        std::vector<Vh_3d> neighbor_handles;
        tria.adjacent_vertices(n, std::back_inserter(neighbor_handles));

        std::vector<SInt> neighbors;
        for(const auto & n : neighbor_handles){
            if(!tria.is_infinite(n)){
                neighbors.push_back(n->info());
            }
        }

        adjacency_io_.PushEdge(n->info(), neighbors);
    }*/
  }

private:

    void SortCellVertices(std::vector<Vertex> & vertices){

        struct LessX {
            bool operator()(const Vertex& p, const Vertex& q) const
            {
                return std::get<0>(p) < std::get<0>(q);
            }
        };
        struct LessY {
            bool operator()(const Vertex& p, const Vertex& q) const
            {
                return std::get<1>(p) < std::get<1>(q);
            }
        };
        struct LessZ {
            bool operator()(const Vertex& p, const Vertex& q) const
            {
                return std::get<2>(p) < std::get<2>(q);
            }
        };

        struct SpatialSortingTraits {
            typedef Vertex Point_3;
            typedef LessX Less_x_3;
            typedef LessY Less_y_3;
            typedef LessZ Less_z_3;

            Less_x_3 less_x_3_object() const
            {
                return Less_x_3();
            }
            Less_y_3 less_y_3_object() const
            {
                return Less_y_3();
            }
            Less_z_3 less_z_3_object() const
            {
                return Less_z_3();
            }
        };

        SpatialSortingTraits sst;
        CGAL::spatial_sort(vertices.begin(), vertices.end(), sst);
    }
};

#endif
