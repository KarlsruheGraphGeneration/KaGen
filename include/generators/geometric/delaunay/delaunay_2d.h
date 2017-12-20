/*******************************************************************************
 * include/generators/geometric/delaunay/delaunay_2d.h
 *
 * Copyright (C) 2017 Daniel Funke <funke@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _DELAUNAY_2D_H_
#define _DELAUNAY_2D_H_

#include "geometric/geometric_2d.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/spatial_sort.h>

#include <sys/stat.h>

using K_2d = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb_2d =
    CGAL::Triangulation_vertex_base_with_info_2<SInt, K_2d>; // attach an ID to
                                                             // each point
using Fb_2d = CGAL::Triangulation_face_base_2<K_2d>;
using Tds_2d = CGAL::Triangulation_data_structure_2<Vb_2d, Fb_2d>;
using Dt_2d = CGAL::Delaunay_triangulation_2<K_2d, Tds_2d>;
using Point_2d = Dt_2d::Point;
using Fh_2d = Dt_2d::Face_handle;
using Vh_2d = Dt_2d::Vertex_handle;

using Points_2d = std::vector<std::pair<Point_2d, SInt>>; // needed to add
                                                          // points with ID to
                                                          // CGAL

using Circ_2d = CGAL::Circle_2<K_2d>;
using Box_2d = CGAL::Bbox_2;

#ifdef DEBUG_CAIRO
#include "tools/cairo_drawer.h"

class CGALCairoPainter : public CairoPainter {
public:
  CGALCairoPainter(const Point_2d &a, const Point_2d &b)
      : CairoPainter(a[0], a[1], b[0], b[1]) {}
  CGALCairoPainter(float x1 = 0, float y1 = 0, float x2 = 1, float y2 = 1)
      : CairoPainter(x1, y1, x2, y2) {}

  CGALCairoPainter(const CGALCairoPainter &a) : CairoPainter(a) {}

  void drawPoint(const Point_2d &p, uint id = ~((uint)0)) {
    CairoPainter::drawPoint(p[0], p[1], id);
  }

  void drawPoints(const Points_2d &ps) {
    for (const auto &p : ps) {
      drawPoint(p.first, p.second);
    }
  }

  void drawLine(const Point_2d &a, const Point_2d &b, bool dashed = false) {
    CairoPainter::drawLine(a[0], a[1], b[0], b[1], dashed);
  }

  void drawTriangle(const Fh_2d &f, bool dashed = false) {
    drawLine(f->vertex(0)->point(), f->vertex(1)->point(), dashed);
    drawLine(f->vertex(1)->point(), f->vertex(2)->point(), dashed);
    drawLine(f->vertex(2)->point(), f->vertex(0)->point(), dashed);
  }

  void drawTriangles(const Dt_2d &dt, bool dashed = false) {
    for (auto it = dt.finite_faces_begin(); it != dt.finite_faces_end(); ++it) {
      drawTriangle(it, dashed);
    }
  }

  void drawRectangle(const Box_2d &b, bool dashed = false) {
    CairoPainter::drawRectangle(b.xmin(), b.ymin(), b.xmax(), b.ymax(), dashed);
  }

  void drawCircle(const Circ_2d &c, bool dashed = false) {
    CairoPainter::drawCircle(c.center()[0], c.center()[1],
                             std::sqrt(c.squared_radius()), dashed);
  }
};
#endif

class Delaunay2D : public Geometric2D {
public:
  Delaunay2D(const PGeneratorConfig &config, const PEID rank)
      : Geometric2D(config, rank), point_io_(config), edge_io_(config) {
    // Chunk variables
    total_chunks_ = config_.k;
    chunks_per_dim_ = sqrt(total_chunks_);
    chunk_size_ = 1.0 / chunks_per_dim_;

    // Cell variables

    // upper bound for expected distance to nearest neighbor, according to paper
    LPFloat edNN = 4.8 / sqrt(config.n);
    LPFloat lCellSize = std::min(edNN, chunk_size_);

    // use chunks_per_dim_ as upper bound for cells_per_dim_
    cells_per_dim_ = std::min(SInt(128), (SInt)floor(chunk_size_ / lCellSize));
    cells_per_chunk_ = cells_per_dim_ * cells_per_dim_;
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
    // adjacency_io_.OutputEdges();
    // rng_.output();
  }

  inline SInt NumberOfEdges() const override { return edge_io_.NumEdges(); }

private:
  // I/O
  GeneratorIO<std::tuple<LPFloat, LPFloat>> point_io_;
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

    //fprintf(fout, "STATS begin\n");
//    fprintf(fout, "%llu %llu %llu %f %llu %llu %f %llu\n", config_.n,
//            total_chunks_, chunks_per_dim_, chunk_size_, cells_per_chunk_,
//            cells_per_dim_, cell_size_, max_radius_);

    for (const auto &s : radius_stats_) {
//      SInt x, y;
//      Decode(s.first, x, y);
      fprintf(fout, "%u %llu %llu %llu %f %llu %llu %f %llu %llu %i %llu %llu %llu\n", 2, config_.n,
              total_chunks_, chunks_per_dim_, chunk_size_, cells_per_chunk_,
              cells_per_dim_, cell_size_, max_radius_, s.first,
              s.second.radius, s.second.n_p_cell, s.second.n_p_tria,
              s.second.n_simplices);
    }

    //fprintf(fout, "STATS end\n");

    fclose(fout);
  };
#endif

  SInt max_radius_;
  static constexpr SInt COPY_FLAG = SInt(1) << (sizeof(SInt) * CHAR_BIT - 1);

  void GenerateEdges(const SInt chunk_row, const SInt chunk_column) override {
    SInt chunk_id = Encode(chunk_column, chunk_row);
    SInt id_low = std::get<4>(chunks_[chunk_id]);
    SInt id_high = id_low + std::get<0>(chunks_[chunk_id]);
    //            CGALCairoPainter basePainter;

    Dt_2d tria;
    Fh_2d hint;

    // bounding box of own chunk
    Box_2d bbChunk(chunk_row * chunk_size_, chunk_column * chunk_size_,
                   (chunk_row + 1) * chunk_size_,
                   (chunk_column + 1) * chunk_size_);

    //            basePainter.setColor(0,0,0);
    //            basePainter.drawRectangle(bbChunk);

    //            printf("[%llu] chunk (%llu, %llu) radius %llu\n",
    //                   chunk_id, chunk_row, chunk_column, cell_radius);

    // Iterate chunk cells
    for (SInt cell_row = 0; cell_row < cells_per_dim_; ++cell_row) {
      for (SInt cell_column = 0; cell_column < cells_per_dim_; ++cell_column) {
        SInt cell_id = cell_row * cells_per_dim_ + cell_column;
        SInt cellIdx = ComputeGlobalCellId(chunk_id, cell_id);

        //                    printf("[%llu] adding %lu points from own cell
        //                    %llu (%llu, %llu)\n",
        //                           chunk_id, vertices_[cellIdx].size(),
        //                           cell_id, cell_row, cell_column);

        SortCellVertices(vertices_[cellIdx]);
        for (const auto &v : vertices_[cellIdx]) {
          Point_2d p(std::get<0>(v), std::get<1>(v));
          assert(bbChunk.xmin() <= p.x() && p.x() <= bbChunk.xmax());
          assert(bbChunk.ymin() <= p.y() && p.y() <= bbChunk.ymax());

          auto vh = tria.insert(p, hint);
          vh->info() = std::get<2>(v);
          hint = vh->face();
        }
      }
    }

    //            printf("[%llu] %lu own points\n",
    //                   chunk_id, points.size());
    //            basePainter.drawPoints(points);

    bool conflictFree = false;
    SInt cell_radius = 1;
    for (; !conflictFree && cell_radius <= max_radius_; ++cell_radius) {
      // bounding box of neighborhood
      Box_2d bbNH(chunk_row * chunk_size_ - cell_radius * cell_size_,
                  chunk_column * chunk_size_ - cell_radius * cell_size_,
                  (chunk_row + 1) * chunk_size_ + cell_radius * cell_size_,
                  (chunk_column + 1) * chunk_size_ + cell_radius * cell_size_);

      //                CGALCairoPainter painter = basePainter;
      //                painter.setColor(1, 0, 0);
      //                painter.drawRectangle(bbNH, true);

      SInt chunk_radius =
          (SInt)std::ceil(cell_radius / static_cast<LPFloat>((cells_per_dim_)));
      for (SSInt chunk_row_diff = -chunk_radius;
           chunk_row_diff <= (SSInt)chunk_radius; ++chunk_row_diff) {
        for (SSInt chunk_col_diff = -chunk_radius;
             chunk_col_diff <= (SSInt)chunk_radius; ++chunk_col_diff) {
          if ((SInt)std::max(std::abs(chunk_row_diff),
                             std::abs(chunk_col_diff)) < chunk_radius)
            continue; // skip chunks in inner region (already added)

          SSInt neighbor_chunk_col = chunk_column + chunk_col_diff;
          SSInt neighbor_chunk_row = chunk_row + chunk_row_diff;

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

          // Get correct grid cells
          neighbor_chunk_row =
              (neighbor_chunk_row % (SSInt)chunks_per_dim_ + chunks_per_dim_) %
              chunks_per_dim_;
          neighbor_chunk_col =
              (neighbor_chunk_col % (SSInt)chunks_per_dim_ + chunks_per_dim_) %
              chunks_per_dim_;
          SInt neighbor_chunk_id =
              Encode((SInt)neighbor_chunk_col, (SInt)neighbor_chunk_row);

          // now we have our neighbor chunk -> dive into its cells
          SInt row_cell_diff = std::min(
              (cell_radius - (std::abs(chunk_row_diff) - 1) * cells_per_dim_),
              cells_per_dim_);
          SInt col_cell_diff = std::min(
              (cell_radius - (std::abs(chunk_col_diff) - 1) * cells_per_dim_),
              cells_per_dim_);

          for (SInt dcr = 0; dcr < row_cell_diff; ++dcr) {
            for (SInt dcc = 0; dcc < col_cell_diff; ++dcc) {
              SInt neighbor_cell_row =
                  (chunk_row_diff < 0 ? cells_per_dim_ - row_cell_diff : 0) +
                  dcr;
              SInt neighbor_cell_col =
                  (chunk_col_diff < 0 ? cells_per_dim_ - col_cell_diff : 0) +
                  dcc;

              SInt total_row_diff =
                  chunk_row_diff != 0
                      ? (chunk_row_diff < 0 ? cells_per_dim_ - neighbor_cell_row
                                            : neighbor_cell_row + 1) +
                            (std::abs(chunk_row_diff) - 1) * cells_per_dim_
                      : 0;
              SInt total_col_diff =
                  chunk_col_diff != 0
                      ? (chunk_col_diff < 0 ? cells_per_dim_ - neighbor_cell_col
                                            : neighbor_cell_col + 1) +
                            (std::abs(chunk_col_diff) - 1) * cells_per_dim_
                      : 0;

              if (std::max(total_row_diff, total_col_diff) < cell_radius)
                continue; // skip cell in inner region (already added)

              SInt neighbor_cell_id =
                  neighbor_cell_row * cells_per_dim_ + neighbor_cell_col;

              // Check if vertices not generated
              SInt neighbor_cell_offset =
                  ComputeGlobalCellId(neighbor_chunk_id, neighbor_cell_id);
              // lazily generate vertices
              GenerateVertices(neighbor_chunk_id, neighbor_cell_id);

              // Gather vertices
              //                            printf("[%llu] adding %lu points
              //                            from chunk %llu (%llu,%llu) cell
              //                            %llu (%llu, %llu) with offset
              //                            (%f,%f)\n",
              //                                   chunk_id,
              //                                   vertices_[neighbor_cell_offset].size(),
              //                                   neighbor_chunk_id,
              //                                   neighbor_chunk_row,
              //                                   neighbor_chunk_col,
              //                                   neighbor_cell_id,
              //                                   neighbor_cell_row,
              //                                   neighbor_cell_col, x_offset,
              //                                   y_offset);

              SortCellVertices(vertices_[neighbor_cell_offset]);
              for (const auto &v : vertices_[neighbor_cell_offset]) {
                Point_2d p(std::get<0>(v) + x_offset,
                           std::get<1>(v) + y_offset);
                assert(bbNH.xmin() <= p.x() && p.x() <= bbNH.xmax());
                assert(bbNH.ymin() <= p.y() && p.y() <= bbNH.ymax());

                auto vh = tria.insert(p, hint);
                vh->info() = std::get<2>(v) + COPY_FLAG;
                hint = vh->face();
              }
            }
          }
        }
      }
      //            printf("[%llu] level %llu %lu total points\n",
      //                   chunk_id, cell_radius, points.size());

      //                painter.setColor(1, 0, 0, .5);
      //                painter.drawPoints(points);

      //            printf("[%llu] Triangulation of %lu points with %lu
      //            faces\n",
      //                   chunk_id, tria.number_of_vertices(),
      //                   tria.number_of_faces());

      // painter.drawTriangles(tria);

      conflictFree = true;
      for (auto f = tria.finite_faces_begin(); f != tria.finite_faces_end(); ++f) {
        bool touchesChunk =    (!(f->vertex(0)->info() & COPY_FLAG))
                            || (!(f->vertex(1)->info() & COPY_FLAG))
                            || (!(f->vertex(2)->info() & COPY_FLAG));

        if (touchesChunk) {
          Circ_2d c(f->vertex(0)->point(),
                    f->vertex(1)->point(),
                    f->vertex(2)->point());

          if (!boxContains(bbNH, c)) {
            //                            printf("simplex (%f, %f) (%f, %f) (%f,
            //                            %f), circumsphere (%f, %f) - %f\n",
            //                                   f->vertex(0)->point().x(),
            //                                   f->vertex(0)->point().y(),
            //                                   f->vertex(1)->point().x(),
            //                                   f->vertex(1)->point().y(),
            //                                   f->vertex(2)->point().x(),
            //                                   f->vertex(2)->point().y(),
            //                                   c.center().x(), c.center().y(),
            //                                   std::sqrt(c.squared_radius()));

            conflictFree = false;
            break;
          }
        }
      }

      //                printf("[%llu] %lu conflicting faces with radius
      //                %llu\n",
      //                       chunk_id, conflicts.size(), cell_radius);

      //                painter.save(std::to_string(chunk_id) + "_level_" +
      //                std::to_string(cell_radius));
    }

#ifdef DEL_STATS
    radius_stats_.emplace(chunk_id,
                          del_stats(conflictFree ? cell_radius - 1 : -1,
                                    id_high - id_low, tria.number_of_vertices(),
                                    tria.number_of_faces()));
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
      auto v2 = e->first->vertex(Dt_2d::cw(e->second));

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
          edge_io_.PushEdge(v1->info(), v2->info());
      }
#endif
    }

    /*adjacency_io_.ReserveEdges(id_high - id_low);
    for (auto n = tria.finite_vertices_begin(); n != tria.finite_vertices_end();
    ++n) {
        // we only save outgoing edges from a vertices within our chunk

        if(!(id_low <= n->info() && n->info() < id_high))
            continue;

        std::vector<SInt> neighbors;
        auto circulator = n->incident_vertices(), done(circulator);
        do {
            if(!tria.is_infinite(circulator)){
                neighbors.push_back(circulator->info());
            }
        } while(++circulator != done);

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

        struct SpatialSortingTraits {
            typedef Vertex Point_2;
            typedef LessX Less_x_2;
            typedef LessY Less_y_2;

            Less_x_2 less_x_2_object() const
            {
                return Less_x_2();
            }
            Less_y_2 less_y_2_object() const
            {
                return Less_y_2();
            }
        };

        SpatialSortingTraits sst;
        CGAL::spatial_sort(vertices.begin(), vertices.end(), sst);
    }
};

#endif
