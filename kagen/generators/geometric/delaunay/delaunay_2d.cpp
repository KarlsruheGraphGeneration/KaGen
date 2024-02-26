#include "kagen/generators/geometric/delaunay/delaunay_2d.h"

#include "kagen/tools/geometry.h"

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/spatial_sort.h>

namespace kagen {
using K_2d  = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb_2d = CGAL::Triangulation_vertex_base_with_info_2<kagen::SInt, K_2d>; // attach an ID to
                                                                              // each point
using Fb_2d    = CGAL::Triangulation_face_base_2<K_2d>;
using Tds_2d   = CGAL::Triangulation_data_structure_2<Vb_2d, Fb_2d>;
using Dt_2d    = CGAL::Delaunay_triangulation_2<K_2d, Tds_2d>;
using Point_2d = Dt_2d::Point;
using Fh_2d    = Dt_2d::Face_handle;

using Circ_2d = CGAL::Circle_2<K_2d>;
using Box_2d  = CGAL::Bbox_2;

Delaunay2D::Delaunay2D(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : Geometric2D(config, rank, size) {
    // Chunk variables
    total_chunks_   = config_.k;
    chunks_per_dim_ = std::sqrt(total_chunks_);
    chunk_size_     = 1.0 / chunks_per_dim_;

    // Cell variables

    // upper bound for expected distance to nearest neighbor, according to paper
    LPFloat edNN      = 4.8 / std::sqrt(config.n);
    LPFloat lCellSize = std::min(edNN, chunk_size_);

    // use chunks_per_dim_ as upper bound for cells_per_dim_
    cells_per_dim_   = std::min(SInt(128), (SInt)floor(chunk_size_ / lCellSize));
    cells_per_chunk_ = cells_per_dim_ * cells_per_dim_;
    cell_size_       = chunk_size_ / cells_per_dim_;

    max_radius_ = (cells_per_dim_ * chunks_per_dim_);

    InitDatastructures();
}

void Delaunay2D::GenerateEdges(const SInt chunk_row, const SInt chunk_column) {
    SInt chunk_id = Encode(chunk_column, chunk_row);

    Dt_2d tria;
    Fh_2d hint;

    // bounding box of own chunk
    Box_2d bbChunk(
        chunk_row * chunk_size_, chunk_column * chunk_size_, (chunk_row + 1) * chunk_size_,
        (chunk_column + 1) * chunk_size_);

    // Iterate chunk cells
    for (SInt cell_row = 0; cell_row < cells_per_dim_; ++cell_row) {
        for (SInt cell_column = 0; cell_column < cells_per_dim_; ++cell_column) {
            SInt cell_id = cell_row * cells_per_dim_ + cell_column;
            SInt cellIdx = ComputeGlobalCellId(chunk_id, cell_id);

            SortCellVertices(vertices_[cellIdx]);
            for (const auto& v: vertices_[cellIdx]) {
                Point_2d p(std::get<0>(v), std::get<1>(v));
                assert(bbChunk.xmin() <= p.x() && p.x() <= bbChunk.xmax());
                assert(bbChunk.ymin() <= p.y() && p.y() <= bbChunk.ymax());

                auto vh    = tria.insert(p, hint);
                vh->info() = std::get<2>(v);
                hint       = vh->face();
            }
        }
    }

    bool conflictFree = false;
    SInt cell_radius  = 1;
    for (; !conflictFree && cell_radius <= max_radius_; ++cell_radius) {
        // bounding box of neighborhood
        Box_2d bbNH(
            chunk_row * chunk_size_ - cell_radius * cell_size_, chunk_column * chunk_size_ - cell_radius * cell_size_,
            (chunk_row + 1) * chunk_size_ + cell_radius * cell_size_,
            (chunk_column + 1) * chunk_size_ + cell_radius * cell_size_);

        SInt chunk_radius = (SInt)std::ceil(cell_radius / static_cast<LPFloat>((cells_per_dim_)));
        for (SSInt chunk_row_diff = -chunk_radius; chunk_row_diff <= (SSInt)chunk_radius; ++chunk_row_diff) {
            for (SSInt chunk_col_diff = -chunk_radius; chunk_col_diff <= (SSInt)chunk_radius; ++chunk_col_diff) {
                if ((SInt)std::max(std::abs(chunk_row_diff), std::abs(chunk_col_diff)) < chunk_radius)
                    continue; // skip chunks in inner region (already added)

                SSInt neighbor_chunk_col = chunk_column + chunk_col_diff;
                SSInt neighbor_chunk_row = chunk_row + chunk_row_diff;

                if (!config_.periodic) {
                    if (neighbor_chunk_row < 0 || static_cast<SInt>(neighbor_chunk_row) >= chunks_per_dim_) {
                        continue;
                    }
                    if (neighbor_chunk_col < 0 || static_cast<SInt>(neighbor_chunk_col) >= chunks_per_dim_) {
                        continue;
                    }
                }

                LPFloat x_offset = 0;
                if (neighbor_chunk_row < 0) {
                    x_offset = -1.0 * (std::ceil(std::abs(neighbor_chunk_row) / static_cast<LPFloat>(chunks_per_dim_)));
                } else {
                    x_offset = 1.0 * (std::floor(std::abs(neighbor_chunk_row) / static_cast<LPFloat>(chunks_per_dim_)));
                }

                LPFloat y_offset = 0;
                if (neighbor_chunk_col < 0) {
                    y_offset = -1.0 * (std::ceil(std::abs(neighbor_chunk_col) / static_cast<LPFloat>(chunks_per_dim_)));
                } else {
                    y_offset = 1.0 * (std::floor(std::abs(neighbor_chunk_col) / static_cast<LPFloat>(chunks_per_dim_)));
                }

                // Get correct grid cells
                neighbor_chunk_row = (neighbor_chunk_row % (SSInt)chunks_per_dim_ + chunks_per_dim_) % chunks_per_dim_;
                neighbor_chunk_col = (neighbor_chunk_col % (SSInt)chunks_per_dim_ + chunks_per_dim_) % chunks_per_dim_;
                SInt neighbor_chunk_id = Encode((SInt)neighbor_chunk_col, (SInt)neighbor_chunk_row);

                // now we have our neighbor chunk -> dive into its cells
                SInt row_cell_diff =
                    std::min((cell_radius - (std::abs(chunk_row_diff) - 1) * cells_per_dim_), cells_per_dim_);
                SInt col_cell_diff =
                    std::min((cell_radius - (std::abs(chunk_col_diff) - 1) * cells_per_dim_), cells_per_dim_);

                for (SInt dcr = 0; dcr < row_cell_diff; ++dcr) {
                    for (SInt dcc = 0; dcc < col_cell_diff; ++dcc) {
                        SInt neighbor_cell_row = (chunk_row_diff < 0 ? cells_per_dim_ - row_cell_diff : 0) + dcr;
                        SInt neighbor_cell_col = (chunk_col_diff < 0 ? cells_per_dim_ - col_cell_diff : 0) + dcc;

                        SInt total_row_diff =
                            chunk_row_diff != 0
                                ? (chunk_row_diff < 0 ? cells_per_dim_ - neighbor_cell_row : neighbor_cell_row + 1)
                                      + (std::abs(chunk_row_diff) - 1) * cells_per_dim_
                                : 0;
                        SInt total_col_diff =
                            chunk_col_diff != 0
                                ? (chunk_col_diff < 0 ? cells_per_dim_ - neighbor_cell_col : neighbor_cell_col + 1)
                                      + (std::abs(chunk_col_diff) - 1) * cells_per_dim_
                                : 0;

                        if (std::max(total_row_diff, total_col_diff) < cell_radius)
                            continue; // skip cell in inner region (already added)

                        SInt neighbor_cell_id = neighbor_cell_row * cells_per_dim_ + neighbor_cell_col;

                        // Check if vertices not generated
                        SInt neighbor_cell_offset = ComputeGlobalCellId(neighbor_chunk_id, neighbor_cell_id);
                        // lazily generate vertices
                        GenerateVertices(neighbor_chunk_id, neighbor_cell_id, false);

                        SortCellVertices(vertices_[neighbor_cell_offset]);
                        for (const auto& v: vertices_[neighbor_cell_offset]) {
                            Point_2d p(std::get<0>(v) + x_offset, std::get<1>(v) + y_offset);
                            assert(bbNH.xmin() <= p.x() && p.x() <= bbNH.xmax());
                            assert(bbNH.ymin() <= p.y() && p.y() <= bbNH.ymax());

                            auto vh    = tria.insert(p, hint);
                            vh->info() = std::get<2>(v) + COPY_FLAG;
                            hint       = vh->face();
                        }
                    }
                }
            }
        }

        conflictFree = true;
        for (auto f = tria.finite_faces_begin(); f != tria.finite_faces_end(); ++f) {
            bool touchesChunk = (!(f->vertex(0)->info() & COPY_FLAG)) || (!(f->vertex(1)->info() & COPY_FLAG))
                                || (!(f->vertex(2)->info() & COPY_FLAG));

            if (touchesChunk) {
                Circ_2d c(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point());

                if (!boxContains(bbNH, c)) {
                    conflictFree = false;
                    break;
                }
            }
        }
    }

    if (config_.periodic && !conflictFree) {
        fprintf(stderr, "[%llu] EXCEPTION: triangulation did not converge\n", chunk_id);
        return;
    }

    auto get_vertex_id = [](auto& handle) {
        return (handle->info() & COPY_FLAG) ? handle->info() - COPY_FLAG : handle->info();
    };

    for (auto vh: tria.finite_vertex_handles()) {
        // Only consider vertex if it is within this PE's bounding box
        auto p = vh->point();
        bool owned =
            bbChunk.xmin() <= p.x() && p.x() <= bbChunk.xmax() && bbChunk.ymin() <= p.y() && p.y() <= bbChunk.ymax();
        if (!owned) {
            continue;
        }

        // Iterate over incident vertices
        auto vc   = tria.incident_vertices(vh);
        auto done = decltype(vc)(vc);
        do {
            if (tria.is_infinite(vc)) {
                continue;
            }
            PushEdge(get_vertex_id(vh), get_vertex_id(vc));
        } while (++vc != done);
    }

    // we have a conflict free triangulation, output edges
    //    for (auto e = tria.finite_edges_begin(); e != tria.finite_edges_end(); ++e) {
    //        // we only save outgoing edges from a vertices within our chunk
    //
    //        if (tria.is_infinite(e->first))
    //            continue;
    //
    //        auto v1 = e->first->vertex(e->second);
    //        auto v2 = e->first->vertex(Dt_2d::cw(e->second));
    //
    //        if (v1 == nullptr || v2 == nullptr)
    //            continue;
    //
    //        // bool touches = true;
    //        if (!tria.is_infinite(v1) && !(v1->info() & COPY_FLAG)) {
    //            // v1 is in chunk we save the edge
    //            PushEdge(
    //                (v1->info() & COPY_FLAG) ? v1->info() - COPY_FLAG : v1->info(),
    //                (v2->info() & COPY_FLAG) ? v2->info() - COPY_FLAG : v2->info());
    //        }
    //        if (!tria.is_infinite(v2) && !(v2->info() & COPY_FLAG)) {
    //            // v1 is not in chunk but v2 is in chunk we save the reverse edge
    //            PushEdge(
    //                (v2->info() & COPY_FLAG) ? v2->info() - COPY_FLAG : v2->info(),
    //                (v1->info() & COPY_FLAG) ? v1->info() - COPY_FLAG : v1->info());
    //        }
    //    }
}

void Delaunay2D::SortCellVertices(std::vector<Vertex>& vertices) {
    struct LessX {
        bool operator()(const Vertex& p, const Vertex& q) const {
            return std::get<0>(p) < std::get<0>(q);
        }
    };
    struct LessY {
        bool operator()(const Vertex& p, const Vertex& q) const {
            return std::get<1>(p) < std::get<1>(q);
        }
    };

    struct SpatialSortingTraits {
        typedef Vertex Point_2;
        typedef LessX  Less_x_2;
        typedef LessY  Less_y_2;

        Less_x_2 less_x_2_object() const {
            return Less_x_2();
        }
        Less_y_2 less_y_2_object() const {
            return Less_y_2();
        }
    };

    SpatialSortingTraits sst;
    CGAL::spatial_sort(vertices.begin(), vertices.end(), sst);
}
} // namespace kagen
