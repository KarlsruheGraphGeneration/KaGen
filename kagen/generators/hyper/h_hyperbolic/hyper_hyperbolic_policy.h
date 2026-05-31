#pragma once

#include "kagen/generators/geometric/geometric_util.h"
#include "kagen/generators/hyper/h_hyperbolic/hyper_hyperbolic.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"

#include <CGAL/number_utils_classes.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace kagen {

template <typename Double>
class HyperbolicGeometryPolicy {
public:
    using GeneratorT = Hyper_Hyperbolic<Double>;
    using Center     = HyperbolicHyperedgeCenter<Double>;
    using Vertex     = typename GeneratorT::Vertex;

    struct Cell {
        SInt annulus_id;
        SInt chunk_id;
        SInt cell_id;
    };

    explicit HyperbolicGeometryPolicy(GeneratorT& generator, const SInt annulus_id, const SInt chunk_id)
        : gen_(generator),
          start_annulus_id_(annulus_id),
          start_chunk_id_(chunk_id),
          start_cell_id_(0) {}

    void SetStartCell(const SInt annulus_id, const SInt chunk_id, const SInt cell_id) {
        start_annulus_id_ = annulus_id;
        start_chunk_id_   = chunk_id;
        start_cell_id_    = cell_id;
    }

    void AddCenter(const Center&, std::vector<SInt>&) const {}

    Double Radius(const Center&) const {
        return gen_.current_hyperedge_radius_;
    }

    std::vector<Cell> CandidateCells(const Center& center, const Double /*radius*/) {
        std::vector<Cell> cells;
        QueryCollect(start_annulus_id_, start_chunk_id_, start_cell_id_, center, true, cells);

        if (gen_.config_.query_both && start_annulus_id_ > 0) {
            auto& chunk = gen_.chunks_[start_chunk_id_];

            const Double min_chunk_phi = std::get<1>(chunk);
            const Double max_chunk_phi = std::get<2>(chunk);

            const SInt grid_size = gen_.GridSizeForAnnulus(start_annulus_id_ - 1);

            const Double grid_phi = (max_chunk_phi - min_chunk_phi) / grid_size;

            const SInt computed_cell_id = static_cast<SInt>(std::floor((center.phi - min_chunk_phi) / grid_phi));

            const SInt next_cell_id = std::clamp<SInt>(computed_cell_id, 0, grid_size - 1);
            QueryCollect(start_annulus_id_ - 1, start_chunk_id_, next_cell_id, center, false, cells);
        }

        return cells;
    }
    Double AngularDistance(Double a, Double b) const {
        a = NormalizePhi(a);
        b = NormalizePhi(b);

        Double d = std::abs(a - b);
        return std::min(d, Double{2.0 * M_PI} - d);
    }

    Double CoshDistanceWithDelta(Double center_r, Double query_r, Double delta_phi) const {
        return (std::cosh(center_r) * std::cosh(query_r))
               - (std::sinh(center_r) * std::sinh(query_r) * std::cos(delta_phi));
    }

    Double MinAngularDistanceToInterval(Double phi, Double min_phi, Double max_phi) const {
        if (AngleInInterval(phi, min_phi, max_phi)) {
            return 0.0;
        }

        return std::min(AngularDistance(phi, min_phi), AngularDistance(phi, max_phi));
    }

    bool AngleInInterval(Double phi, Double min_phi, Double max_phi) const {
        phi     = NormalizePhi(phi);
        min_phi = NormalizePhi(min_phi);
        max_phi = NormalizePhi(max_phi);

        if (min_phi <= max_phi) {
            return min_phi <= phi && phi <= max_phi;
        }

        // interval wraps around 2π, e.g. [5.9, 0.2]
        return phi >= min_phi || phi <= max_phi;
    }

    bool CellsOutside(
        const Center& center, const Double hyperball_radius, const Double min_r, const Double max_r,
        const Double min_phi, const Double max_phi) const {
        const Double min_delta = MinAngularDistanceToInterval(center.phi, min_phi, max_phi);
        Double       best_r    = center.r;
        const Double candidate = std::tanh(center.r) * std::cos(min_delta);

        if (candidate > 0.0 && candidate < 1.0) {
            best_r = std::atanh(candidate);
        }

        best_r = std::clamp(best_r, min_r, max_r);

        const Double min_cosh_dist = CoshDistanceWithDelta(center.r, best_r, min_delta);

        return min_cosh_dist > std::cosh(hyperball_radius);
    }

    CellBallRelation ClassifyCell(const Center& center, const Double radius, const Cell& cell) const {
        const auto& annulus = gen_.annuli_.find(gen_.ComputeGlobalChunkId(cell.annulus_id, cell.chunk_id))->second;
        const auto& stored_cell =
            gen_.cells_.find(gen_.ComputeGlobalCellId(cell.annulus_id, cell.chunk_id, cell.cell_id))->second;

        const Double min_r   = std::get<1>(annulus);
        const Double max_r   = std::get<2>(annulus);
        const Double min_phi = std::get<1>(stored_cell);
        const Double max_phi = std::get<2>(stored_cell);

        if (CellsOutside(center, radius, min_r, max_r, min_phi, max_phi)) {
            return CellBallRelation::OUTSIDE;
        }

        const bool inside =
            IsPointInside(center, min_r, min_phi, radius) && IsPointInside(center, min_r, max_phi, radius)
            && IsPointInside(center, max_r, min_phi, radius) && IsPointInside(center, max_r, max_phi, radius);

        return inside ? CellBallRelation::INSIDE : CellBallRelation::PARTIAL;
    }

    Double NormalizePhi(Double phi) const {
        const Double two_pi = 2.0 * M_PI;

        while (phi < 0.0) {
            phi += two_pi;
        }

        while (phi >= two_pi) {
            phi -= two_pi;
        }

        return phi;
    }

    Double IntervalOverlap(Double a_begin, Double a_end, Double b_begin, Double b_end) const {
        const Double begin = std::max(a_begin, b_begin);
        const Double end   = std::min(a_end, b_end);

        return std::max<Double>(0.0, end - begin);
    }

    std::vector<std::pair<Double, Double>> SplitCircularInterval(Double begin, Double end) const {
        const Double two_pi = 2.0 * M_PI;

        begin = NormalizePhi(begin);
        end   = NormalizePhi(end);

        if (begin <= end) {
            return {{begin, end}};
        }

        return {{begin, two_pi}, {0.0, end}};
    }

    Double AllowedHalfAngleAtRadius(const Double center_r, const Double query_r, const Double hyperball_radius) const {
        const Double eps = std::numeric_limits<Double>::epsilon();

        if (std::abs(center_r) < eps || std::abs(query_r) < eps) {
            const Double distance = std::abs(center_r - query_r);
            return distance <= hyperball_radius ? M_PI : 0.0;
        }

        const Double numerator = std::cosh(center_r) * std::cosh(query_r) - std::cosh(hyperball_radius);

        const Double denominator = std::sinh(center_r) * std::sinh(query_r);

        Double arg = numerator / denominator;

        if (arg <= -1.0) {
            return M_PI;
        }

        if (arg >= 1.0) {
            return 0.0;
        }

        return std::acos(arg);
    }

    Double CircularIntervalOverlap(Double a_begin, Double a_end, Double b_begin, Double b_end) const {
        Double overlap = 0.0;

        const auto a_parts = SplitCircularInterval(a_begin, a_end);
        const auto b_parts = SplitCircularInterval(b_begin, b_end);

        for (const auto& [ab, ae]: a_parts) {
            for (const auto& [bb, be]: b_parts) {
                overlap += IntervalOverlap(ab, ae, bb, be);
            }
        }

        return overlap;
    }

    double CellCoverage(const Center& center, const Double hyperball_radius, const Cell& cell) const {
        const CellBallRelation relation = ClassifyCell(center, hyperball_radius, cell);

        if (relation == CellBallRelation::OUTSIDE) {
            return 0.0;
        }

        if (relation == CellBallRelation::INSIDE) {
            return 1.0;
        }

        const auto& annulus = gen_.annuli_.find(gen_.ComputeGlobalChunkId(cell.annulus_id, cell.chunk_id))->second;
        const auto& stored_cell =
            gen_.cells_.find(gen_.ComputeGlobalCellId(cell.annulus_id, cell.chunk_id, cell.cell_id))->second;

        const Double min_r   = std::get<1>(annulus);
        const Double max_r   = std::get<2>(annulus);
        const Double min_phi = std::get<1>(stored_cell);
        const Double max_phi = std::get<2>(stored_cell);

        const Double center_phi = center.phi;
        const Double center_r   = center.r;

        const Double cell_phi_width = max_phi - min_phi;

        if (cell_phi_width <= 0.0 || max_r <= min_r) {
            return 0.0;
        }

        constexpr int num_samples = 32;

        Double weighted_inside = 0.0;
        Double weighted_total  = 0.0;

        const Double dr = (max_r - min_r) / static_cast<Double>(num_samples);

        for (int i = 0; i < num_samples; ++i) {
            const Double radius = min_r + ((static_cast<Double>(i) + 0.5) * dr);

            const Double half_angle = AllowedHalfAngleAtRadius(center_r, radius, hyperball_radius);

            Double overlap_phi = 0.0;

            if (half_angle >= M_PI) {
                overlap_phi = cell_phi_width;
            } else if (half_angle > 0.0) {
                overlap_phi =
                    CircularIntervalOverlap(min_phi, max_phi, center_phi - half_angle, center_phi + half_angle);
            }

            const Double alpha  = (gen_.config_.plexp - 1.0) / 2.0;
            const Double weight = std::sinh(alpha * radius);

            weighted_inside += overlap_phi * weight * dr;
            weighted_total += cell_phi_width * weight * dr;
        }

        if (weighted_total <= 0.0) {
            return 0.0;
        }

        return std::clamp(static_cast<double>(weighted_inside / weighted_total), 0.0, 1.0);
    }

    SInt AddWholeCell(const Cell& cell, std::vector<PinRange>& ranges) const {
        const SInt global_cell_id = gen_.ComputeGlobalCellId(cell.annulus_id, cell.chunk_id, cell.cell_id);

        if (gen_.cells_.find(global_cell_id) == gen_.cells_.end()) {
            return 0;
        }

        const auto& stored_cell = gen_.cells_[global_cell_id];

        const SInt size   = std::get<0>(stored_cell);
        const SInt offset = std::get<4>(stored_cell);

        if (size > 0) {
            ranges.push_back({.begin = offset, .end = offset + size});
        }
        return size;
    }

    SInt AddPartialCell(
        const Center& /*center*/, const Cell& cell, const double coverage, std::vector<PinRange>& ranges) const {
        const SInt global_cell_id = gen_.ComputeGlobalCellId(cell.annulus_id, cell.chunk_id, cell.cell_id);

        if (gen_.cells_.find(global_cell_id) == gen_.cells_.end()) {
            return 0;
        }

        const auto& stored_cell = gen_.cells_[global_cell_id];

        SInt size       = std::get<0>(stored_cell);
        SInt offset     = std::get<4>(stored_cell);
        SInt range_size = static_cast<SInt>(std::floor(static_cast<double>(size) * coverage));
        SInt seed       = gen_.config_.seed + (131 * start_cell_id_) + std::floor(coverage * 100);
        ranges.push_back(getRandomPinRange(size, range_size, offset, seed, gen_.config_));
        return range_size;
    }

    SInt AddPartialCellExact(const Center& center, const Double /*radius*/, const Cell& cell, std::vector<SInt>& pins) {
        gen_.GenerateVertices(cell.annulus_id, cell.chunk_id, cell.cell_id);

        const SInt global_cell_id = gen_.ComputeGlobalCellId(cell.annulus_id, cell.chunk_id, cell.cell_id);

        if (gen_.vertices_.find(global_cell_id) == gen_.vertices_.end()) {
            return 0;
        }

        const auto& vertices       = gen_.vertices_[global_cell_id];
        SInt        vertex_counter = 0;
        for (const Vertex& vertex: vertices) {
            if (HyperbolicDistance(center, vertex) <= gen_.current_hyperedge_pdm_radius_) {
                pins.push_back(std::get<5>(vertex));
                vertex_counter++;
            }
        }
        return vertex_counter;
    }

    void EmitHyperedge(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges) {
        gen_.PushHyperedgeCompressed(pins, ranges);
    }

    std::string CenterToString(const Center& center) const {
        std::ostringstream out;

        out << "phi=" << center.phi << ";r=" << center.r;

        return out.str();
    }

private:
    GeneratorT& gen_;
    SInt        start_annulus_id_;
    SInt        start_chunk_id_;
    SInt        start_cell_id_;

    void
    PushCandidateCell(const SInt annulus_id, const SInt chunk_id, const SInt cell_id, std::vector<Cell>& cells) const {
        cells.push_back(Cell{.annulus_id = annulus_id, .chunk_id = chunk_id, .cell_id = cell_id});
    }

    void QueryCollect(
        const SInt annulus_id, const SInt chunk_id, const SInt cell_id, const Center& center, const bool search_down,
        std::vector<Cell>& cells) {
        auto& annulus = gen_.annuli_[gen_.ComputeGlobalChunkId(annulus_id, chunk_id)];
        auto  bounds  = GetBoundaryPhis(center.phi, center.r, annulus_id);

        gen_.current_min_phi_ = bounds.first;
        gen_.current_max_phi_ = bounds.second;

        const auto& stored_cell = gen_.cells_[gen_.ComputeGlobalCellId(annulus_id, chunk_id, cell_id)];

        const Double min_cell_phi = std::get<1>(stored_cell);
        const Double max_cell_phi = std::get<2>(stored_cell);

        gen_.right_processed_chunk_ = chunk_id;
        gen_.right_processed_cell_  = cell_id;

        if (search_down || !gen_.IsLocalChunk(chunk_id)) {
            PushCandidateCell(annulus_id, chunk_id, cell_id, cells);
        }

        if (std::get<1>(annulus) >= gen_.clique_thres_
            && std::max(gen_.TotalGridSizeForAnnulus(annulus_id), gen_.config_.k) > 1) {
            if (gen_.current_min_phi_ < min_cell_phi
                || (gen_.OutOfBounds(gen_.current_min_phi_) && !(std::abs(min_cell_phi - 0.0) < gen_.cell_eps_))) {
                SInt next_chunk_id = chunk_id;

                if (cell_id == 0) {
                    next_chunk_id = (chunk_id + gen_.config_.k - 1) % gen_.config_.k;
                }

                const SInt next_cell_id =
                    (cell_id + gen_.GridSizeForAnnulus(annulus_id) - 1) % gen_.GridSizeForAnnulus(annulus_id);

                gen_.GenerateCells(annulus_id, next_chunk_id);

                QueryRightNeighborCollect(
                    annulus_id, next_chunk_id, next_cell_id, center, std::abs(min_cell_phi - 0.0) < gen_.cell_eps_,
                    search_down, cells);
            }

            if (gen_.current_max_phi_ > max_cell_phi
                || (gen_.OutOfBounds(gen_.current_max_phi_)
                    && !(std::abs(max_cell_phi - (2 * M_PI)) < gen_.cell_eps_))) {
                SInt next_chunk_id = chunk_id;

                if (cell_id == gen_.GridSizeForAnnulus(annulus_id) - 1) {
                    next_chunk_id = (chunk_id + gen_.config_.k + 1) % gen_.config_.k;
                }

                const SInt next_cell_id =
                    (cell_id + gen_.GridSizeForAnnulus(annulus_id) + 1) % gen_.GridSizeForAnnulus(annulus_id);

                gen_.GenerateCells(annulus_id, next_chunk_id);

                QueryLeftNeighborCollect(
                    annulus_id, next_chunk_id, next_cell_id, center,
                    std::abs(max_cell_phi - (2 * M_PI)) < gen_.cell_eps_, search_down, cells);
            }
        }

        const SSInt next_annulus =
            search_down ? static_cast<SSInt>(annulus_id + 1) : static_cast<SSInt>(annulus_id - 1);

        if (next_annulus < 0 || next_annulus >= static_cast<SSInt>(gen_.total_annuli_)) {
            return;
        }

        auto& chunk = gen_.chunks_[chunk_id];

        const Double min_chunk_phi = std::get<1>(chunk);
        const Double max_chunk_phi = std::get<2>(chunk);
        const Double grid_phi      = (max_chunk_phi - min_chunk_phi) / gen_.GridSizeForAnnulus(next_annulus);

        const SInt next_cell_id = static_cast<SInt>(std::floor((center.phi - min_chunk_phi) / grid_phi));

        QueryCollect(next_annulus, chunk_id, next_cell_id, center, search_down, cells);
    }

    void QueryRightNeighborCollect(
        const SInt annulus_id, SInt chunk_id, SInt cell_id, const Center& /*center*/, bool phase,
        const bool search_down, std::vector<Cell>& cells) {
        while (true) {
            if (phase && gen_.current_min_phi_ < 0.0) {
                gen_.current_min_phi_ += 2 * M_PI;
            }

            if (phase && gen_.OutOfBounds(gen_.current_min_phi_)) {
                return;
            }

            auto& stored_cell = gen_.cells_[gen_.ComputeGlobalCellId(annulus_id, chunk_id, cell_id)];

            const Double min_cell_phi = std::get<1>(stored_cell);

            gen_.right_processed_chunk_ = chunk_id;
            gen_.right_processed_cell_  = cell_id;

            if (search_down || !gen_.IsLocalChunk(chunk_id)) {
                PushCandidateCell(annulus_id, chunk_id, cell_id, cells);
            }

            phase = phase || std::abs(min_cell_phi - 0.0) < gen_.cell_eps_;

            if (gen_.current_min_phi_ < min_cell_phi || gen_.OutOfBounds(gen_.current_min_phi_)) {
                SInt next_chunk_id = chunk_id;

                if (cell_id == 0) {
                    next_chunk_id = (chunk_id + gen_.config_.k - 1) % gen_.config_.k;
                }

                const SInt next_cell_id =
                    (cell_id + gen_.GridSizeForAnnulus(annulus_id) - 1) % gen_.GridSizeForAnnulus(annulus_id);

                gen_.GenerateCells(annulus_id, next_chunk_id);

                cell_id  = next_cell_id;
                chunk_id = next_chunk_id;
                continue;
            }

            return;
        }
    }

    void QueryLeftNeighborCollect(
        const SInt annulus_id, SInt chunk_id, SInt cell_id, const Center& /*center*/, bool phase,
        const bool search_down, std::vector<Cell>& cells) {
        while (true) {
            if (phase && gen_.current_max_phi_ >= 2 * M_PI) {
                gen_.current_max_phi_ -= 2 * M_PI;
            }

            if (phase && gen_.OutOfBounds(gen_.current_max_phi_)) {
                return;
            }

            if (chunk_id == gen_.right_processed_chunk_ && cell_id == gen_.right_processed_cell_) {
                return;
            }

            auto& stored_cell = gen_.cells_[gen_.ComputeGlobalCellId(annulus_id, chunk_id, cell_id)];

            const Double max_cell_phi = std::get<2>(stored_cell);

            if (search_down || !gen_.IsLocalChunk(chunk_id)) {
                PushCandidateCell(annulus_id, chunk_id, cell_id, cells);
            }

            phase = phase || std::abs(max_cell_phi - (2 * M_PI)) < gen_.cell_eps_;

            if (gen_.current_max_phi_ > max_cell_phi || gen_.OutOfBounds(gen_.current_max_phi_)) {
                SInt next_chunk_id = chunk_id;

                if (cell_id == gen_.GridSizeForAnnulus(annulus_id) - 1) {
                    next_chunk_id = (chunk_id + gen_.config_.k + 1) % gen_.config_.k;
                }

                const SInt next_cell_id =
                    (cell_id + gen_.GridSizeForAnnulus(annulus_id) + 1) % gen_.GridSizeForAnnulus(annulus_id);

                gen_.GenerateCells(annulus_id, next_chunk_id);

                cell_id  = next_cell_id;
                chunk_id = next_chunk_id;
                continue;
            }

            return;
        }
    }

    bool IsPointInside(const Center& center, Double r, Double phi, Double radius) const {
        Double delta_phi = std::abs(center.phi - phi);
        delta_phi        = std::min(delta_phi, Double{2.0 * M_PI} - delta_phi);

        const Double cosh_dist =
            std::cosh(center.r) * std::cosh(r) - std::sinh(center.r) * std::sinh(r) * std::cos(delta_phi);

        return cosh_dist <= std::cosh(radius);
    }

    Double HyperbolicDistance(const Center& center, const Vertex& vertex) const {
        Double delta_phi = std::abs(center.phi - std::get<0>(vertex));
        delta_phi        = std::min(delta_phi, Double{2.0 * M_PI} - delta_phi);

        const Double cosh_dist = std::cosh(center.r) * std::cosh(std::get<1>(vertex))
                                 - std::sinh(center.r) * std::sinh(std::get<1>(vertex)) * std::cos(delta_phi);

        return (cosh_dist - 1.0) / 2.0;
    }

    std::pair<Double, Double>
    GetBoundaryPhis(const Double boundary_phi, const Double boundary_r, const SInt annulus_id) const {
        const auto& boundary = gen_.boundaries_[annulus_id];

        const Double cosh_radius = std::cosh(gen_.current_hyperedge_radius_);
        const Double cosh_min_r  = std::get<0>(boundary);
        const Double sinh_min_r  = std::get<1>(boundary);
        if (std::abs(std::sinh(boundary_r) * sinh_min_r) <= std::numeric_limits<Double>::epsilon()) {
            return {boundary_phi - M_PI, boundary_phi + M_PI};
        }

        Double arg = ((std::cosh(boundary_r) * cosh_min_r) - cosh_radius) / (std::sinh(boundary_r) * sinh_min_r);
        arg        = std::clamp(arg, Double{-1.0}, Double{1.0});

        const Double diff = std::acos(arg);

        return {boundary_phi - diff, boundary_phi + diff};
    }
};

} // namespace kagen
