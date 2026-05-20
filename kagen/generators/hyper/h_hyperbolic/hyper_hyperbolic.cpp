#include "kagen/generators/hyper/h_hyperbolic/hyper_hyperbolic.h"

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"
#include "kagen/tools/geometry.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace kagen {

constexpr double M_N_D_RATIO                     = 2.0;
constexpr double HIGH_PRECISION_THRESHOLD_LOG2_N = 29.0;

// Use epsilon as 0.1% of angular resolution to stabilize boundary comparisons.
template <typename Double>
constexpr Double EPSILON_SCALE = 1000.0;

constexpr double ZERO_TOLERANCE = 1e-8;

PGeneratorConfig
Hyper_HyperbolicFactory::NormalizeParameters(PGeneratorConfig config, PEID rank, PEID size, const bool output) const {
    config.setChunkSizeIfMissing(size);

    EnsureOneChunkPerPE(config, size);

    if (config.avg_degree == 0) {
        if (config.m == 0 || config.n == 0) {
            throw ConfigurationError("at least two parameters out of {n, m, d} must be nonzero");
        }

        config.avg_degree = M_N_D_RATIO * static_cast<double>(config.m) / static_cast<double>(config.n);

        if (output) {
            std::cout << "Setting average degree to " << config.avg_degree << '\n';
        }
    } else if (config.n == 0) {
        if (config.avg_degree == 0 || config.m == 0) {
            throw ConfigurationError("at least two parameters out of {n, m, d} must be nonzero");
        }

        config.n = static_cast<SInt>(M_N_D_RATIO * static_cast<double>(config.m) / config.avg_degree);

        if (output) {
            std::cout << "Setting number of nodes to " << config.n << '\n';
        }
    }

    config.is_hypergraph = true;

    RandomRadiusChecks(config);

    if (config.streaming) {
        if (config.k < 1) {
            throw ConfigurationError("Number of chunks must be at least 1");
        }

        if (config.k > config.n) {
            throw ConfigurationError("Number of chunks must not exceed number of nodes");
        }
    }

    const HPFloat alpha = (config.plexp - 1) / 2;

    if (!PGGeometry<HPFloat>::TestTargetRadius(
            config.n, static_cast<double>(config.n) + (config.avg_degree / 2), alpha)) {
        using namespace std::string_literals;

        throw ConfigurationError(
            "generator configuration with n="s + std::to_string(config.n) + ", avg_degree="
            + std::to_string(config.avg_degree) + " and gamma=" + std::to_string(config.plexp) + " is infeasible");
    }

    if (config.hp_floats == 0 && std::log2(config.n) > HIGH_PRECISION_THRESHOLD_LOG2_N) {
        if (output) {
            std::cout << "Enabling high-precision FP for hyperbolic hypergraph generator" << '\n';
        }

        config.hp_floats = 1;
    }

    return config;
}

std::unique_ptr<Generator>
Hyper_HyperbolicFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    if (config.hp_floats != 0) {
        return std::make_unique<HighPrecisionHyperHyperbolic>(config, rank, size);
    }

    return std::make_unique<LowPrecisionHyperHyperbolic>(config, rank, size);
}

template <typename Double>
Hyper_Hyperbolic<Double>::Hyper_Hyperbolic(const PGeneratorConfig& config, PEID rank, PEID size)
    : CSROnlyGenerator(),
      config_(config),
      rank_(rank),
      size_(size),
      rng_(config) {
    alpha_    = (config_.plexp - 1) / 2;
    target_r_ = PGGeometry<Double>::GetTargetRadius(config_.n, config_.n * config_.avg_degree / 2, alpha_);

    cosh_target_r_ = std::cosh(target_r_);
    pdm_target_r_  = (cosh_target_r_ - 1) / 2;

    current_hyperedge_radius_     = static_cast<Double>(config_.r);
    current_hyperedge_pdm_radius_ = (std::cosh(current_hyperedge_radius_) - 1) / 2;

    clique_thres_ = 0;

    total_annuli_        = std::floor(alpha_ * target_r_ / std::log(2));
    SInt chunks_per_pe   = config_.k / size_;
    SInt leftover_chunks = config_.k % size_;
    local_chunks_        = chunks_per_pe + static_cast<SInt>((SInt)rank_ < leftover_chunks);

    local_chunk_start_ = (local_chunks_ * rank_) + ((SInt)rank_ >= leftover_chunks ? leftover_chunks : 0);
    local_chunk_end_   = local_chunk_start_ + local_chunks_;

    Double phi_per_chunk = 2 * M_PI / config_.k;
    pe_min_phi_          = local_chunk_start_ * phi_per_chunk;
    pe_max_phi_          = local_chunk_end_ * phi_per_chunk;

    annuli_.set_empty_key(total_annuli_ * config_.k);
    chunks_.set_empty_key(config_.k);
    boundaries_.resize(total_annuli_);

    SInt total_cells = 0;
    cells_per_annulus_.resize(total_annuli_, std::numeric_limits<SInt>::max());

    for (SInt i = 0; i < total_annuli_; ++i) {
        global_cell_ids_.push_back(total_cells);
        total_cells += GridSizeForAnnulus(i) * size_;
    }

    cells_.set_empty_key(total_cells + 1);
    vertices_.set_empty_key(total_cells + 1);

    chunk_eps_ = phi_per_chunk / EPSILON_SCALE<Double>;
    cell_eps_  = (2 * M_PI / GridSizeForAnnulus(total_annuli_ - 1)) / EPSILON_SCALE<Double>;
    point_eps_ = std::numeric_limits<Double>::epsilon();

    num_nodes_ = 0;
}

template <typename Double>
void Hyper_Hyperbolic<Double>::ComputeChunk(const SInt chunk_id) {
    ComputeChunk(chunk_id, config_.n, config_.k, 0, 2 * M_PI, 0, 1, 0);
}

template <typename Double>
void Hyper_Hyperbolic<Double>::ComputeChunk(
    const SInt chunk_id, const SInt n, const SInt num_of_chunks, const Double min_phi, const Double max_phi,
    const SInt chunk_start, const SInt level, const SInt offset) {
    if (num_of_chunks == 1) {
        chunks_[chunk_id] = std::make_tuple(n, min_phi, max_phi, offset);
        return;
    }

    SInt midk = (num_of_chunks + 1) / 2;

    SInt hash_value       = sampling::Spooky::hash(config_.seed + (level * config_.k) + chunk_start);
    SInt splitter_variate = rng_.GenerateBinomial(hash_value, n, (Double)midk / num_of_chunks);

    Double middlePhi = ((max_phi - min_phi) * ((Double)midk / num_of_chunks)) + min_phi;

    if (ZERO_TOLERANCE < middlePhi && middlePhi <= 0.0) {
        middlePhi = 0;
    }

    if (chunk_id < chunk_start + midk) {
        ComputeChunk(chunk_id, splitter_variate, midk, min_phi, middlePhi, chunk_start, level + 1, offset);
    } else {
        ComputeChunk(
            chunk_id, n - splitter_variate, num_of_chunks - midk, middlePhi, max_phi, chunk_start + midk, level + 1,
            offset + splitter_variate);
    }
}

template <typename Double>
void Hyper_Hyperbolic<Double>::ComputeAnnuli(const SInt chunk_id) {
    SInt size   = std::get<0>(chunks_[chunk_id]);
    SInt offset = std::get<3>(chunks_[chunk_id]);

    Double min_r      = 0;
    Double total_area = PGGeometry<Double>::RadiusToHyperbolicArea(alpha_ * target_r_);

    for (SInt i = 1; i < total_annuli_ + 1; i++) {
        Double max_r     = i * target_r_ / total_annuli_;
        Double ring_area = PGGeometry<Double>::RadiusToHyperbolicArea(alpha_ * max_r)
                           - PGGeometry<Double>::RadiusToHyperbolicArea(alpha_ * min_r);

        SInt hash_value =
            sampling::Spooky::hash(config_.seed + (total_annuli_ * config_.k) + (chunk_id * total_annuli_) + i);

        SInt n_annulus =
            rng_.GenerateBinomial(hash_value, size, std::clamp(ring_area / total_area, Double{0.0}, Double{1.0}));

        annuli_[ComputeGlobalChunkId(i - 1, chunk_id)] = std::make_tuple(n_annulus, min_r, max_r, false, offset);
        boundaries_[i - 1]                             = std::make_pair(std::cosh(min_r), std::sinh(min_r));

        min_r = max_r;
        size -= n_annulus;
        offset += n_annulus;
        total_area -= ring_area;
    }
}

template <typename Double>
void Hyper_Hyperbolic<Double>::GenerateCells(const SInt annulus_id, SInt chunk_id) {
    bool clique = false;

    if (chunks_.find(chunk_id) == std::end(chunks_)) {
        ComputeChunk(chunk_id);
        ComputeAnnuli(chunk_id);
    }

    auto& chunk   = chunks_[chunk_id];
    auto& annulus = annuli_[ComputeGlobalChunkId(annulus_id, chunk_id)];

    if (std::get<3>(annulus)) {
        return;
    }

    SInt   size;
    SInt   offset;
    SInt   seed = 0;
    Double min_phi;
    Double max_phi;

    if (clique) {
        size    = std::get<0>(annulus);
        offset  = std::get<3>(annulus);
        min_phi = 0.0;
        max_phi = 2 * M_PI;
        seed    = config_.seed + annulus_id + config_.n;
    } else {
        size    = std::get<0>(annulus);
        offset  = std::get<4>(annulus);
        min_phi = std::get<1>(chunk);
        max_phi = std::get<2>(chunk);
    }

    if (ZERO_TOLERANCE < min_phi && min_phi <= 0.0) {
        min_phi = 0;
    }

    Double total_phi = max_phi - min_phi;
    Double grid_phi  = total_phi / GridSizeForAnnulus(annulus_id);

    for (SInt i = 0; i < GridSizeForAnnulus(annulus_id); ++i) {
        if (!clique) {
            seed = config_.seed + (annulus_id * config_.k) + chunk_id + i + config_.n;
        }

        SInt hash_value = sampling::Spooky::hash(seed);

        SInt n_cell =
            rng_.GenerateBinomial(hash_value, size, std::clamp(grid_phi / total_phi, Double{0.0}, Double{1.0}));

        SInt global_cell_id = ComputeGlobalCellId(annulus_id, chunk_id, i);

        cells_[global_cell_id] =
            std::make_tuple(n_cell, min_phi + (grid_phi * i), min_phi + (grid_phi * (i + 1)), false, offset);

        size -= n_cell;
        offset += n_cell;
        total_phi -= grid_phi;
    }

    std::get<3>(annulus) = true;
}

template <typename Double>
void Hyper_Hyperbolic<Double>::GenerateVertices(const SInt annulus_id, SInt chunk_id, const SInt cell_id) {
    bool clique = false;

    if (chunks_.find(chunk_id) == std::end(chunks_)) {
        ComputeChunk(chunk_id);
        ComputeAnnuli(chunk_id);
    }

    auto& annulus = annuli_[ComputeGlobalChunkId(annulus_id, chunk_id)];

    if (!std::get<3>(annulus)) {
        GenerateCells(annulus_id, chunk_id);
    }

    SInt  global_cell_id = ComputeGlobalCellId(annulus_id, chunk_id, cell_id);
    auto& cell           = cells_[global_cell_id];

    if (std::get<3>(cell)) {
        return;
    }

    SInt   size    = std::get<0>(cell);
    SInt   offset  = std::get<4>(cell);
    Double min_phi = std::get<1>(cell);
    Double max_phi = std::get<2>(cell);
    Double min_r   = std::get<1>(annulus);
    Double max_r   = std::get<2>(annulus);

    SInt seed = 0;

    if (clique) {
        seed = config_.seed + (annulus_id * config_.k * GridSizeForAnnulus(annulus_id));
    } else {
        seed = config_.seed + (annulus_id * config_.k * GridSizeForAnnulus(annulus_id))
               + (chunk_id * GridSizeForAnnulus(annulus_id)) + cell_id + config_.n;
    }

    SInt hash_value = sampling::Spooky::hash(seed);
    mersenne.RandomInit(hash_value);
    sorted_mersenne.RandomInit(hash_value, size);

    const Double mincdf = std::cosh(alpha_ * min_r);
    const Double maxcdf = std::cosh(alpha_ * max_r);

    std::vector<Vertex>& cell_vertices = vertices_[global_cell_id];
    cell_vertices.reserve(size);

    for (SInt i = 0; i < size; i++) {
        Double angle  = (sorted_mersenne.Random() * (max_phi - min_phi)) + min_phi;
        Double radius = std::acosh((mersenne.Random() * (maxcdf - mincdf)) + mincdf) / alpha_;

        Double inv_len    = (std::cosh(radius) + 1.0) / 2.0;
        Double pdm_radius = std::sqrt(1.0 - (1.0 / inv_len));
        Double x          = pdm_radius * std::sin(angle);
        Double y          = pdm_radius * std::cos(angle);
        Double gamma      = 1.0 / (1.0 - (pdm_radius * pdm_radius));

        cell_vertices.emplace_back(angle, radius, x, y, gamma, offset + i);

        if (pe_min_phi_ <= angle && pe_max_phi_ > angle) {
            num_nodes_++;

            if (config_.coordinates) {
                PushCoordinate(x, y);
            }
        }
    }

    std::get<3>(cell) = true;
}

template <typename Double>
inline std::pair<Double, Double> Hyper_Hyperbolic<Double>::GetBoundaryPhis(
    const Double boundary_phi, const Double boundary_r, const SInt annulus_id) const {
    const auto&  boundary    = boundaries_[annulus_id];
    const Double cosh_radius = std::cosh(current_hyperedge_radius_);
    const Double cosh_min_r  = std::get<0>(boundary);
    const Double sinh_min_r  = std::get<1>(boundary);

    Double arg = ((std::cosh(boundary_r) * cosh_min_r) - cosh_radius) / (std::sinh(boundary_r) * sinh_min_r);
    arg        = std::clamp(arg, Double{-1.0}, Double{1.0});

    const Double diff        = std::acos(arg);
    const Double lower_bound = boundary_phi - diff;
    const Double upper_bound = boundary_phi + diff;
    const Double min_phi     = std::min(lower_bound, upper_bound);
    const Double max_phi     = std::max(lower_bound, upper_bound);

    return std::make_pair(min_phi, max_phi);
}

template <typename Double>
inline SInt Hyper_Hyperbolic<Double>::ComputeGlobalCellId(const SInt annulus, const SInt chunk, const SInt cell) {
    return global_cell_ids_[annulus] + (chunk * GridSizeForAnnulus(annulus)) + cell;
}

template <typename Double>
inline SInt Hyper_Hyperbolic<Double>::GridSizeForAnnulus(const SInt annulus_id) {
    return std::max<SInt>(1, TotalGridSizeForAnnulus(annulus_id) / size_);
}

template <typename Double>
void Hyper_Hyperbolic<Double>::GenerateCSR() {
    for (SInt i = local_chunk_start_; i < local_chunk_end_; ++i) {
        ComputeChunk(i);
        ComputeAnnuli(i);
    }

    for (SInt i = local_chunk_start_; i < local_chunk_end_; ++i) {
        for (SInt j = 0; j < total_annuli_; ++j) {
            GenerateCells(j, i);

            for (SInt k = 0; k < GridSizeForAnnulus(j); ++k) {
                GenerateVertices(j, i, k);
            }
        }
    }

    for (SInt i = local_chunk_start_; i < local_chunk_end_; ++i) {
        for (SInt j = 0; j < total_annuli_; ++j) {
            GenerateHyperedges(j, i);
        }
    }

    const SInt start_node = std::get<3>(chunks_[local_chunk_start_]);
    SetVertexRange(start_node, start_node + num_nodes_);
}

template <typename Double>
void Hyper_Hyperbolic<Double>::BeginHyperedge(
    const HyperbolicHyperedgeCenter<Double>& /*center*/, const SInt sampled_center_id) {
    current_hyperedge_pins_.clear();
    current_hyperedge_ranges_.clear();

    current_hyperedge_radius_ = static_cast<Double>(getSampledOrConstantRadius(config_, sampled_center_id)) * target_r_;

    current_hyperedge_pdm_radius_ = (std::cosh(current_hyperedge_radius_) - 1.0) / 2.0;
}

template <typename Double>
bool Hyper_Hyperbolic<Double>::EndHyperedge() {
    auto normalized =
        NormalizeCurrentHyperedge(std::move(current_hyperedge_pins_), std::move(current_hyperedge_ranges_), 4);

    current_hyperedge_pins_   = std::move(normalized.first);
    current_hyperedge_ranges_ = std::move(normalized.second);

    SInt range_pin_count = 0;

    for (const PinRange& range: current_hyperedge_ranges_) {
        range_pin_count += range.end - range.begin;
    }

    if (current_hyperedge_pins_.size() + range_pin_count >= 2) {
        PushHyperedgeCompressed(current_hyperedge_pins_, current_hyperedge_ranges_);
        return true;
    }
    return false;
}

template <typename Double>
void Hyper_Hyperbolic<Double>::PushHyperedgeCompressed(
    const std::vector<SInt>& pins, const std::vector<PinRange>& ranges) {
    if (graph_.hyperedge_offsets.empty()) {
        graph_.hyperedge_offsets.push_back(0);
    }

    if (graph_.hyperedge_range_offsets.empty()) {
        graph_.hyperedge_range_offsets.push_back(0);
    }

    graph_.hyperedge_pins.insert(graph_.hyperedge_pins.end(), pins.begin(), pins.end());

    graph_.hyperedge_offsets.push_back(graph_.hyperedge_pins.size());

    graph_.hyperedge_ranges.insert(graph_.hyperedge_ranges.end(), ranges.begin(), ranges.end());

    graph_.hyperedge_range_offsets.push_back(graph_.hyperedge_ranges.size());
}

template <typename Double>
void Hyper_Hyperbolic<Double>::GenerateHyperedges(const SInt annulus_id, const SInt chunk_id) {
    current_annulus_ = annulus_id;
    current_chunk_   = chunk_id;

    const SInt total_cells = global_cell_ids_.back() + (GridSizeForAnnulus(total_annuli_ - 1) * size_ * config_.k);

    for (SInt cell_id = 0; cell_id < GridSizeForAnnulus(annulus_id); ++cell_id) {
        current_cell_ = cell_id;

        const SInt global_cell_id = ComputeGlobalCellId(annulus_id, chunk_id, cell_id);

        const SInt base_m    = config_.m / total_cells;
        const SInt remainder = config_.m % total_cells;

        const SInt cell_m = base_m + static_cast<SInt>(global_cell_id < remainder);

        const SInt first_sample_id = global_cell_id * base_m + std::min<SInt>(global_cell_id, remainder);

        auto& annulus = annuli_[ComputeGlobalChunkId(annulus_id, chunk_id)];
        auto& cell    = cells_[global_cell_id];

        const Double min_r   = std::get<1>(annulus);
        const Double max_r   = std::get<2>(annulus);
        const Double min_phi = std::get<1>(cell);
        const Double max_phi = std::get<2>(cell);

        const Double mincdf = std::cosh(alpha_ * min_r);
        const Double maxcdf = std::cosh(alpha_ * max_r);

        // Sample attempts
        SInt emitted  = 0;
        SInt attempts = 0;

        const SInt max_attempts = std::max<SInt>(1000, 1000 * cell_m);

        while (emitted < cell_m && attempts < max_attempts) {
            const SInt sampled_center_id = (global_cell_id * max_attempts) + attempts;

            const SInt seed_phi = sampling::Spooky::hash(config_.seed + (17 * sampled_center_id));
            const SInt seed_r   = sampling::Spooky::hash(config_.seed + (31 * sampled_center_id));

            const Double u_phi = rng_.GenerateUniform<Double>(seed_phi);
            const Double u_r   = rng_.GenerateUniform<Double>(seed_r);

            const HyperbolicHyperedgeCenter<Double> center{
                .phi = min_phi + (u_phi * (max_phi - min_phi)),
                .r   = std::acosh((u_r * (maxcdf - mincdf)) + mincdf) / alpha_};

            BeginHyperedge(center, sampled_center_id);
            QueryHyperBallBoth(annulus_id, chunk_id, cell_id, center);

            if (EndHyperedge()) {
                ++emitted;
            }

            ++attempts;
        }
    }
}

template <typename Double>
void Hyper_Hyperbolic<Double>::QueryHyperBallBoth(
    const SInt annulus_id, const SInt chunk_id, const SInt cell_id, const HyperbolicHyperedgeCenter<Double>& center) {
    QueryHyperBall(annulus_id, chunk_id, cell_id, center);

    if (config_.query_both && annulus_id > 0) {
        auto&  chunk         = chunks_[chunk_id];
        Double min_chunk_phi = std::get<1>(chunk);
        Double max_chunk_phi = std::get<2>(chunk);
        Double grid_phi      = (max_chunk_phi - min_chunk_phi) / GridSizeForAnnulus(annulus_id - 1);

        SInt next_cell_id = std::floor((center.phi - min_chunk_phi) / grid_phi);

        QueryHyperBall(annulus_id - 1, chunk_id, next_cell_id, center, false);
    }
}

template <typename Double>
void Hyper_Hyperbolic<Double>::QueryHyperBall(
    const SInt annulus_id, const SInt chunk_id, const SInt cell_id, const HyperbolicHyperedgeCenter<Double>& center,
    bool search_down) {
    auto& annulus        = annuli_[ComputeGlobalChunkId(annulus_id, chunk_id)];
    auto  current_bounds = GetBoundaryPhis(center.phi, center.r, annulus_id);

    current_min_phi_ = std::get<0>(current_bounds);
    current_max_phi_ = std::get<1>(current_bounds);

    Double min_cell_phi = std::get<1>(cells_[ComputeGlobalCellId(annulus_id, chunk_id, cell_id)]);
    Double max_cell_phi = std::get<2>(cells_[ComputeGlobalCellId(annulus_id, chunk_id, cell_id)]);

    right_processed_chunk_ = chunk_id;
    right_processed_cell_  = cell_id;

    if (search_down || !IsLocalChunk(chunk_id)) {
        GenerateGridPins(annulus_id, chunk_id, cell_id, center);
    }

    bool found_nonlocal_chunk = false;

    if (std::get<1>(annulus) >= clique_thres_ && std::max(TotalGridSizeForAnnulus(annulus_id), config_.k) > 1) {
        if (current_min_phi_ < min_cell_phi
            || (OutOfBounds(current_min_phi_) && !(std::abs(min_cell_phi - 0.0) < cell_eps_))) {
            SInt next_chunk_id = chunk_id;

            if (cell_id == 0) {
                next_chunk_id = (chunk_id + config_.k - 1) % config_.k;
            }

            SInt next_cell_id = (cell_id + GridSizeForAnnulus(annulus_id) - 1) % GridSizeForAnnulus(annulus_id);

            GenerateVertices(annulus_id, next_chunk_id, next_cell_id);

            found_nonlocal_chunk |= QueryHyperBallRightNeighbor(
                annulus_id, next_chunk_id, next_cell_id, center, std::abs(min_cell_phi - 0.0) < cell_eps_, search_down);
        }

        if (current_max_phi_ > max_cell_phi
            || (OutOfBounds(current_max_phi_) && !(std::abs(max_cell_phi - (2 * M_PI)) < cell_eps_))) {
            SInt next_chunk_id = chunk_id;

            if (cell_id == GridSizeForAnnulus(annulus_id) - 1) {
                next_chunk_id = (chunk_id + config_.k + 1) % config_.k;
            }

            SInt next_cell_id = (cell_id + GridSizeForAnnulus(annulus_id) + 1) % GridSizeForAnnulus(annulus_id);

            GenerateVertices(annulus_id, next_chunk_id, next_cell_id);

            found_nonlocal_chunk |= QueryHyperBallLeftNeighbor(
                annulus_id, next_chunk_id, next_cell_id, center, std::abs(max_cell_phi - (2 * M_PI)) < cell_eps_,
                search_down);
        }
    }

    SInt next_annulus;

    if (search_down) {
        next_annulus = annulus_id + 1;
    } else {
        next_annulus = annulus_id - 1;
    }

    if (next_annulus >= total_annuli_ || (SSInt)next_annulus < 0) {
        return;
    }

    auto&  chunk         = chunks_[chunk_id];
    Double min_chunk_phi = std::get<1>(chunk);
    Double max_chunk_phi = std::get<2>(chunk);
    Double grid_phi      = (max_chunk_phi - min_chunk_phi) / GridSizeForAnnulus(next_annulus);

    SInt next_cell_id = std::floor((center.phi - min_chunk_phi) / grid_phi);

    QueryHyperBall(next_annulus, chunk_id, next_cell_id, center, search_down);
}

template <typename Double>
bool Hyper_Hyperbolic<Double>::QueryHyperBallRightNeighbor(
    const SInt annulus_id, SInt chunk_id, SInt cell_id, const HyperbolicHyperedgeCenter<Double>& center, bool phase,
    bool search_down) {
    bool found_nonlocal_chunk = false;

    while (true) {
        if (phase && current_min_phi_ < 0.0) {
            current_min_phi_ += 2 * M_PI;
        }

        if (phase && OutOfBounds(current_min_phi_)) {
            return found_nonlocal_chunk;
        }

        auto&  cell            = cells_[ComputeGlobalCellId(annulus_id, chunk_id, cell_id)];
        Double min_cell_phi    = std::get<1>(cell);
        right_processed_chunk_ = chunk_id;
        right_processed_cell_  = cell_id;

        if (search_down || !IsLocalChunk(chunk_id)) {
            if (!IsLocalChunk(chunk_id)) {
                found_nonlocal_chunk = true;
            }

            GenerateGridPins(annulus_id, chunk_id, cell_id, center);
        }

        phase = phase || std::abs(min_cell_phi - 0.0) < cell_eps_;

        if (current_min_phi_ < min_cell_phi || OutOfBounds(current_min_phi_)) {
            SInt next_chunk_id = chunk_id;

            if (cell_id == 0) {
                next_chunk_id = (chunk_id + config_.k - 1) % config_.k;
            }

            SInt next_cell_id = (cell_id + GridSizeForAnnulus(annulus_id) - 1) % GridSizeForAnnulus(annulus_id);

            GenerateVertices(annulus_id, next_chunk_id, next_cell_id);

            cell_id  = next_cell_id;
            chunk_id = next_chunk_id;

            continue;
        }

        return found_nonlocal_chunk;
    }
}

template <typename Double>
bool Hyper_Hyperbolic<Double>::QueryHyperBallLeftNeighbor(
    const SInt annulus_id, SInt chunk_id, SInt cell_id, const HyperbolicHyperedgeCenter<Double>& center, bool phase,
    bool search_down) {
    bool found_nonlocal_chunk = false;

    while (true) {
        if (phase && current_max_phi_ >= 2 * M_PI) {
            current_max_phi_ -= 2 * M_PI;
        }

        if (phase && OutOfBounds(current_max_phi_)) {
            return found_nonlocal_chunk;
        }

        if (chunk_id == right_processed_chunk_ && cell_id == right_processed_cell_) {
            return found_nonlocal_chunk;
        }

        auto&  cell         = cells_[ComputeGlobalCellId(annulus_id, chunk_id, cell_id)];
        Double max_cell_phi = std::get<2>(cell);

        if (search_down || !IsLocalChunk(chunk_id)) {
            if (!search_down) {
                found_nonlocal_chunk = true;
            }

            GenerateGridPins(annulus_id, chunk_id, cell_id, center);
        }

        phase = phase || std::abs(max_cell_phi - (2 * M_PI)) < cell_eps_;

        if (current_max_phi_ > max_cell_phi || OutOfBounds(current_max_phi_)) {
            SInt next_chunk_id = chunk_id;

            if (cell_id == GridSizeForAnnulus(annulus_id) - 1) {
                next_chunk_id = (chunk_id + config_.k + 1) % config_.k;
            }

            SInt next_cell_id = (cell_id + GridSizeForAnnulus(annulus_id) + 1) % GridSizeForAnnulus(annulus_id);

            GenerateVertices(annulus_id, next_chunk_id, next_cell_id);

            cell_id  = next_cell_id;
            chunk_id = next_chunk_id;

            continue;
        }

        return found_nonlocal_chunk;
    }
}

template <typename Double>
void Hyper_Hyperbolic<Double>::AddWholeCellPins(SInt global_cell_id) {
    const auto& cell = cells_[global_cell_id];

    const SInt cell_size     = std::get<0>(cell);
    const SInt vertex_offset = std::get<4>(cell);

    if (cell_size > 0) {
        current_hyperedge_ranges_.push_back({vertex_offset, vertex_offset + cell_size});
    }
}

template <typename Double>
void Hyper_Hyperbolic<Double>::GenerateGridPins(
    const SInt annulus_id, const SInt chunk_id, const SInt cell_id, const HyperbolicHyperedgeCenter<Double>& center) {
    const SInt  global_cell_id = ComputeGlobalCellId(annulus_id, chunk_id, cell_id);
    const auto& cell           = cells_[global_cell_id];

    if (std::get<0>(cell) == 0) {
        return;
    }

    const auto relation = ClassifyCellAgainstHyperball(annulus_id, chunk_id, cell_id, center);

    if (relation == CellBallRelation::OUTSIDE) {
        return;
    }

    if (relation == CellBallRelation::INSIDE) {
        AddWholeCellPins(global_cell_id);
        return;
    }

    GenerateVertices(annulus_id, chunk_id, cell_id);

    const std::vector<Vertex>& cell_vertices = vertices_[global_cell_id];

    for (const Vertex& vertex_to_check: cell_vertices) {
        if (HyperbolicDistance(center, vertex_to_check) <= current_hyperedge_pdm_radius_) {
            current_hyperedge_pins_.push_back(std::get<5>(vertex_to_check));
        }
    }
}

template <typename Double>
CellBallRelation Hyper_Hyperbolic<Double>::ClassifyCellAgainstHyperball(
    const SInt annulus_id, const SInt chunk_id, const SInt cell_id,
    const HyperbolicHyperedgeCenter<Double>& hyperball_center) {
    const auto& annulus = annuli_.find(ComputeGlobalChunkId(annulus_id, chunk_id))->second;
    const auto& cell    = cells_.find(ComputeGlobalCellId(annulus_id, chunk_id, cell_id))->second;

    const Double min_r   = std::get<1>(annulus);
    const Double max_r   = std::get<2>(annulus);
    const Double min_phi = std::get<1>(cell);
    const Double max_phi = std::get<2>(cell);

    const Double q_phi = hyperball_center.phi;
    const Double q_r   = hyperball_center.r;

    const bool all_corners_inside =
        IsPointInsideHyperball(q_r, q_phi, min_r, min_phi) && IsPointInsideHyperball(q_r, q_phi, min_r, max_phi)
        && IsPointInsideHyperball(q_r, q_phi, max_r, min_phi) && IsPointInsideHyperball(q_r, q_phi, max_r, max_phi);

    if (all_corners_inside) {
        return CellBallRelation::INSIDE;
    }

    return CellBallRelation::PARTIAL;
}

template <typename Double>
bool Hyper_Hyperbolic<Double>::IsPointInsideHyperball(
    const Double center_r, const Double center_phi, const Double r, const Double phi) const {
    Double delta_phi = std::abs(center_phi - phi);
    delta_phi        = std::min(delta_phi, Double{2.0 * M_PI} - delta_phi);

    const Double cosh_dist =
        (std::cosh(center_r) * std::cosh(r)) - (std::sinh(center_r) * std::sinh(r) * std::cos(delta_phi));

    return cosh_dist <= std::cosh(current_hyperedge_radius_);
}

template <typename Double>
Double
Hyper_Hyperbolic<Double>::HyperbolicDistance(const HyperbolicHyperedgeCenter<Double>& center, const Vertex& vertex) {
    Double delta_phi = std::abs(center.phi - std::get<0>(vertex));
    delta_phi        = std::min(delta_phi, Double{2.0 * M_PI} - delta_phi);

    const Double cosh_dist = std::cosh(center.r) * std::cosh(std::get<1>(vertex))
                             - std::sinh(center.r) * std::sinh(std::get<1>(vertex)) * std::cos(delta_phi);

    return (cosh_dist - 1.0) / 2.0;
}

template <typename Double>
SInt Hyper_Hyperbolic<Double>::TotalGridSizeForAnnulus(const SInt annulus_id) {
    if (cells_per_annulus_[annulus_id] != std::numeric_limits<SInt>::max()) {
        return cells_per_annulus_[annulus_id];
    }

    Double min_r      = annulus_id * target_r_ / total_annuli_;
    Double max_r      = (annulus_id + 1) * target_r_ / total_annuli_;
    Double ring_area  = PGGeometry<Double>::RadiusToHyperbolicArea(alpha_ * max_r)
                        - PGGeometry<Double>::RadiusToHyperbolicArea(alpha_ * min_r);
    Double total_area = PGGeometry<Double>::RadiusToHyperbolicArea(alpha_ * target_r_);

    SInt exp_points = config_.n * ring_area / total_area;
    SInt cells      = exp_points / config_.hyp_base;

    SInt result                    = std::max<SInt>(1, cells);
    cells_per_annulus_[annulus_id] = result;

    return result;
}

template <typename Double>
inline bool Hyper_Hyperbolic<Double>::OutOfBounds(const Double num) const {
    return std::isnan(num) || num < -2 * M_PI || num > 2 * M_PI;
}

template <typename Double>
inline SInt Hyper_Hyperbolic<Double>::ComputeGlobalChunkId(const SInt annulus, const SInt chunk) const {
    return (chunk * total_annuli_) + annulus;
}

template <typename Double>
inline bool Hyper_Hyperbolic<Double>::IsLocalChunk(const SInt chunk_id) const {
    return chunk_id >= local_chunk_start_ && chunk_id < local_chunk_end_;
}

template class Hyper_Hyperbolic<LPFloat>;
template class Hyper_Hyperbolic<HPFloat>;

} // namespace kagen