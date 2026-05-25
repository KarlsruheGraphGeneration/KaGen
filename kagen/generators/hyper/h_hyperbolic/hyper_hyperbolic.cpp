#include "kagen/generators/hyper/h_hyperbolic/hyper_hyperbolic.h"

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/generators/hyper/h_hyperbolic/hyper_hyperbolic_policy.h"
#include "kagen/hypergraph/hyperedge_builder.h"
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

template <typename Double>
constexpr Double EPSILON_SCALE = 1000.0;

constexpr double ZERO_TOLERANCE = 1e-8;

PGeneratorConfig Hyper_HyperbolicFactory::NormalizeParameters(
    PGeneratorConfig config, PEID /*rank*/, PEID size, const bool output) const {
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

    center_chunks_.set_empty_key(config_.k);
    center_annuli_.set_empty_key(total_annuli_ * config_.k);

    SInt total_cells = 0;
    cells_per_annulus_.resize(total_annuli_, std::numeric_limits<SInt>::max());

    for (SInt i = 0; i < total_annuli_; ++i) {
        global_cell_ids_.push_back(total_cells);
        total_cells += GridSizeForAnnulus(i) * size_;
    }

    cells_.set_empty_key(total_cells + 1);
    vertices_.set_empty_key(total_cells + 1);
    center_cells_.set_empty_key(total_cells + 1);

    chunk_eps_ = phi_per_chunk / EPSILON_SCALE<Double>;
    cell_eps_  = (2 * M_PI / GridSizeForAnnulus(total_annuli_ - 1)) / EPSILON_SCALE<Double>;
    point_eps_ = std::numeric_limits<Double>::epsilon();

    num_nodes_ = 0;
}
template <typename Double>
void Hyper_Hyperbolic<Double>::ComputeChunk(const SInt chunk_id) {
    ComputeChunkInto(chunks_, chunk_id, config_.n, config_.k, 0, 2 * M_PI, 0, 1, 0, 0);
}

template <typename Double>
void Hyper_Hyperbolic<Double>::ComputeCenterChunk(const SInt chunk_id) {
    ComputeChunkInto(center_chunks_, chunk_id, config_.m, config_.k, 0, 2 * M_PI, 0, 1, 0, 9001);
}

template <typename Double>
void Hyper_Hyperbolic<Double>::ComputeAnnuli(const SInt chunk_id) {
    ComputeAnnuliInto(chunks_, annuli_, chunk_id, 0);
}

template <typename Double>
void Hyper_Hyperbolic<Double>::ComputeCenterAnnuli(const SInt chunk_id) {
    ComputeAnnuliInto(center_chunks_, center_annuli_, chunk_id, 9101);
}

template <typename Double>
void Hyper_Hyperbolic<Double>::GenerateCells(const SInt annulus_id, SInt chunk_id) {
    if (chunks_.find(chunk_id) == std::end(chunks_)) {
        ComputeChunk(chunk_id);
        ComputeAnnuli(chunk_id);
    }

    GenerateCellsInto(annulus_id, chunk_id, chunks_, annuli_, cells_, 0);
}
template <typename Double>
void Hyper_Hyperbolic<Double>::GenerateCenterCells(const SInt annulus_id, SInt chunk_id) {
    if (center_chunks_.find(chunk_id) == std::end(center_chunks_)) {
        ComputeCenterChunk(chunk_id);
        ComputeCenterAnnuli(chunk_id);
    }

    GenerateCellsInto(annulus_id, chunk_id, center_chunks_, center_annuli_, center_cells_, 9201);
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

            if (config_.partial_cell_mode == PartialCellMode::GenerateAndCheck) {
                for (SInt k = 0; k < GridSizeForAnnulus(j); ++k) {
                    GenerateVertices(j, i, k);
                }
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
    const HyperbolicHyperedgeCenter<Double>& /*center*/, const SInt sampled_center_id, const Double lower_bound) {
    current_hyperedge_pins_.clear();
    current_hyperedge_ranges_.clear();

    current_hyperedge_radius_ =
        static_cast<Double>(getSampledOrConstantRadius(config_, sampled_center_id), lower_bound) * target_r_;

    current_hyperedge_pdm_radius_ = (std::cosh(current_hyperedge_radius_) - 1.0) / 2.0;
}

template <typename Double>
void Hyper_Hyperbolic<Double>::EndHyperedge() {
    auto normalized =
        NormalizeCurrentHyperedge(std::move(current_hyperedge_pins_), std::move(current_hyperedge_ranges_), 4);

    current_hyperedge_pins_   = std::move(normalized.first);
    current_hyperedge_ranges_ = std::move(normalized.second);

    SInt range_pin_count = 0;

    for (const PinRange& range: current_hyperedge_ranges_) {
        range_pin_count += range.end - range.begin;
    }

    PushHyperedgeCompressed(current_hyperedge_pins_, current_hyperedge_ranges_);
}

template <typename Double>
void kagen::Hyper_Hyperbolic<Double>::PushHyperedgeCompressed(
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
template <typename ChunkMap>
void Hyper_Hyperbolic<Double>::ComputeChunkInto(
    ChunkMap& chunks, const SInt chunk_id, const SInt total_objects, const SInt num_chunks, const Double min_phi,
    const Double max_phi, const SInt chunk_start, const SInt level, const SInt offset, const SInt seed_offset) {
    if (num_chunks == 1) {
        chunks[chunk_id] = std::make_tuple(total_objects, min_phi, max_phi, offset);
        return;
    }

    const SInt midk = (num_chunks + 1) / 2;

    const SInt hash_value = sampling::Spooky::hash(config_.seed + seed_offset + level * config_.k + chunk_start);

    const SInt splitter_variate =
        rng_.GenerateBinomial(hash_value, total_objects, static_cast<Double>(midk) / num_chunks);

    const Double middle_phi = ((max_phi - min_phi) * (static_cast<Double>(midk) / num_chunks)) + min_phi;

    if (chunk_id < chunk_start + midk) {
        ComputeChunkInto(
            chunks, chunk_id, splitter_variate, midk, min_phi, middle_phi, chunk_start, level + 1, offset, seed_offset);
    } else {
        ComputeChunkInto(
            chunks, chunk_id, total_objects - splitter_variate, num_chunks - midk, middle_phi, max_phi,
            chunk_start + midk, level + 1, offset + splitter_variate, seed_offset);
    }
}

template <typename Double>
template <typename ChunkMap, typename AnnulusMap>
void Hyper_Hyperbolic<Double>::ComputeAnnuliInto(
    const ChunkMap& chunks, AnnulusMap& annuli, const SInt chunk_id, const SInt seed_offset) {
    const auto& chunk  = chunks.find(chunk_id)->second;
    SInt        size   = std::get<0>(chunk);
    SInt        offset = std::get<3>(chunk);

    Double min_r      = 0;
    Double total_area = PGGeometry<Double>::RadiusToHyperbolicArea(alpha_ * target_r_);

    for (SInt i = 1; i < total_annuli_ + 1; ++i) {
        const Double max_r = i * target_r_ / total_annuli_;

        const Double ring_area = PGGeometry<Double>::RadiusToHyperbolicArea(alpha_ * max_r)
                                 - PGGeometry<Double>::RadiusToHyperbolicArea(alpha_ * min_r);

        const SInt hash_value = sampling::Spooky::hash(
            config_.seed + seed_offset + total_annuli_ * config_.k + chunk_id * total_annuli_ + i);

        const SInt n_annulus =
            rng_.GenerateBinomial(hash_value, size, std::clamp(ring_area / total_area, Double{0.0}, Double{1.0}));

        annuli[ComputeGlobalChunkId(i - 1, chunk_id)] = std::make_tuple(n_annulus, min_r, max_r, false, offset);

        min_r = max_r;
        size -= n_annulus;
        offset += n_annulus;
        total_area -= ring_area;
    }
}
template <typename Double>
Double NormalizePhi(Double phi) {
    const Double two_pi = 2.0 * M_PI;

    while (phi < 0.0) {
        phi += two_pi;
    }

    while (phi >= two_pi) {
        phi -= two_pi;
    }

    return phi;
}

template <typename Double>
Double CircularPhiWidth(Double min_phi, Double max_phi) {
    const Double two_pi = 2.0 * M_PI;

    min_phi = NormalizePhi(min_phi);
    max_phi = NormalizePhi(max_phi);

    if (min_phi <= max_phi) {
        return max_phi - min_phi;
    }

    return (two_pi - min_phi) + max_phi;
}
template <typename Double>
Double HyperbolicCellArea(const Double min_r, const Double max_r, const Double min_phi, const Double max_phi) {
    const Double phi_width = CircularPhiWidth(min_phi, max_phi);

    if (max_r <= min_r || phi_width <= 0.0) {
        return 0.0;
    }

    return phi_width * (std::cosh(max_r) - std::cosh(min_r));
}
template <typename Double>
Double ExpectedRadiusForKVertices(const Double density, const SInt k) {
    if (density <= 0.0 || k == 0) {
        return std::numeric_limits<Double>::infinity();
    }

    const Double arg = 1.0 + static_cast<Double>(k) / (2.0 * M_PI * density);

    return std::acosh(arg);
}

template <typename Double>
template <typename ChunkMap, typename AnnulusMap, typename CellMap>
void Hyper_Hyperbolic<Double>::GenerateCellsInto(
    const SInt annulus_id, SInt chunk_id, ChunkMap& chunks, AnnulusMap& annuli, CellMap& cells,
    const SInt seed_offset) {
    if (chunks.find(chunk_id) == std::end(chunks)) {
        // This part cannot be fully generic unless you pass callbacks,
        // so usually call ComputeChunk/ComputeAnnuli before this.
    }

    auto& chunk   = chunks[chunk_id];
    auto& annulus = annuli[ComputeGlobalChunkId(annulus_id, chunk_id)];

    if (std::get<3>(annulus)) {
        return;
    }

    SInt size   = std::get<0>(annulus);
    SInt offset = std::get<4>(annulus);

    const Double min_phi = std::get<1>(chunk);
    const Double max_phi = std::get<2>(chunk);

    Double       total_phi = max_phi - min_phi;
    const Double grid_phi  = total_phi / GridSizeForAnnulus(annulus_id);

    for (SInt i = 0; i < GridSizeForAnnulus(annulus_id); ++i) {
        const SInt seed = config_.seed + seed_offset + annulus_id * config_.k + chunk_id + i + size;

        const SInt h = sampling::Spooky::hash(seed);

        const SInt n_cell = rng_.GenerateBinomial(h, size, std::clamp(grid_phi / total_phi, Double{0.0}, Double{1.0}));

        cells[ComputeGlobalCellId(annulus_id, chunk_id, i)] =
            std::make_tuple(n_cell, min_phi + grid_phi * i, min_phi + grid_phi * (i + 1), false, offset);

        size -= n_cell;
        offset += n_cell;
        total_phi -= grid_phi;
    }

    std::get<3>(annulus) = true;
}

template <typename Double>
void Hyper_Hyperbolic<Double>::GenerateHyperedges(const SInt annulus_id, const SInt chunk_id) {
    current_annulus_ = annulus_id;
    current_chunk_   = chunk_id;

    GenerateCenterCells(annulus_id, chunk_id);

    HyperbolicGeometryPolicy<Double>                   geometry(*this, annulus_id, chunk_id);
    HyperedgeBuilder<HyperbolicGeometryPolicy<Double>> builder(geometry, config_.partial_cell_mode);

    for (SInt cell_id = 0; cell_id < GridSizeForAnnulus(annulus_id); ++cell_id) {
        current_cell_ = cell_id;
        geometry.SetStartCell(annulus_id, chunk_id, cell_id);

        const SInt global_cell_id = ComputeGlobalCellId(annulus_id, chunk_id, cell_id);

        auto& annulus = center_annuli_[ComputeGlobalChunkId(annulus_id, chunk_id)];
        auto& cell    = center_cells_[global_cell_id];

        const SInt cell_m = std::get<0>(cell);
        if (cell_m == 0) {
            continue;
        }

        const Double min_r         = std::get<1>(annulus);
        const Double max_r         = std::get<2>(annulus);
        const Double min_phi       = std::get<1>(cell);
        const Double max_phi       = std::get<2>(cell);
        const SInt   center_offset = std::get<4>(cell);

        const Double mincdf = std::cosh(alpha_ * min_r);
        const Double maxcdf = std::cosh(alpha_ * max_r);

        const SInt   size        = std::get<0>(cell);
        const Double area        = HyperbolicCellArea(min_r, max_r, min_phi, max_phi);
        const Double density     = area > 0.0 ? static_cast<Double>(size) / area : 0.0;
        const Double lower_bound = ExpectedRadiusForKVertices(density, 2);

        for (SInt emitted = 0; emitted < cell_m; ++emitted) {
            const SInt sampled_center_id = center_offset + emitted;

            const SInt seed_phi = sampling::Spooky::hash(config_.seed + (17 * sampled_center_id));

            const SInt seed_r = sampling::Spooky::hash(config_.seed + (31 * sampled_center_id));

            const Double u_phi = rng_.GenerateUniform<Double>(seed_phi);
            const Double u_r   = rng_.GenerateUniform<Double>(seed_r);

            const HyperbolicHyperedgeCenter<Double> center{
                .phi = min_phi + (u_phi * (max_phi - min_phi)),
                .r   = std::acosh((u_r * (maxcdf - mincdf)) + mincdf) / alpha_};

            BeginHyperedge(center, sampled_center_id, lower_bound);
            builder.Build(center);
        }
    }
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

template <typename Double>
void Hyper_Hyperbolic<Double>::PushHyperedgePin(const SInt pin) {
    current_hyperedge_pins_.push_back(pin);
}

template <typename Double>
void Hyper_Hyperbolic<Double>::PushHyperedgeRange(const SInt begin, const SInt end) {
    if (begin < end) {
        current_hyperedge_ranges_.push_back({begin, end});
    }
}

template class Hyper_Hyperbolic<LPFloat>;
template class Hyper_Hyperbolic<HPFloat>;

} // namespace kagen