#include "kagen/generators/hyper/h_hyperbolic/hyper_hyperbolic.h"

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/generators/hyper/h_hyperbolic/hyper_hyperbolic_policy.h"
#include "kagen/generators/hyper/h_hyperbolic/poincare_geometry.h"
#include "kagen/hypergraph/hyperedge_builder.h"
#include "kagen/hypergraph/hyperedge_management.h"
#include "kagen/hypergraph/hypergraph_utils.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"
#include "kagen/tools/geometry.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>
namespace kagen {

constexpr double M_N_D_RATIO                     = 2.0;
constexpr double HIGH_PRECISION_THRESHOLD_LOG2_N = 29.0;

template <typename Double>
constexpr Double EPSILON_SCALE = 1000.0;

PGeneratorConfig
Hyper_HyperbolicFactory::NormalizeParameters(PGeneratorConfig config, PEID rank, PEID size, const bool output) const {
    config.setChunkSizeIfMissing(size);

    if (config.k < static_cast<SInt>(size)) {
        throw ConfigurationError("Number of chunks must be at least the number of PEs");
    }

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

    if (config.size_dist_pin_budget > 0) {
        if (config.n <= 0 || config.m <= 0) {
            throw ConfigurationError("expected total pins requires n > 0 and m > 0");
        }

        if (!config.random_radius) {
            throw ConfigurationError("expected total pins requires random_radius=true");
        }

        config.hyperedge_radius_exponent = SolveHyperbolicRadiusExponentForExpectedPins(config);

        if (rank == 0) {
            std::cout << " Chosen radius exponent = " << config.hyperedge_radius_exponent << '\n';
        }
    }

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

template <typename Double>
void Hyper_Hyperbolic<Double>::PrecomputeMinimumRadii() {
    minimum_radii_by_center_annulus_.resize(total_annuli_);

    constexpr Double desired_pins = Double{2.0};

    for (SInt a = 0; a < total_annuli_; ++a) {
        const Double min_r = a * target_r_ / total_annuli_;
        const Double max_r = (a + 1) * target_r_ / total_annuli_;
        const Double mid_r = (min_r + max_r) / Double{2.0};

        const HyperbolicHyperedgeCenter<Double> center{
            .phi        = Double{0.0},
            .r          = mid_r,
            .sampled_id = 0,
            .annulus_id = a,
        };

        minimum_radii_by_center_annulus_[a] = FindRadiusForExpectedPins(center, desired_pins);
    }
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
      rng_(config),
      alpha_((config_.plexp - 1) / 2),
      target_r_(PGGeometry<Double>::GetTargetRadius(config_.n, config_.n * config_.avg_degree / 2, alpha_)),
      cosh_target_r_(std::cosh(target_r_)),
      pdm_target_r_((cosh_target_r_ - 1) / 2),
      clique_thres_(0),
      total_annuli_(std::floor(alpha_ * target_r_ / std::numbers::ln2)),
      current_hyperedge_radius_(static_cast<Double>(config_.r)),
      current_hyperedge_pdm_radius_((std::cosh(current_hyperedge_radius_) - 1) / 2) {
    const SInt chunks_per_pe   = config_.k / size_;
    const SInt leftover_chunks = config_.k % size_;

    local_chunks_ = chunks_per_pe + static_cast<SInt>(static_cast<SInt>(rank_) < leftover_chunks);

    local_chunk_start_ =
        (static_cast<SInt>(rank_) * chunks_per_pe) + std::min<SInt>(static_cast<SInt>(rank_), leftover_chunks);

    local_chunk_end_ = local_chunk_start_ + local_chunks_;

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
        total_cells += GridSizeForAnnulus(i) * config_.k;
    }

    cells_.set_empty_key(total_cells + 1);
    vertices_.set_empty_key(total_cells + 1);
    center_cells_.set_empty_key(total_cells + 1);

    chunk_eps_ = phi_per_chunk / EPSILON_SCALE<Double>;
    cell_eps_  = (2 * M_PI / GridSizeForAnnulus(total_annuli_ - 1)) / EPSILON_SCALE<Double>;
    point_eps_ = std::numeric_limits<Double>::epsilon();

    num_nodes_ = 0;

    if (config_.random_radius && config_.min_hyperedge_radius == -1.0) {
        PrecomputeMinimumRadii();
    }

    annulus_min_r_.resize(total_annuli_);
    annulus_max_r_.resize(total_annuli_);
    annulus_min_cosh_.resize(total_annuli_);
    annulus_min_sinh_.resize(total_annuli_);
    annulus_max_cosh_.resize(total_annuli_);
    annulus_max_sinh_.resize(total_annuli_);

    for (SInt a = 0; a < total_annuli_; ++a) {
        const Double min_r = a * target_r_ / total_annuli_;
        const Double max_r = (a + 1) * target_r_ / total_annuli_;

        annulus_min_r_[a] = min_r;
        annulus_max_r_[a] = max_r;

        annulus_min_cosh_[a] = std::cosh(min_r);
        annulus_min_sinh_[a] = std::sinh(min_r);
        annulus_max_cosh_[a] = std::cosh(max_r);
        annulus_max_sinh_[a] = std::sinh(max_r);
    }

    if (config_.debug) {
        debug_logger_.emplace(MakeDebugFilename(), true);
    }
}
template <typename Double>
std::string Hyper_Hyperbolic<Double>::MakeDebugFilename() const {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::string output = config_.output_graph.filename + "_" + std::to_string(config_.n) + "_"
                         + std::to_string(config_.m) + "_" + std::to_string(config_.hyperedge_radius_exponent)
                         + "_debug_rank_" + std::to_string(rank) + ".csv";
    return output;
}

template <typename Double>
Double Hyper_Hyperbolic<Double>::FindRadiusForExpectedPins(
    const HyperbolicHyperedgeCenter<Double>& center, const Double desired_pins) {
    auto expected_pins = [&](const Double radius) {
        Double expected = 0.0;

        const Double total_area = PGGeometry<Double>::RadiusToHyperbolicArea(alpha_ * target_r_);

        for (SInt a = 0; a < total_annuli_; ++a) {
            const Double min_r = a * target_r_ / total_annuli_;
            const Double max_r = (a + 1) * target_r_ / total_annuli_;

            const Double ring_area = PGGeometry<Double>::RadiusToHyperbolicArea(alpha_ * max_r)
                                     - PGGeometry<Double>::RadiusToHyperbolicArea(alpha_ * min_r);

            const Double expected_annulus_n = static_cast<Double>(config_.n) * ring_area / total_area;

            const Double query_r = std::clamp(center.r, min_r, max_r);
            const Double denom   = std::sinh(center.r) * std::sinh(query_r);

            Double angular_reach;

            if (denom <= Double{0.0}) {
                angular_reach = M_PI;
            } else {
                const Double x = ((std::cosh(center.r) * std::cosh(query_r)) - std::cosh(radius)) / denom;

                if (x <= Double{-1.0}) {
                    angular_reach = M_PI;
                } else if (x >= Double{1.0}) {
                    angular_reach = Double{0.0};
                } else {
                    angular_reach = std::acos(x);
                }
            }

            const Double angular_fraction = std::min<Double>(Double{1.0}, angular_reach / Double{M_PI});

            expected += expected_annulus_n * angular_fraction;
        }

        return expected;
    };

    Double lo = Double{0.0};
    Double hi = center.r + target_r_;

    for (SInt i = 0; i < 40; ++i) {
        const Double mid = (lo + hi) / Double{2.0};

        if (expected_pins(mid) >= desired_pins) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    return hi;
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
SInt Hyper_Hyperbolic<Double>::VertexCellSeed(SInt annulus_id, SInt chunk_id, SInt cell_id) {
    return config_.seed + (annulus_id * config_.k * GridSizeForAnnulus(annulus_id))
           + (chunk_id * GridSizeForAnnulus(annulus_id)) + cell_id + config_.n;
}

template <typename Double>
Hyper_Hyperbolic<Double>::SampledVertex
Hyper_Hyperbolic<Double>::SampleVertex(Double min_phi, Double max_phi, Double min_cdf, Double max_cdf) {
    const Double phi = (sorted_mersenne.Random() * (max_phi - min_phi)) + min_phi;
    const Double r   = std::acosh((mersenne.Random() * (max_cdf - min_cdf)) + min_cdf) / alpha_;

    const Double cosh_r = std::cosh(r);
    const Double sinh_r = std::sinh(r);

    const Double inv_len    = (cosh_r + Double{1.0}) / Double{2.0};
    const Double pdm_radius = std::sqrt(Double{1.0} - (Double{1.0} / inv_len));
    const Double sin_phi    = std::sin(phi);
    const Double cos_phi    = std::cos(phi);

    return {
        .r       = r,
        .phi     = phi,
        .cosh_r  = cosh_r,
        .sinh_r  = sinh_r,
        .cos_phi = cos_phi,
        .sin_phi = sin_phi,
        .x       = pdm_radius * sin_phi,
        .y       = pdm_radius * cos_phi,
        .gamma   = cosh_r,
    };
}
template <typename Double>
void Hyper_Hyperbolic<Double>::AppendVertex(
    VertexBlock& block, const SInt id, const Hyper_Hyperbolic<Double>::SampledVertex& vertex) {
    block.r.push_back(vertex.r);
    block.id.push_back(id);

    block.x.push_back(vertex.x);
    block.y.push_back(vertex.y);
    block.gamma.push_back(vertex.gamma);

    block.cosh_r.push_back(vertex.cosh_r);
    block.sinh_r.push_back(vertex.sinh_r);
    block.cos_phi.push_back(vertex.cos_phi);
    block.sin_phi.push_back(vertex.sin_phi);
    block.phi.push_back(vertex.phi);
}

template <typename Double>
void Hyper_Hyperbolic<Double>::GenerateVertices(const SInt annulus_id, SInt chunk_id, const SInt cell_id) {
    const SInt global_cell_id = ComputeGlobalCellId(annulus_id, chunk_id, cell_id);
    auto&      cell           = cells_[global_cell_id];

    if (std::get<3>(cell)) {
        return;
    }

    VertexBlock& cell_vertices = vertices_[global_cell_id];
    GenerateVertices(annulus_id, chunk_id, cell_id, cell_vertices);

    std::get<3>(cell) = true;
}

template <typename Double>
void Hyper_Hyperbolic<Double>::GenerateVertices(
    const SInt annulus_id, SInt chunk_id, const SInt cell_id, VertexBlock& out) {
    if (chunks_.find(chunk_id) == std::end(chunks_)) {
        ComputeChunk(chunk_id);
        ComputeAnnuli(chunk_id);
    }

    auto& annulus = annuli_[ComputeGlobalChunkId(annulus_id, chunk_id)];

    if (!std::get<3>(annulus)) {
        GenerateCells(annulus_id, chunk_id);
    }

    const SInt global_cell_id = ComputeGlobalCellId(annulus_id, chunk_id, cell_id);
    auto&      cell           = cells_[global_cell_id];

    out.clear();

    const SInt   size    = std::get<0>(cell);
    const SInt   offset  = std::get<4>(cell);
    const Double min_phi = std::get<1>(cell);
    const Double max_phi = std::get<2>(cell);
    const Double min_r   = std::get<1>(annulus);
    const Double max_r   = std::get<2>(annulus);

    const SInt seed = VertexCellSeed(annulus_id, chunk_id, cell_id);

    const SInt hash_value = sampling::Spooky::hash(seed);
    mersenne.RandomInit(hash_value);
    sorted_mersenne.RandomInit(hash_value, size);

    const Double mincdf = std::cosh(alpha_ * min_r);
    const Double maxcdf = std::cosh(alpha_ * max_r);

    out.reserve(size);

    for (SInt i = 0; i < size; ++i) {
        const auto vertex = SampleVertex(min_phi, max_phi, mincdf, maxcdf);
        AppendVertex(out, offset + i, vertex);

        if (config_.coordinates && pe_min_phi_ <= vertex.phi && vertex.phi < pe_max_phi_) {
            PushCoordinate(vertex.x, vertex.y);
        }
    }
}

template <typename Double>
inline SInt Hyper_Hyperbolic<Double>::ComputeGlobalCellId(const SInt annulus, const SInt chunk, const SInt cell) {
    return global_cell_ids_[annulus] + (chunk * GridSizeForAnnulus(annulus)) + cell;
}

template <typename Double>
inline SInt Hyper_Hyperbolic<Double>::GridSizeForAnnulus(const SInt annulus_id) {
    return std::max<SInt>(1, TotalGridSizeForAnnulus(annulus_id) / config_.k);
}

template <typename Double>
void Hyper_Hyperbolic<Double>::GenerateCSR() {
    graph_.hyperedge_offsets.reserve(config_.m + 1);
    graph_.hyperedge_range_offsets.reserve(config_.m + 1);
    for (SInt i = local_chunk_start_; i < local_chunk_end_; ++i) {
        ComputeChunk(i);
        ComputeAnnuli(i);
    }

    for (SInt i = local_chunk_start_; i < local_chunk_end_; ++i) {
        num_nodes_ += std::get<0>(chunks_[i]);
    }
    const SInt start_node = std::get<3>(chunks_[local_chunk_start_]);
    SetVertexRange(start_node, start_node + num_nodes_);

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
    HyperbolicGeometryPolicy<Double>                   geometry(*this, 0, 0);
    HyperedgeBuilder<HyperbolicGeometryPolicy<Double>> builder(
        geometry, config_, debug_logger_ ? &*debug_logger_ : nullptr);

    for (SInt i = local_chunk_start_; i < local_chunk_end_; ++i) {
        for (SInt j = 0; j < total_annuli_; ++j) {
            GenerateHyperedges(j, i, builder);
        }
    }

    if (config_.debug) {
        geometry.PrintExactCacheStats();
    }
}

template <typename Double>
void Hyper_Hyperbolic<Double>::BeginHyperedge(const HyperbolicHyperedgeCenter<Double>& center, Mersenne& mersenne) {
    current_hyperedge_pins_.clear();
    current_hyperedge_ranges_.clear();

    current_hyperedge_radius_ = Radius(center, mersenne);

    current_hyperedge_pdm_radius_ = (std::cosh(current_hyperedge_radius_) - 1.0) / 2.0;
}

template <typename Double>
void Hyper_Hyperbolic<Double>::EndHyperedge() {
    auto normalized =
        NormalizeCurrentHyperedge(std::move(current_hyperedge_pins_), std::move(current_hyperedge_ranges_), 4);

    current_hyperedge_pins_   = std::move(normalized.first);
    current_hyperedge_ranges_ = std::move(normalized.second);

    PushHyperedgeCompressed(current_hyperedge_pins_, current_hyperedge_ranges_);
}

template <typename Double>
void Hyper_Hyperbolic<Double>::ObserveHyperedgeAndMaybeReserve(std::size_t pins, std::size_t ranges) {
    const std::size_t needed_pins = graph_.hyperedge_pins.size() + pins;
    if (needed_pins > graph_.hyperedge_pins.capacity()) {
        const std::size_t slack = std::min<std::size_t>(needed_pins / 8, 1'000'000);
        graph_.hyperedge_pins.reserve(needed_pins + slack);
    }

    const std::size_t needed_ranges = graph_.hyperedge_ranges.size() + ranges;
    if (needed_ranges > graph_.hyperedge_ranges.capacity()) {
        const std::size_t slack = std::min<std::size_t>(needed_ranges / 8, 100'000);
        graph_.hyperedge_ranges.reserve(needed_ranges + slack);
    }
}

template <typename Double>
void kagen::Hyper_Hyperbolic<Double>::PushHyperedgeCompressed(
    const std::vector<SInt>& pins, const std::vector<PinRange>& ranges) {
    local_memory_stats_.max_pins_per_hyperedge =
        std::max(local_memory_stats_.max_pins_per_hyperedge, static_cast<SInt>(pins.size()));

    local_memory_stats_.max_ranges_per_hyperedge =
        std::max(local_memory_stats_.max_ranges_per_hyperedge, static_cast<SInt>(ranges.size()));

    ObserveHyperedgeAndMaybeReserve(pins.size(), ranges.size());
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

    const SInt hash_value = sampling::Spooky::hash(config_.seed + seed_offset + (level * config_.k) + chunk_start);

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
            config_.seed + seed_offset + (total_annuli_ * config_.k) + (chunk_id * total_annuli_) + i);

        const SInt n_annulus =
            rng_.GenerateBinomial(hash_value, size, std::clamp(ring_area / total_area, Double{0.0}, Double{1.0}));

        annuli[ComputeGlobalChunkId(i - 1, chunk_id)] = std::make_tuple(n_annulus, min_r, max_r, false, offset);

        if (seed_offset == 0) {
            boundaries_[i - 1] = std::make_pair(std::cosh(min_r), std::sinh(min_r));
        }

        min_r = max_r;
        size -= n_annulus;
        offset += n_annulus;
        total_area -= ring_area;
    }
}

template <typename Double>
template <typename ChunkMap, typename AnnulusMap, typename CellMap>
void Hyper_Hyperbolic<Double>::GenerateCellsInto(
    const SInt annulus_id, SInt chunk_id, ChunkMap& chunks, AnnulusMap& annuli, CellMap& cells,
    const SInt seed_offset) {
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
        const SInt seed = config_.seed + seed_offset + (annulus_id * config_.k) + chunk_id + i + size;

        const SInt hash_value = sampling::Spooky::hash(seed);

        const SInt n_cell =
            rng_.GenerateBinomial(hash_value, size, std::clamp(grid_phi / total_phi, Double{0.0}, Double{1.0}));

        const SInt global_cell_id = ComputeGlobalCellId(annulus_id, chunk_id, i);

        if constexpr (std::is_same_v<typename CellMap::mapped_type, Cell>) {
            const Double cell_min_phi = min_phi + (grid_phi * i);
            const Double cell_max_phi = min_phi + (grid_phi * (i + 1));

            const Double annulus_min_r = std::get<1>(annulus);
            const Double annulus_max_r = std::get<2>(annulus);

            const auto box =
                poincare_geometry::ComputeCellAABB(annulus_min_r, annulus_max_r, cell_min_phi, cell_max_phi, cell_eps_);

            cells[global_cell_id] = std::make_tuple(
                n_cell, cell_min_phi, cell_max_phi, false, offset, box.min_x, box.max_x, box.min_y, box.max_y);
        } else {
            cells[global_cell_id] =
                std::make_tuple(n_cell, min_phi + (grid_phi * i), min_phi + (grid_phi * (i + 1)), false, offset);
        }
        size -= n_cell;
        offset += n_cell;
        total_phi -= grid_phi;
    }

    std::get<3>(annulus) = true;
}
template <typename Double>
void Hyper_Hyperbolic<Double>::SeedHyperedgeRNG(const SInt sampled_center_id) {
    const SInt seed = sampling::Spooky::hash(config_.seed + sampled_center_id);
    mersenne.RandomInit(seed);
}

template <typename Double>
HyperbolicHyperedgeCenter<Double>
Hyper_Hyperbolic<Double>::SampleCenter(SInt annulus_id, SInt sampled_center_id, const CenterSamplingRegion& region) {
    const Double u_phi = mersenne.Random();
    const Double u_r   = mersenne.Random();

    return {
        .phi        = region.min_phi + (u_phi * (region.max_phi - region.min_phi)),
        .r          = std::acosh((u_r * (region.max_cdf - region.min_cdf)) + region.min_cdf) / alpha_,
        .sampled_id = sampled_center_id,
        .annulus_id = annulus_id,
    };
}

template <typename Double>
void Hyper_Hyperbolic<Double>::GenerateHyperedges(
    const SInt annulus_id, const SInt chunk_id, HyperedgeBuilder<HyperbolicGeometryPolicy<Double>>& builder) {
    current_annulus_ = annulus_id;
    current_chunk_   = chunk_id;

    GenerateCenterCells(annulus_id, chunk_id);

    for (SInt cell_id = 0; cell_id < GridSizeForAnnulus(annulus_id); ++cell_id) {
        current_cell_ = cell_id;

        const SInt global_cell_id = ComputeGlobalCellId(annulus_id, chunk_id, cell_id);

        auto& center_annulus = center_annuli_[ComputeGlobalChunkId(annulus_id, chunk_id)];
        auto& center_cell    = center_cells_[global_cell_id];

        const auto region = BuildCenterSamplingRegion(center_annulus, center_cell);

        if (region.count == 0) {
            continue;
        }

        for (SInt emitted = 0; emitted < region.count; ++emitted) {
            const SInt sampled_center_id = region.offset + emitted;

            SeedHyperedgeRNG(sampled_center_id);

            const auto center = SampleCenter(annulus_id, sampled_center_id, region);

            BeginHyperedge(center, mersenne);
            builder.Build(center);
        }
    }
}

template <typename Double>
Hyper_Hyperbolic<Double>::CenterSamplingRegion Hyper_Hyperbolic<Double>::BuildCenterSamplingRegion(
    const CenterAnnulus& center_annulus, const CenterCell& center_cell) const {
    const Double min_r = std::get<1>(center_annulus);
    const Double max_r = std::get<2>(center_annulus);

    return {
        .min_phi = std::get<1>(center_cell),
        .max_phi = std::get<2>(center_cell),
        .min_cdf = std::cosh(alpha_ * min_r),
        .max_cdf = std::cosh(alpha_ * max_r),
        .offset  = std::get<4>(center_cell),
        .count   = std::get<0>(center_cell),
    };
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

template <typename Double>
Double Hyper_Hyperbolic<Double>::MinimumRadius(const HyperbolicHyperedgeCenter<Double>& center) const {
    if (config_.min_hyperedge_radius != -1.0) {
        return static_cast<Double>(config_.min_hyperedge_radius);
    }

    return minimum_radii_by_center_annulus_[center.annulus_id];
}

template <typename Double>
Double Hyper_Hyperbolic<Double>::MaximumRadius(const HyperbolicHyperedgeCenter<Double>& center) const {
    if (config_.max_hyperedge_radius != -1.0) {
        return static_cast<Double>(config_.max_hyperedge_radius);
    }

    return center.r + target_r_;
}

template <typename Double>
Double Hyper_Hyperbolic<Double>::Radius(const HyperbolicHyperedgeCenter<Double>& center, Mersenne& mersenne) {
    if (!config_.random_radius) {
        return static_cast<Double>(config_.r);
    }

    const Double lower = MinimumRadius(center);
    const Double upper = std::max(lower, MaximumRadius(center));

    return static_cast<Double>(
        SampleHyperedgeRadius(config_, static_cast<double>(lower), static_cast<double>(upper), mersenne));
}

template class Hyper_Hyperbolic<LPFloat>;
template class Hyper_Hyperbolic<HPFloat>;

} // namespace kagen