#include "kagen/generators/hyper/h_geometric/hyper_rgg2d.h"

#include "kagen/generators/generator.h"
#include "kagen/generators/hyper/h_geometric/hyper_rgg2d_policy.h"
#include "kagen/hypergraph/hyperedge_builder.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>
namespace kagen {

HyperRGG2D::HyperRGG2D(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : SpatialGrid2D(config, rank, size) {
    if (config_.debug) {
        debug_logger_.emplace(MakeDebugFilename(), true);
    }
}

HyperRGG2D::~HyperRGG2D() = default;

void HyperRGG2D::GenerateCells(const SInt chunk_id) {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    ++weak_scaling_instrumentation_.generate_cells_requests;
#endif

    if (chunks_.find(chunk_id) == chunks_.end()) {
        ComputeChunk(chunk_id);
    }

    const auto chunk_it = chunks_.find(chunk_id);

    if (chunk_it == chunks_.end()) {
        return;
    }

    auto& chunk = chunk_it->second;

    if (std::get<3>(chunk)) {
        return;
    }

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    const std::size_t records_before = cells_.size();
    const double      start          = MPI_Wtime();
#endif

    GenerateCellMetadataSubtree(chunk_id, 0, cells_per_chunk_, std::get<0>(chunk), std::get<4>(chunk), 1);

    std::get<3>(chunk) = true;

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    weak_scaling_instrumentation_.chunk_materialization_seconds += MPI_Wtime() - start;

    ++weak_scaling_instrumentation_.chunks_materialized;

    weak_scaling_instrumentation_.cell_slots_scanned += static_cast<std::uint64_t>(cells_per_chunk_);

    weak_scaling_instrumentation_.cell_records_inserted += static_cast<std::uint64_t>(cells_.size() - records_before);
#endif
}

#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION

void HyperRGG2D::PrintWeakScalingInstrumentation(MPI_Comm comm) const {
    int rank = 0;
    int size = 1;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    auto reduce_counter = [&](const char* name, const std::uint64_t local) {
        std::uint64_t global_sum = 0;
        std::uint64_t global_max = 0;

        MPI_Reduce(&local, &global_sum, 1, MPI_UINT64_T, MPI_SUM, 0, comm);
        MPI_Reduce(&local, &global_max, 1, MPI_UINT64_T, MPI_MAX, 0, comm);

        if (rank == 0) {
            const double mean_per_rank = static_cast<double>(global_sum) / static_cast<double>(size);

            std::cout << "HRGG2D_PROFILE"
                      << ";metric=" << name << ";sum=" << global_sum << ";mean_per_rank=" << mean_per_rank
                      << ";max_per_rank=" << global_max << '\n';
        }
    };

    reduce_counter("generate_cells_requests", weak_scaling_instrumentation_.generate_cells_requests);

    reduce_counter("chunks_materialized", weak_scaling_instrumentation_.chunks_materialized);

    reduce_counter("cell_slots_scanned", weak_scaling_instrumentation_.cell_slots_scanned);

    reduce_counter("cell_records_inserted", weak_scaling_instrumentation_.cell_records_inserted);

    const double local_seconds = weak_scaling_instrumentation_.chunk_materialization_seconds;

    double minimum_seconds = 0.0;
    double sum_seconds     = 0.0;
    double maximum_seconds = 0.0;

    MPI_Reduce(&local_seconds, &minimum_seconds, 1, MPI_DOUBLE, MPI_MIN, 0, comm);

    MPI_Reduce(&local_seconds, &sum_seconds, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

    MPI_Reduce(&local_seconds, &maximum_seconds, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    if (rank == 0) {
        const double mean_seconds = sum_seconds / static_cast<double>(size);

        std::cout << "HRGG2D_PROFILE"
                  << ";metric=chunk_materialization_seconds"
                  << ";min=" << minimum_seconds << ";mean=" << mean_seconds << ";max=" << maximum_seconds << '\n';
    }
}

#endif

std::string HyperRGG2D::MakeDebugFilename() const {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::string output = config_.output_graph.filename + "_" + std::to_string(config_.n) + "_"
                         + std::to_string(config_.m) + "_" + std::to_string(config_.hyperedge_radius_exponent)
                         + "_debug_rank_" + std::to_string(rank) + ".csv";
    return output;
}

SInt HyperRGG2D::SafeTotalCellsPerDim() const {
    if (cells_per_dim_ <= 0 || chunks_per_dim_ <= 0) {
        throw ConfigurationError("Invalid grid: cells_per_dim_ or chunks_per_dim_ is zero");
    }

    if (chunks_per_dim_ > std::numeric_limits<SInt>::max() / cells_per_dim_) {
        throw ConfigurationError("Invalid grid: total_cells_per_dim overflow");
    }

    return chunks_per_dim_ * cells_per_dim_;
}

HyperRGG2D::CellPosition HyperRGG2D::GetCellPosition(
    const SInt chunk_row, const SInt chunk_column, const SInt cell_row, const SInt cell_column,
    const SInt total_cells_per_dim) const {
    const SInt cell_id = EncodeCell(cell_column, cell_row);

    const SInt global_cell_x  = (chunk_column * cells_per_dim_) + cell_column;
    const SInt global_cell_y  = (chunk_row * cells_per_dim_) + cell_row;
    const SInt global_cell_id = (global_cell_y * total_cells_per_dim) + global_cell_x;

    return {
        .cell_id        = cell_id,
        .global_cell_x  = global_cell_x,
        .global_cell_y  = global_cell_y,
        .global_cell_id = global_cell_id,
    };
}

void HyperRGG2D::GenerateEdges(const SInt chunk_row, const SInt chunk_column) {
    const SInt chunk_id = Encode(chunk_column, chunk_row);

    if (!IsLocalChunk(chunk_id)) {
        return;
    }

    const SInt total_cells_per_dim = SafeTotalCellsPerDim();

    if (total_cells_per_dim > std::numeric_limits<SInt>::max() / total_cells_per_dim) {
        throw ConfigurationError("Invalid grid: total_cells overflow");
    }

    const SInt total_cells = total_cells_per_dim * total_cells_per_dim;

    if (!policy_) {
        throw std::logic_error("HRGG2D policy was not initialized");
    }

    HyperedgeBuilder<HyperRGG2DPolicy> builder(*policy_, config_, debug_logger_ ? &*debug_logger_ : nullptr);

    LPFloat lower_bound = static_cast<LPFloat>(policy_->MinimumRadius({}));
    LPFloat upper_bound = LPFloat{1.0};

    if (config_.min_hyperedge_radius != -1.0) {
        lower_bound = static_cast<LPFloat>(config_.min_hyperedge_radius);
    }

    if (config_.max_hyperedge_radius != -1.0) {
        upper_bound = static_cast<LPFloat>(config_.max_hyperedge_radius);
    }

    const SInt    base_m     = config_.m / total_cells;
    const SInt    remainder  = config_.m % total_cells;
    const LPFloat cell_width = LPFloat{1.0} / static_cast<LPFloat>(total_cells_per_dim);

    CenterSampler sampler(*this, cell_width, lower_bound, upper_bound);

    constexpr SInt dump_cell_x = 10;
    constexpr SInt dump_cell_y = 10;

    for (SInt cell_row = 0; cell_row < cells_per_dim_; ++cell_row) {
        for (SInt cell_column = 0; cell_column < cells_per_dim_; ++cell_column) {
            const auto cell = GetCellPosition(chunk_row, chunk_column, cell_row, cell_column, total_cells_per_dim);

            const bool dump_this_cell = cell.global_cell_x == dump_cell_x && cell.global_cell_y == dump_cell_y;

            const SInt cell_m = base_m + static_cast<SInt>(cell.global_cell_id < remainder);

            // ----------------------------------------------------
            // Dump vertices belonging to this cell
            // ----------------------------------------------------

            if (dump_this_cell) {
                std::vector<Vertex> vertices;

                GenerateVertices(chunk_id, cell.cell_id, vertices);

                std::ofstream out("/tmp/cell_vertices.csv");

                out << std::setprecision(17);

                out << "vertex_id,x,y\n";

                for (const auto& vertex: vertices) {
                    out << std::get<2>(vertex) << "," << std::get<0>(vertex) << "," << std::get<1>(vertex) << "\n";
                }

                std::cout << "[geometry dump] dumped " << vertices.size() << " vertices from cell (" << dump_cell_x
                          << ", " << dump_cell_y << ")\n";
            }

            // ----------------------------------------------------
            // Generate centers normally
            // ----------------------------------------------------

            std::ofstream centers_out;

            if (dump_this_cell) {
                centers_out.open("/tmp/cell_centers.csv");
                centers_out << std::setprecision(17);

                centers_out << "sampled_id,x,y,radius\n";
            }

            for (SInt emitted = 0; emitted < cell_m; emitted++) {
                const Center center = sampler.Sample(chunk_id, cell, emitted, base_m, remainder);

                if (dump_this_cell) {
                    centers_out << center.sampled_id << "," << center.x << "," << center.y << "," << center.radius
                                << "\n";
                }

                // IMPORTANT:
                // use the already sampled center.
                // Do NOT call sampler.Sample() again.
                builder.Build(center);
            }
        }
    }
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    policy_->PrintExactCacheStats();
#endif
}

void HyperRGG2D::ObserveHyperedgeAndMaybeReserve(std::size_t pins, std::size_t ranges) {
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

void HyperRGG2D::PushHyperedgeCompressed(const std::vector<SInt>& pins, const std::vector<PinRange>& ranges) {
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

void HyperRGG2D::GenerateCSR() {
    const std::size_t expected_local_hyperedges = (static_cast<std::size_t>(config_.m) + size_ - 1) / size_;

    graph_.hyperedge_offsets.reserve(expected_local_hyperedges + 1);
    graph_.hyperedge_range_offsets.reserve(expected_local_hyperedges + 1);

    if (config_.size_dist_pin_budget > 0) {
        const std::size_t expected_local_pins =
            (static_cast<std::size_t>(config_.size_dist_pin_budget) + size_ - 1) / size_;

        graph_.hyperedge_pins.reserve(expected_local_pins);
    }
    policy_ = std::make_unique<HyperRGG2DPolicy>(*this);

    GenerateGeometry();

    policy_.reset();
    if (config_.debug) {
        std::cerr << "Cells per dimension: " << cells_per_dim_ << '\n';
    }
}

void HyperRGG2D::FinalizeCSR(MPI_Comm comm) {
#ifdef KAGEN_ENABLE_HYPER_INSTRUMENTATION
    PrintWeakScalingInstrumentation(comm);
#endif
    if (!config_.debug) {
        return;
    }

    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    const SInt local_pins                 = static_cast<SInt>(graph_.hyperedge_pins.size());
    const SInt local_max_hyperedge_pins   = local_memory_stats_.max_pins_per_hyperedge;
    const SInt local_max_hyperedge_ranges = local_memory_stats_.max_ranges_per_hyperedge;

    SInt global_max_pins             = 0;
    SInt global_sum_pins             = 0;
    SInt global_max_hyperedge_pins   = 0;
    SInt global_max_hyperedge_ranges = 0;

    MPI_Reduce(&local_pins, &global_max_pins, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, comm);
    MPI_Reduce(&local_pins, &global_sum_pins, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
    MPI_Reduce(&local_max_hyperedge_pins, &global_max_hyperedge_pins, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, comm);
    MPI_Reduce(&local_max_hyperedge_ranges, &global_max_hyperedge_ranges, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, comm);

    if (rank == 0) {
        std::cerr << "[HyperRGG2D memory] "
                  << "sum_pins=" << global_sum_pins << " max_local_pins=" << global_max_pins
                  << " max_hyperedge_pins=" << global_max_hyperedge_pins
                  << " max_hyperedge_ranges=" << global_max_hyperedge_ranges << '\n';
    }
}
HyperRGG2D::CenterSampler::CenterSampler(HyperRGG2D& gen, LPFloat cell_width, LPFloat lower_bound, LPFloat upper_bound)
    : gen_(gen),
      cell_width_(cell_width),
      lower_bound_(lower_bound),
      upper_bound_(upper_bound) {}

HyperRGG2D::Center
HyperRGG2D::CenterSampler::Sample(SInt chunk_id, const CellPosition& cell, SInt emitted, SInt base_m, SInt remainder) {
    const SInt sampled_center_id = (cell.global_cell_id * base_m) + std::min(cell.global_cell_id, remainder) + emitted;

    const SInt center_seed = sampling::Spooky::hash(gen_.config_.seed + (17 * sampled_center_id));

    gen_.mersenne.RandomInit(center_seed);

    const LPFloat random_x = static_cast<LPFloat>(gen_.mersenne.Random());
    const LPFloat random_y = static_cast<LPFloat>(gen_.mersenne.Random());

    const LPFloat center_x = (static_cast<LPFloat>(cell.global_cell_x) + random_x) * cell_width_;
    const LPFloat center_y = (static_cast<LPFloat>(cell.global_cell_y) + random_y) * cell_width_;

    const LPFloat radius =
        gen_.config_.random_radius
            ? static_cast<LPFloat>(SampleHyperedgeRadius(gen_.config_, lower_bound_, upper_bound_, gen_.mersenne))
            : static_cast<LPFloat>(gen_.config_.r);

    return {
        .x          = center_x,
        .y          = center_y,
        .radius     = radius,
        .sampled_id = sampled_center_id,
        .chunk_id   = chunk_id,
        .cell_id    = cell.cell_id,
    };
}

SInt HyperRGG2D::SampleLeftCellCount(
    const SInt chunk_id, const SInt tree_node, const SInt vertices, const SInt left_cells, const SInt total_cells) {
    if (vertices == 0 || left_cells == 0) {
        return 0;
    }

    if (left_cells == total_cells) {
        return vertices;
    }

    const SInt seed = sampling::Spooky::hash(
        config_.seed ^ sampling::Spooky::hash(chunk_id + 0x9e3779b97f4a7c15ULL)
        ^ sampling::Spooky::hash(tree_node + 0xd1b54a32d192ed03ULL));

    const LPFloat probability = static_cast<LPFloat>(left_cells) / static_cast<LPFloat>(total_cells);

    return rng_.GenerateBinomial(seed, vertices, std::clamp(probability, LPFloat{0.0}, LPFloat{1.0}));
}

HyperRGG2D::CellMetadata HyperRGG2D::ReconstructCellMetadata(const SInt chunk_id, const SInt cell_id) {
    if (cell_id >= cells_per_chunk_) {
        throw std::out_of_range("HRGG2D cell ID is out of range");
    }

    if (chunks_.find(chunk_id) == chunks_.end()) {
        ComputeChunk(chunk_id);
    }

    const auto chunk_it = chunks_.find(chunk_id);

    if (chunk_it == chunks_.end()) {
        return {
            .size    = 0,
            .offset  = 0,
            .start_x = 0.0,
            .start_y = 0.0,
        };
    }

    const auto& chunk = chunk_it->second;

    SInt begin     = 0;
    SInt end       = cells_per_chunk_;
    SInt vertices  = std::get<0>(chunk);
    SInt offset    = std::get<4>(chunk);
    SInt tree_node = 1;

    while (end - begin > 1) {
        const SInt middle      = begin + ((end - begin) / 2);
        const SInt left_cells  = middle - begin;
        const SInt total_cells = end - begin;

        const SInt left_vertices = SampleLeftCellCount(chunk_id, tree_node, vertices, left_cells, total_cells);

        if (cell_id < middle) {
            end       = middle;
            vertices  = left_vertices;
            tree_node = 2 * tree_node;
        } else {
            begin = middle;
            offset += left_vertices;
            vertices -= left_vertices;
            tree_node = (2 * tree_node) + 1;
        }
    }

    SInt cell_x = 0;
    SInt cell_y = 0;

    DecodeCell(cell_id, cell_x, cell_y);

    return {
        .size    = vertices,
        .offset  = offset,
        .start_x = std::get<1>(chunk) + (static_cast<LPFloat>(cell_x) * cell_size_),
        .start_y = std::get<2>(chunk) + (static_cast<LPFloat>(cell_y) * cell_size_),
    };
}

void HyperRGG2D::GenerateCellMetadataSubtree(
    const SInt chunk_id, const SInt begin_cell, const SInt end_cell, const SInt vertices, const SInt offset,
    const SInt tree_node) {
    if (end_cell - begin_cell == 1) {
        if (vertices == 0) {
            return;
        }

        SInt cell_x = 0;
        SInt cell_y = 0;

        DecodeCell(begin_cell, cell_x, cell_y);

        const auto& chunk = chunks_[chunk_id];

        cells_[ComputeGlobalCellId(chunk_id, begin_cell)] = std::make_tuple(
            vertices, std::get<1>(chunk) + (static_cast<LPFloat>(cell_x) * cell_size_),
            std::get<2>(chunk) + (static_cast<LPFloat>(cell_y) * cell_size_), false, offset);

        return;
    }

    const SInt middle      = begin_cell + ((end_cell - begin_cell) / 2);
    const SInt left_cells  = middle - begin_cell;
    const SInt total_cells = end_cell - begin_cell;

    const SInt left_vertices = SampleLeftCellCount(chunk_id, tree_node, vertices, left_cells, total_cells);

    GenerateCellMetadataSubtree(chunk_id, begin_cell, middle, left_vertices, offset, 2 * tree_node);

    GenerateCellMetadataSubtree(
        chunk_id, middle, end_cell, vertices - left_vertices, offset + left_vertices, (2 * tree_node) + 1);
}

void HyperRGG2D::GenerateRemoteCellVertices(const SInt chunk_id, const SInt cell_id, std::vector<Vertex>& buffer) {
    const CellMetadata cell = ReconstructCellMetadata(chunk_id, cell_id);

    buffer.clear();
    buffer.reserve(cell.size);

    const SInt seed = config_.seed + (chunk_id * cells_per_chunk_) + cell_id;

    const SInt hash = sampling::Spooky::hash(seed);

    mersenne.RandomInit(hash);

    for (SInt i = 0; i < cell.size; ++i) {
        const LPFloat x = (mersenne.Random() * cell_size_) + cell.start_x;

        const LPFloat y = (mersenne.Random() * cell_size_) + cell.start_y;

        buffer.emplace_back(x, y, cell.offset + i);
    }
}

} // namespace kagen
