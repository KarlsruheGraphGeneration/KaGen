#include "kagen/generators/grid/grid_2d.h"

#include "kagen/generators/generator.h"
#include "kagen/sampling/hash.hpp"

namespace kagen {
PGeneratorConfig
Grid2DFactory::NormalizeParameters(PGeneratorConfig config, PEID, const PEID size, const bool output) const {
    EnsureSquarePowerOfTwoChunkSize(config, size, output);
    if ((config.grid_x == 0 || config.grid_y == 0) && config.n == 0) { 
        throw ConfigurationError("Grid dimension 0 not allowed.");
    } else if (config.grid_x == 0 && config.grid_y == 0 && config.n != 0) {
        const SInt sqrt_n = std::sqrt(config.n);
        const bool is_perfect_square = (sqrt_n * sqrt_n == config.n);
        if (!is_perfect_square && output) {
            std::cerr << "Warning: n = " << config.n << " is not a perfect square; using grid " << sqrt_n << " x " << sqrt_n << " (" << sqrt_n * sqrt_n << " vertices)\n";
        } 
        config.grid_x     = sqrt_n;
        config.grid_y     = sqrt_n;
    } else if (config.n == 0 && config.grid_x > 0 && config.grid_y > 0) {
       config.n = config.grid_x * config.grid_y;
    }
    if (config.n != config.grid_x * config.grid_y) {
        throw ConfigurationError("Dimensions do not fit number of vertices.");
    }
    if (config.p == 0) {
        if (config.m == 0) {
            throw ConfigurationError("if p is not given, m must be nonzero.");
        } 
        
        if (config.grid_x == 1 && config.grid_y == 1) {
            throw ConfigurationError("p is not given, and the resulting graph would have zero edges.");
        }

        // directed
        SInt max_directed_edges = 0;
        
        auto axis_neighbors = [](const SInt L) -> SInt {
            if (L <= 1) return 0;
            if (L == 2) return 1;
            return 2;
        };
        
        if (config.periodic) {
            const SInt deg = axis_neighbors(config.grid_x) + axis_neighbors(config.grid_y);
            max_directed_edges = deg * config.n;
        } else {
            if (config.grid_x == 1 || config.grid_y == 1) {
                    max_directed_edges = (config.grid_x == 1) ? 2 * (config.grid_y - 1) : 2 * (config.grid_x - 1);
            } else {
                const SInt num_deg2_vertices = 4;
                const SInt num_deg3_vertices = 2 * config.grid_x + 2 * config.grid_y - 8;
                const SInt num_deg4_vertices = config.grid_x * config.grid_y - num_deg2_vertices - num_deg3_vertices;
                max_directed_edges = (4 * num_deg4_vertices + 3 * num_deg3_vertices + 2 * num_deg2_vertices);
            }
        }

        if (max_directed_edges % 2 != 0) throw std::logic_error("Sum of degrees (directed edges) must be even.");

        if (max_directed_edges == 0) throw ConfigurationError("If p is not given, the maximum number of possible edges must be non-zero.");

        config.p = 2.0 * config.m / max_directed_edges;
        if (output) {
            std::cout << "Setting edge probability to " << config.p << std::endl;
            if (config.p > 1) {
                throw ConfigurationError("Configuration infeasible, too many edges");
            }
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

    return config;
}

std::unique_ptr<Generator>
Grid2DFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<Grid2D>(config, rank, size);
}

Grid2D::Grid2D(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size),
      rng_(config) {}

void Grid2D::GenerateEdgeList() {
    // Init dimensions
    // @todo Only tested for cube PEs and one chunk per PE
    total_rows_       = config_.grid_y;
    total_cols_       = config_.grid_x;
    edge_probability_ = config_.p;

    // Init chunks
    total_chunks_   = config_.k;
    chunks_per_dim_ = sqrt(total_chunks_);

    SInt leftover_chunks = total_chunks_ % size_;
    SInt num_chunks      = (total_chunks_ / size_) + ((SInt)rank_ < leftover_chunks);
    SInt start_chunk     = rank_ * num_chunks + ((SInt)rank_ >= leftover_chunks ? leftover_chunks : 0);
    SInt end_chunk       = start_chunk + num_chunks;

    // Chunk distribution
    rows_per_chunk_ = total_rows_ / chunks_per_dim_;
    remaining_rows_ = total_rows_ % chunks_per_dim_;

    cols_per_chunk_ = total_cols_ / chunks_per_dim_;
    remaining_cols_ = total_cols_ % chunks_per_dim_;

    vertices_per_chunk_ = rows_per_chunk_ * cols_per_chunk_;

    start_node_ = OffsetForChunk(start_chunk);
    end_node_   = OffsetForChunk(end_chunk);
    num_nodes_  = end_node_ - start_node_;

    for (SInt i = 0; i < num_chunks; i++) {
        GenerateChunk(start_chunk + i);
    }

    if (config_.coordinates) {
        for (SInt y = 0; y < total_rows_; ++y) {
            for (SInt x = 0; x < total_cols_; ++x) {
                PushCoordinate(1.0 * x / total_cols_, 1.0 * y / total_rows_);
            }
        }
    }

    SetVertexRange(start_node_, start_node_ + num_nodes_);
}

void Grid2D::GenerateChunk(const SInt chunk) {
    SInt start_vertex = OffsetForChunk(chunk);
    SInt end_vertex   = OffsetForChunk(chunk + 1);
    for (SInt i = start_vertex; i < end_vertex; i++) {
        GenerateEdges(chunk, i);
    }
}

void Grid2D::GenerateEdges(const SInt chunk, const SInt vertex) {
    QueryInDirection(chunk, vertex, Direction::Right);
    QueryInDirection(chunk, vertex, Direction::Left);
    QueryInDirection(chunk, vertex, Direction::Up);
    QueryInDirection(chunk, vertex, Direction::Down);
}

void Grid2D::QueryInDirection(const SInt chunk, const SInt vertex, Direction direction) {
    SInt offset       = OffsetForChunk(chunk);
    SInt local_vertex = vertex - offset;

    SInt chunk_row, chunk_col;
    Decode(chunk, chunk_row, chunk_col);

    SInt rows = rows_per_chunk_ + (chunk_row < remaining_rows_);
    SInt cols = cols_per_chunk_ + (chunk_col < remaining_cols_);

    SInt local_row = local_vertex / cols;
    SInt local_col = local_vertex % cols;

    SInt local_neighbor_row = local_row + DirectionRow(direction);
    SInt local_neighbor_col = local_col + DirectionColumn(direction);

    if (IsLocalVertex(local_neighbor_row, local_neighbor_col, rows, cols)) {
        SInt neighbor_vertex = offset + (local_neighbor_row * cols + local_neighbor_col);
        GenerateEdge(vertex, neighbor_vertex);
    } else {
        // Determine neighboring chunk
        SSInt neighbor_chunk_row = (SSInt)chunk_row + DirectionRow(direction);
        SSInt neighbor_chunk_col = (SSInt)chunk_col + DirectionColumn(direction);
        if (config_.periodic) {
            neighbor_chunk_row = (neighbor_chunk_row + chunks_per_dim_) % chunks_per_dim_;
            neighbor_chunk_col = (neighbor_chunk_col + chunks_per_dim_) % chunks_per_dim_;
        }
        if (!IsValidChunk(neighbor_chunk_row, neighbor_chunk_col))
            return;

        SInt neighbor_chunk  = Encode(neighbor_chunk_row, neighbor_chunk_col);
        SInt neighbor_vertex = LocateVertexInChunk(neighbor_chunk, local_row, local_col, direction);
        GenerateEdge(vertex, neighbor_vertex);
    }
}

bool Grid2D::IsLocalVertex(const SInt local_row, const SInt local_col, const SInt rows, const SInt cols) const {
    if (local_row >= rows)
        return false;
    if (local_col >= cols)
        return false;
    return true;
}

bool Grid2D::IsValidChunk(const SInt chunk_row, const SInt chunk_col) const {
    if (chunk_row >= chunks_per_dim_)
        return false;
    if (chunk_col >= chunks_per_dim_)
        return false;
    return true;
}

SInt Grid2D::LocateVertexInChunk(
    const SInt chunk, const SInt local_row, const SInt local_col, Direction direction) const {
    SInt offset = OffsetForChunk(chunk);

    SInt chunk_row, chunk_col;
    Decode(chunk, chunk_row, chunk_col);

    SInt rows = rows_per_chunk_ + (chunk_row < remaining_rows_);
    SInt cols = cols_per_chunk_ + (chunk_col < remaining_cols_);

    SInt local_neighbor_row = 0, local_neighbor_col = 0;
    switch (direction) {
        case Right:
            local_neighbor_row = local_row;
            local_neighbor_col = 0;
            break;
        case Left:
            local_neighbor_row = local_row;
            local_neighbor_col = cols - 1;
            break;
        case Up:
            local_neighbor_row = rows - 1;
            local_neighbor_col = local_col;
            break;
        case Down:
            local_neighbor_row = 0;
            local_neighbor_col = local_col;
            break;
        default:
            break;
    }
    return offset + (local_neighbor_row * cols + local_neighbor_col);
}

void Grid2D::GenerateEdge(const SInt source, const SInt target) {
    SInt edge_seed = std::min(source, target) * total_rows_ * total_cols_ + std::max(source, target);
    SInt h         = sampling::Spooky::hash(config_.seed + edge_seed);
    if (rng_.GenerateBinomial(h, 1, edge_probability_)) {
        PushEdge(source, target);
    }
}

inline SSInt Grid2D::DirectionRow(Direction direction) const {
    switch (direction) {
        case Up:
            return -1;
        case Down:
            return 1;
        default:
            return 0;
    }
}

inline SSInt Grid2D::DirectionColumn(Direction direction) const {
    switch (direction) {
        case Left:
            return -1;
        case Right:
            return 1;
        default:
            return 0;
    }
}

SInt Grid2D::OffsetForChunk(const SInt chunk) const {
    SInt chunk_row, chunk_col;
    Decode(chunk, chunk_row, chunk_col);

    // Compute start vertex coordinates from chunk
    SInt vertex_row = chunk_row * rows_per_chunk_ + std::min(chunk_row, remaining_rows_);
    SInt vertex_col = chunk_col * cols_per_chunk_ + std::min(chunk_col, remaining_cols_);

    // Compute offset of start vertex
    SInt upper_rectangle     = vertex_row * total_cols_;
    SInt upperleft_rectangle = vertex_col * vertex_row;

    SInt next_vertex_row = (chunk_row + 1) * rows_per_chunk_ + std::min(chunk_row + 1, remaining_rows_);
    SInt left_rectangle  = vertex_col * next_vertex_row;

    return upper_rectangle + left_rectangle - upperleft_rectangle;
}

inline void Grid2D::Decode(const SInt id, SInt& x, SInt& y) const {
    x = id / chunks_per_dim_;
    y = id % chunks_per_dim_;
}

// Chunk/vertex coding
inline SInt Grid2D::Encode(const SInt x, const SInt y) const {
    return x * chunks_per_dim_ + y;
}
} // namespace kagen
