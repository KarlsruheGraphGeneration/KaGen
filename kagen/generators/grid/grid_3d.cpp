#include "kagen/generators/grid/grid_3d.h"

#include "kagen/sampling/hash.hpp"

namespace kagen {
std::unique_ptr<Generator>
Grid3DFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<Grid3D>(config, rank, size);
}

PGeneratorConfig
Grid3DFactory::NormalizeParameters(PGeneratorConfig config, PEID, const PEID size, const bool output) const {
    EnsureCubicPowerOfTwoChunkSize(config, size, output);

    if (config.p == 0) {
        if (config.grid_x == 0 || config.grid_y == 0 || config.grid_z == 0 || config.m == 0) {
            throw ConfigurationError("at least two parameters out of {(x, y, z), m, p} must be nonzero");
        }

        const SInt num_deg4_vertices = config.periodic ? 0 : 4 * config.grid_x + 4 * config.grid_y + 4 * config.grid_z;
        const SInt num_deg5_vertices = config.periodic ? 0
                                                       : 2 * (config.grid_x - 1) * (config.grid_y - 1)
                                                             + 2 * (config.grid_y - 1) * (config.grid_z - 1)
                                                             + 2 * (config.grid_x - 1) * (config.grid_z - 1);
        const SInt num_deg6_vertices =
            config.grid_x * config.grid_y * config.grid_z - num_deg4_vertices - num_deg5_vertices;

        config.p = 2.0 * config.m / (6 * num_deg6_vertices + 5 * num_deg5_vertices + 4 * num_deg4_vertices);
        if (output) {
            std::cout << "Setting edge probability to " << config.p << std::endl;
            if (config.p > 1) {
                std::cerr << "Warning: configuration infeasible, too many edges\n";
            }
        }
    }

    return config;
}

Grid3D::Grid3D(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rng_(config),
      rank_(rank),
      size_(size) {}

void Grid3D::GenerateEdgeList() {
    // Init dimensions
    // TODO: Only tested for cube PEs and one chunk per PE
    total_x_          = config_.grid_x;
    total_y_          = config_.grid_y;
    total_z_          = config_.grid_z;
    edge_probability_ = config_.p;

    // Init chunks
    total_chunks_   = config_.k;
    chunks_per_dim_ = std::cbrt(total_chunks_);

    SInt leftover_chunks = total_chunks_ % size_;
    SInt num_chunks      = (total_chunks_ / size_) + ((SInt)rank_ < leftover_chunks);
    SInt start_chunk     = rank_ * num_chunks + ((SInt)rank_ >= leftover_chunks ? leftover_chunks : 0);
    SInt end_chunk       = start_chunk + num_chunks;

    // Chunk distribution
    y_per_chunk_ = total_y_ / chunks_per_dim_;
    remaining_y_ = total_y_ % chunks_per_dim_;

    x_per_chunk_ = total_x_ / chunks_per_dim_;
    remaining_x_ = total_x_ % chunks_per_dim_;

    z_per_chunk_ = total_z_ / chunks_per_dim_;
    remaining_z_ = total_z_ % chunks_per_dim_;

    vertices_per_chunk_ = x_per_chunk_ * y_per_chunk_ * z_per_chunk_;

    start_node_ = OffsetForChunk(start_chunk);
    end_node_   = OffsetForChunk(end_chunk);
    num_nodes_  = end_node_ - start_node_;

    for (SInt i = 0; i < num_chunks; i++) {
        GenerateChunk(start_chunk + i);
    }

    if (config_.coordinates) {
        for (SInt z = 0; z < total_z_; ++z) {
            for (SInt y = 0; y < total_y_; ++y) {
                for (SInt x = 0; x < total_x_; ++x) {
                    PushCoordinate(1.0 * x / total_x_, 1.0 * y / total_y_, 1.0 * z / total_z_);
                }
            }
        }
    }

    SetVertexRange(start_node_, start_node_ + num_nodes_);
}

void Grid3D::GenerateChunk(const SInt chunk) {
    SInt start_vertex = OffsetForChunk(chunk);
    SInt end_vertex   = OffsetForChunk(chunk + 1);
    for (SInt i = start_vertex; i < end_vertex; i++) {
        GenerateEdges(chunk, i);
    }
}

void Grid3D::GenerateEdges(const SInt chunk, const SInt vertex) {
    QueryInDirection(chunk, vertex, Direction::Right);
    QueryInDirection(chunk, vertex, Direction::Left);
    QueryInDirection(chunk, vertex, Direction::Up);
    QueryInDirection(chunk, vertex, Direction::Down);
    QueryInDirection(chunk, vertex, Direction::Front);
    QueryInDirection(chunk, vertex, Direction::Back);
}

void Grid3D::QueryInDirection(const SInt chunk, const SInt vertex, Direction direction) {
    SInt offset       = OffsetForChunk(chunk);
    SInt local_vertex = vertex - offset;

    SInt chunk_x, chunk_y, chunk_z;
    Decode(chunk, chunk_x, chunk_y, chunk_z);

    SInt xs = x_per_chunk_ + (chunk_x < remaining_x_);
    SInt ys = y_per_chunk_ + (chunk_y < remaining_y_);
    SInt zs = z_per_chunk_ + (chunk_z < remaining_z_);

    SInt local_x = local_vertex % xs;
    SInt local_y = (local_vertex / xs) % ys;
    SInt local_z = local_vertex / (xs * ys);

    SSInt local_neighbor_x = (SSInt)local_x + DirectionX(direction);
    SSInt local_neighbor_y = (SSInt)local_y + DirectionY(direction);
    SSInt local_neighbor_z = (SSInt)local_z + DirectionZ(direction);

    if (IsLocalVertex(local_neighbor_x, local_neighbor_y, local_neighbor_z, xs, ys, zs)) {
        SInt neighbor_vertex = offset + (local_neighbor_x + local_neighbor_y * xs + local_neighbor_z * (xs * ys));
        GenerateEdge(vertex, neighbor_vertex);
    } else {
        // Determine neighboring chunk
        SSInt neighbor_chunk_x = (SSInt)chunk_x + DirectionX(direction);
        SSInt neighbor_chunk_y = (SSInt)chunk_y + DirectionY(direction);
        SSInt neighbor_chunk_z = (SSInt)chunk_z + DirectionZ(direction);
        if (config_.periodic) {
            neighbor_chunk_x = (neighbor_chunk_x + chunks_per_dim_) % chunks_per_dim_;
            neighbor_chunk_y = (neighbor_chunk_y + chunks_per_dim_) % chunks_per_dim_;
            neighbor_chunk_z = (neighbor_chunk_z + chunks_per_dim_) % chunks_per_dim_;
        }
        if (!IsValidChunk(neighbor_chunk_x, neighbor_chunk_y, neighbor_chunk_z))
            return;

        SInt neighbor_chunk  = Encode(neighbor_chunk_x, neighbor_chunk_y, neighbor_chunk_z);
        SInt neighbor_vertex = LocateVertexInChunk(neighbor_chunk, local_x, local_y, local_z, direction);
        GenerateEdge(vertex, neighbor_vertex);
    }
}

bool Grid3D::IsLocalVertex(
    const SSInt local_x, const SSInt local_y, const SSInt local_z, const SInt xs, const SInt ys, const SInt zs) const {
    if (local_x < 0 || local_x >= static_cast<SSInt>(xs))
        return false;
    if (local_y < 0 || local_y >= static_cast<SSInt>(ys))
        return false;
    if (local_z < 0 || local_z >= static_cast<SSInt>(zs))
        return false;
    return true;
}

bool Grid3D::IsValidChunk(const SSInt chunk_x, const SSInt chunk_y, const SSInt chunk_z) const {
    if (chunk_x < 0 || chunk_x >= static_cast<SSInt>(chunks_per_dim_))
        return false;
    if (chunk_y < 0 || chunk_y >= static_cast<SSInt>(chunks_per_dim_))
        return false;
    if (chunk_z < 0 || chunk_z >= static_cast<SSInt>(chunks_per_dim_))
        return false;
    return true;
}

SInt Grid3D::LocateVertexInChunk(
    const SInt chunk, const SInt local_x, const SInt local_y, const SInt local_z, Direction direction) const {
    SInt offset = OffsetForChunk(chunk);

    SInt chunk_x, chunk_y, chunk_z;
    Decode(chunk, chunk_x, chunk_y, chunk_z);

    SInt xs = x_per_chunk_ + (chunk_x < remaining_x_);
    SInt ys = y_per_chunk_ + (chunk_y < remaining_y_);
    SInt zs = z_per_chunk_ + (chunk_z < remaining_z_);

    SInt local_neighbor_x = 0, local_neighbor_y = 0, local_neighbor_z = 0;
    switch (direction) {
        case Right:
            local_neighbor_x = 0;
            local_neighbor_y = local_y;
            local_neighbor_z = local_z;
            break;
        case Left:
            local_neighbor_x = xs - 1;
            local_neighbor_y = local_y;
            local_neighbor_z = local_z;
            break;
        case Up:
            local_neighbor_x = local_x;
            local_neighbor_y = ys - 1;
            local_neighbor_z = local_z;
            break;
        case Down:
            local_neighbor_x = local_x;
            local_neighbor_y = 0;
            local_neighbor_z = local_z;
            break;
        case Front:
            local_neighbor_x = local_x;
            local_neighbor_y = local_y;
            local_neighbor_z = zs - 1;
            break;
        case Back:
            local_neighbor_x = local_x;
            local_neighbor_y = local_y;
            local_neighbor_z = 0;
            break;
        default:
            break;
    }
    return offset + (local_neighbor_x + local_neighbor_y * xs + local_neighbor_z * (xs * ys));
}

void Grid3D::GenerateEdge(const SInt source, const SInt target) {
    SInt edge_seed = std::min(source, target) * total_y_ * total_x_ * total_z_ + std::max(source, target);
    SInt h         = sampling::Spooky::hash(config_.seed + edge_seed);
    if (rng_.GenerateBinomial(h, 1, edge_probability_)) {
        PushEdge(source, target);
    }
}

inline SSInt Grid3D::DirectionX(Direction direction) const {
    switch (direction) {
        case Left:
            return -1;
        case Right:
            return 1;
        default:
            return 0;
    }
}

inline SSInt Grid3D::DirectionY(Direction direction) const {
    switch (direction) {
        case Up:
            return -1;
        case Down:
            return 1;
        default:
            return 0;
    }
}

inline SSInt Grid3D::DirectionZ(Direction direction) const {
    switch (direction) {
        case Front:
            return -1;
        case Back:
            return 1;
        default:
            return 0;
    }
}

SInt Grid3D::OffsetForChunk(const SInt chunk) const {
    SInt chunk_x, chunk_y, chunk_z;
    Decode(chunk, chunk_x, chunk_y, chunk_z);

    // Compute start vertex coordinates from chunk
    SInt vertex_x = chunk_x * x_per_chunk_ + std::min(chunk_x, remaining_x_);
    SInt vertex_y = chunk_y * y_per_chunk_ + std::min(chunk_y, remaining_y_);
    SInt vertex_z = chunk_z * z_per_chunk_ + std::min(chunk_z, remaining_z_);

    SInt next_vertex_y = (chunk_y + 1) * y_per_chunk_ + std::min(chunk_y + 1, remaining_y_);
    SInt next_vertex_z = (chunk_z + 1) * z_per_chunk_ + std::min(chunk_z + 1, remaining_z_);

    // Compute offset of start vertex
    SInt upper_cube        = total_x_ * vertex_y * next_vertex_z;
    SInt frontal_cube      = total_x_ * total_y_ * vertex_z;
    SInt frontal_left_cube = vertex_x * next_vertex_y * next_vertex_z;

    SInt intersect_upper_frontal        = total_x_ * vertex_y * vertex_z;
    SInt intersect_upper_frontal_left   = vertex_x * vertex_y * next_vertex_z;
    SInt intersect_frontal_frontal_left = vertex_x * next_vertex_y * vertex_z;
    SInt intersect_all                  = vertex_x * vertex_y * vertex_z;

    return upper_cube + frontal_cube + frontal_left_cube
           - (intersect_upper_frontal + intersect_upper_frontal_left + intersect_frontal_frontal_left) + intersect_all;
}

inline void Grid3D::Decode(const SInt id, SInt& x, SInt& y, SInt& z) const {
    x = id % chunks_per_dim_;
    y = (id / chunks_per_dim_) % chunks_per_dim_;
    z = id / (chunks_per_dim_ * chunks_per_dim_);
}

// Chunk/vertex coding
inline SInt Grid3D::Encode(const SInt x, const SInt y, const SInt z) const {
    return x + y * chunks_per_dim_ + z * (chunks_per_dim_ * chunks_per_dim_);
}
} // namespace kagen
