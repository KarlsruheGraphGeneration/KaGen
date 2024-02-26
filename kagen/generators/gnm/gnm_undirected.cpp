#include "kagen/generators/gnm/gnm_undirected.h"

#include "kagen/sampling/hash.hpp"

namespace kagen {
std::unique_ptr<Generator>
GNMUndirectedFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    if (config.n <= (1ull << 31)) {
        return std::make_unique<GNMUndirectedSmall>(config, rank, size);
    }
    return std::make_unique<GNMUndirectedBig>(config, rank, size);
}

template <typename BigInt>
GNMUndirected<BigInt>::GNMUndirected(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size),
      rng_(config) {}

template <typename BigInt>
void GNMUndirected<BigInt>::GenerateEdgeList() {
    leftover_chunks_ = config_.k % size_;
    SInt num_chunks  = config_.k / size_ + ((SInt)rank_ < leftover_chunks_);
    SInt row         = rank_ * num_chunks + ((SInt)rank_ >= leftover_chunks_ ? leftover_chunks_ : 0);

    nodes_per_chunk_ = config_.n / config_.k;
    remaining_nodes_ = config_.n % config_.k;

    SInt start_chunk = rank_ * (config_.k / size_) + std::min(leftover_chunks_, (SInt)rank_);
    SInt end_chunk   = start_chunk + num_chunks;

    start_node_ = start_chunk * nodes_per_chunk_ + std::min(remaining_nodes_, start_chunk);
    end_node_   = end_chunk * nodes_per_chunk_ + std::min(remaining_nodes_, end_chunk);
    num_nodes_  = end_node_ - start_node_;

    for (SInt i = 0; i < num_chunks; i++) {
        GenerateChunks(row);
        row++;
    }

    SetVertexRange(start_node_, start_node_ + num_nodes_);
}

template <typename BigInt>
void GNMUndirected<BigInt>::GenerateChunks(const SInt row) {
    QueryTriangular(config_.m, config_.k, config_.k, row, row, 0, 0, 1);
}

template <typename BigInt>
void GNMUndirected<BigInt>::QueryTriangular(
    const SInt m, const SInt num_rows, const SInt num_columns, const SInt row_id, const SInt column_id,
    const SInt offset_row, const SInt offset_column, const SInt level) {
    // Stop if there are no edges left
    if (m <= 0)
        return;

    // Total number of edges;
    SInt    n_row       = NodesInRows(num_rows, offset_row);
    SInt    n_column    = NodesInColumns(num_columns, offset_column);
    HPFloat total_edges = NumTriangleEdges(n_row, n_column, config_.self_loops);

    // Base Case if only one chunk is left
    if (num_rows == 1 && num_columns == 1) {
        GenerateTriangularEdges(m, offset_row, offset_column);
        return;
    }
    // Find splitter
    SInt row_splitter    = (num_rows + 1) / 2;
    SInt column_splitter = (num_columns + 1) / 2;

    // Compute nodes/edges per quadrant
    SInt    ul_nodes_row    = NodesInRows(row_splitter, offset_row);
    SInt    ll_nodes_row    = NodesInRows(num_rows / 2, offset_row + row_splitter);
    SInt    ll_nodes_column = NodesInColumns(column_splitter, offset_column);
    HPFloat ul_edges        = NumTriangleEdges(ul_nodes_row, ul_nodes_row);
    HPFloat ll_edges        = NumRectangleEdges(ll_nodes_row, ll_nodes_column);
    HPFloat lr_edges        = NumTriangleEdges(ll_nodes_row, ll_nodes_row);

    // Generate variate for quadrants
    SInt chunk_start   = ChunkStart(offset_row, offset_column);
    SInt h             = sampling::Spooky::hash(config_.seed + level * config_.n + chunk_start);
    SInt upper_variate = rng_.GenerateHypergeometric(h, ul_edges, m, total_edges);
    SInt ll_variate    = rng_.GenerateHypergeometric(h, ll_edges, m - upper_variate, ll_edges + lr_edges);

    // Recursive calls for quadrants
    // Row within upper half then column will automatically be in left half
    if (row_id < offset_row + row_splitter) {
        QueryTriangular(
            upper_variate, row_splitter, column_splitter, row_id, column_id, offset_row, offset_column, level + 1);
        QueryColumnRectangle(
            ll_variate, num_rows / 2, column_splitter, row_id, offset_row + row_splitter, offset_column, level + 1);
    } else {
        QueryRowRectangle(
            ll_variate, num_rows / 2, column_splitter, row_id, offset_row + row_splitter, offset_column, level + 1);
        QueryTriangular(
            m - upper_variate - ll_variate, num_rows / 2, num_columns / 2, row_id, column_id, offset_row + row_splitter,
            offset_column + column_splitter, level + 1);
    }
}

template <typename BigInt>
void GNMUndirected<BigInt>::QueryRowRectangle(
    const SInt m, const SInt num_rows, const SInt num_columns, const SInt row_id, const SInt offset_row,
    const SInt offset_column, const SInt level) {
    // Stop if there are no edges left
    if (m <= 0)
        return;

    // Total number of edges;
    SInt    n_row       = NodesInRows(num_rows, offset_row);
    SInt    n_column    = NodesInColumns(num_columns, offset_column);
    HPFloat total_edges = NumRectangleEdges(n_row, n_column);

    // Base Case if only one chunk is left
    if (num_rows == 1 && num_columns == 1) {
        GenerateRectangleEdges(m, row_id, offset_column);
        return;
    }

    // Find splitters
    SInt row_splitter    = (num_rows + 1) / 2;
    SInt column_splitter = (num_columns + 1) / 2;

    // Compute nodes/edges per quadrant
    SInt    ul_nodes_row    = NodesInRows(row_splitter, offset_row);
    SInt    ul_nodes_column = NodesInColumns(column_splitter, offset_column);
    SInt    ur_nodes_column = NodesInColumns(num_columns / 2, offset_column + column_splitter);
    HPFloat ul_edges        = NumRectangleEdges(ul_nodes_row, ul_nodes_column);
    HPFloat ur_edges        = NumRectangleEdges(ul_nodes_row, ur_nodes_column);

    // Generate variate for upper/lower half
    SInt chunk_start   = ChunkStart(offset_row, offset_column);
    SInt h             = sampling::Spooky::hash(config_.seed + level * config_.n + chunk_start);
    SInt upper_variate = rng_.GenerateHypergeometric(h, ul_edges + ur_edges, m, total_edges);

    // Recursive calls for quadrants
    // Row in upper half
    if (row_id < offset_row + row_splitter) {
        // Generate variate for upper left quadrant
        SInt ul_variate = rng_.GenerateHypergeometric(h, ul_edges, upper_variate, ul_edges + ur_edges);

        QueryRowRectangle(ul_variate, row_splitter, column_splitter, row_id, offset_row, offset_column, level + 1);
        QueryRowRectangle(
            upper_variate - ul_variate, row_splitter, num_columns / 2, row_id, offset_row,
            offset_column + column_splitter, level + 1);
    } else { // lower half
        // Compute nodes/edges per quadrant
        SInt    ll_nodes_row = NodesInRows(num_rows / 2, offset_row + row_splitter);
        HPFloat ll_edges     = NumRectangleEdges(ll_nodes_row, ul_nodes_column);
        HPFloat lr_edges     = NumRectangleEdges(ll_nodes_row, ur_nodes_column);

        // Generate variate for lower left quadrant
        SInt ll_variate = rng_.GenerateHypergeometric(h, ll_edges, m - upper_variate, ll_edges + lr_edges);

        QueryRowRectangle(
            ll_variate, num_rows / 2, column_splitter, row_id, offset_row + row_splitter, offset_column, level + 1);
        QueryRowRectangle(
            m - upper_variate - ll_variate, num_rows / 2, num_columns / 2, row_id, offset_row + row_splitter,
            offset_column + column_splitter, level + 1);
    }
}

template <typename BigInt>
void GNMUndirected<BigInt>::QueryColumnRectangle(
    const SInt m, const SInt num_rows, const SInt num_columns, const SInt column_id, const SInt offset_row,
    const SInt offset_column, const SInt level) {
    // Stop if there are no edges left
    if (m <= 0)
        return;

    // Total number of edges;
    SInt    n_row       = NodesInRows(num_rows, offset_row);
    SInt    n_column    = NodesInColumns(num_columns, offset_column);
    HPFloat total_edges = (HPFloat)n_row * n_column;

    // Base Case if only one chunk is left
    if (num_rows == 1 && num_columns == 1) {
        if (offset_row == offset_column)
            return;
        GenerateRectangleEdges(m, offset_row, column_id);
        return;
    }

    // Find splitters
    SInt row_splitter    = (num_rows + 1) / 2;
    SInt column_splitter = (num_columns + 1) / 2;

    // Compute nodes/edges per quadrant
    SInt    ul_nodes_row    = NodesInRows(row_splitter, offset_row);
    SInt    ul_nodes_column = NodesInColumns(column_splitter, offset_column);
    SInt    ur_nodes_column = NodesInColumns(num_columns / 2, offset_column + column_splitter);
    SInt    ll_nodes_row    = NodesInRows(num_rows / 2, offset_row + row_splitter);
    HPFloat ul_edges        = NumRectangleEdges(ul_nodes_row, ul_nodes_column);
    HPFloat ur_edges        = NumRectangleEdges(ul_nodes_row, ur_nodes_column);
    HPFloat ll_edges        = NumRectangleEdges(ll_nodes_row, ul_nodes_column);
    HPFloat lr_edges        = NumRectangleEdges(ll_nodes_row, ur_nodes_column);

    // Generate variate for upper/lower half
    SInt chunk_start   = ChunkStart(offset_row, offset_column);
    SInt h             = sampling::Spooky::hash(config_.seed + level * config_.n + chunk_start);
    SInt upper_variate = rng_.GenerateHypergeometric(h, ul_edges + ur_edges, m, total_edges);
    SInt ul_variate    = rng_.GenerateHypergeometric(h, ul_edges, upper_variate, ul_edges + ur_edges);
    SInt ll_variate    = rng_.GenerateHypergeometric(h, ll_edges, m - upper_variate, ll_edges + lr_edges);

    // Recursive calls for quadrants
    // Column in left half
    if (column_id < offset_column + column_splitter) {
        QueryColumnRectangle(
            ul_variate, row_splitter, column_splitter, column_id, offset_row, offset_column, level + 1);
        QueryColumnRectangle(
            ll_variate, num_rows / 2, column_splitter, column_id, offset_row + row_splitter, offset_column, level + 1);
        // Column in right half
    } else {
        QueryColumnRectangle(
            upper_variate - ul_variate, row_splitter, num_columns / 2, column_id, offset_row,
            offset_column + column_splitter, level + 1);
        QueryColumnRectangle(
            m - upper_variate - ll_variate, num_rows / 2, num_columns / 2, column_id, offset_row + row_splitter,
            offset_column + column_splitter, level + 1);
    }
}

template <typename BigInt>
void GNMUndirected<BigInt>::GenerateTriangularEdges(const SInt m, const SInt row_id, const SInt column_id) {
    // No self loops -> skip diagonal entries
    SInt offset_row    = OffsetInRow(row_id);
    SInt offset_column = OffsetInColumn(column_id);
    if (!config_.self_loops)
        offset_row += 1;
    // bool local_row = (offset_row >= start_node_ && offset_row < end_node_);

    // Number edges
    SInt    n_row       = NodesInRow(row_id);
    SInt    n_column    = NodesInColumn(column_id);
    HPFloat total_edges = NumTriangleEdges(n_row, n_column, config_.self_loops);

    // Sample from [1, total_edges]
    SInt h = sampling::Spooky::hash(config_.seed + (((row_id + 1) * row_id) / 2) + column_id);
    rng_.GenerateSample(h, total_edges, m, [&](SInt sample) {
        // Absolute triangular point
        // if (loops) sqr = (sqrt(8*((double)sample-1)+1) - 1)/2 + 1;
        SInt sqr = sqrt(8 * (sample - 1) + 1);
        // TODO: Nasty hack
        while (sqr * sqr > 8 * (sample - 1) + 1)
            sqr--;
        SInt i = (sqr - 1) / 2;
        SInt j = (sample - 1) - i * (i + 1) / 2;

        SInt from = i + offset_row;
        SInt to   = j + offset_column;

        PushEdge(from, to);
        if (to >= start_node_ && to < end_node_) {
            PushEdge(to, from);
        }
    });
}

template <typename BigInt>
void GNMUndirected<BigInt>::GenerateRectangleEdges(const SInt m, const SInt row_id, const SInt column_id) {
    SInt                  offset_row    = OffsetInRow(row_id);
    SInt                  offset_column = OffsetInColumn(column_id);
    [[maybe_unused]] bool local_row     = (offset_row >= start_node_ && offset_row < end_node_);

    // Sample from [1, num_edges]
    SInt    n_row       = NodesInRow(row_id);
    SInt    n_column    = NodesInColumn(column_id);
    HPFloat total_edges = NumRectangleEdges(n_row, n_column);

    // Sample from [1, total_edges]
    SInt h = sampling::Spooky::hash(config_.seed + (((row_id + 1) * row_id) / 2) + column_id);
    rng_.GenerateSample(h, total_edges, m, [&](SInt sample) {
        SInt i = (sample - 1) / n_column;
        SInt j = (sample - 1) % n_column;

        SInt from = i + offset_row;
        SInt to   = j + offset_column;

        if (local_row) {
            PushEdge(from, to);
        }
        if (to >= start_node_ && to < end_node_) {
            PushEdge(to, from);
        }
    });
}

template <typename BigInt>
inline SInt GNMUndirected<BigInt>::NodesInRows(const SInt rows, const SInt offset) const {
    return nodes_per_chunk_ * rows + std::min(remaining_nodes_ - offset, rows);
}

template <typename BigInt>
inline SInt GNMUndirected<BigInt>::NodesInColumns(const SInt columns, const SInt offset) const {
    return nodes_per_chunk_ * columns + std::min(remaining_nodes_ - offset, columns);
}

template <typename BigInt>
inline SInt GNMUndirected<BigInt>::NodesInRow(const SInt row) const {
    return nodes_per_chunk_ + (row < leftover_chunks_);
}

template <typename BigInt>
inline SInt GNMUndirected<BigInt>::NodesInColumn(const SInt column) const {
    return nodes_per_chunk_ + (column < leftover_chunks_);
}

template <typename BigInt>
inline SInt GNMUndirected<BigInt>::OffsetInRow(const SInt row) const {
    return nodes_per_chunk_ * row + std::min(remaining_nodes_, row);
}

template <typename BigInt>
inline SInt GNMUndirected<BigInt>::OffsetInColumn(const SInt column) const {
    return nodes_per_chunk_ * column + std::min(remaining_nodes_, column);
}

template <typename BigInt>
inline SInt GNMUndirected<BigInt>::ChunkStart(const SInt row, const SInt column) const {
    return (((row + 1) * row) / 2) + column;
}

template <typename BigInt>
inline HPFloat GNMUndirected<BigInt>::NumTriangleEdges(const HPFloat row, const HPFloat column, bool loops) const {
    return (loops && config_.self_loops) ? row * (column + 1) / 2 : row * (column - 1) / 2;
}

template <typename BigInt>
inline HPFloat GNMUndirected<BigInt>::NumRectangleEdges(const HPFloat row, const HPFloat column) const {
    return row * column;
}
} // namespace kagen
