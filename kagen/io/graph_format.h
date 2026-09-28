#pragma once

#include "kagen/context.h"
#include "kagen/kagen.h"

#include <mpi.h>

#include <memory>
#include <stdexcept>
#include <string>

namespace kagen {
enum ReaderDeficits {
    //! Indicates that the reader has no deficits, i.e., none of the deficits enumerated is applicable.
    NONE = 0,

    //! Indicates that the reader does not respect the given vertex- or edge range, and that outgoing edges of one
    //! vertex might exist on multiple PEs
    REQUIRES_REDISTRIBUTION = 1,

    //! Indicates that the reader always returns the graph as an edge list.
    EDGE_LIST_ONLY = 2,

    //! Indicates that the reader always returns the graph as CSR.
    CSR_ONLY = 4,

    //! Indicates that the `ReadSize()` function does not return the correct number of vertices.
    UNKNOWN_NUM_VERTICES = 8,

    //! Indicates that the `ReadSize()` function does not return the correct number of edges.
    UNKNOWN_NUM_EDGES = 16,
};

class GraphReader {
public:
    //! Number of vertices (first) and number of edges (second) in the graph.
    using GraphSize = std::pair<SInt, SInt>;

    virtual ~GraphReader() = default;

    /*!
     * Should return the number of vertices and the number of edges in the graph. If the size of the graph can not be
     * determined without reading the graph, the reader should report the `ReaderDeficits::UNKNOWN_NUM_VERTICES` and/or
     * `ReaderDeficits::UNKNOWN_NUM_EDGES` deficits via `Deficits()`.
     * The returned values are used to determine the parameters passed to `Read()`.
     *
     * @return Number of vertices and number of edges in the graph, or any other values appropriate to guide the vertex
     * or edge range passed to `Read()`.
     */
    virtual GraphSize ReadSize() = 0;

    /*!
     * Reads the part of the graph specified by the given vertex- or edge range. If the reader cannot respect the range,
     * it should report the `ReaderDeficits::REQUIRES_REDISTRIBUTION` deficit.
     *
     * @param from_vertex The first vertex ID to be read.
     * @param to_vertex The last vertex ID to be read; this might possibly be larger than the number of vertices in the
     * graph.
     * @param to_edge The last edge ID to be read; this might possibly be larger than the number of edges in the graph.
     * @param representation Whether the graph should be returned as an edge list or in CSR. This parameter can be
     * ignored if indicated via the appropriate deficit.
     * @return The read graph.
     */
    virtual Graph Read(SInt from_vertex, SInt to_vertex, SInt to_edge, GraphRepresentation representation) = 0;

    /*!
     * Finds the tail vertex of the given edge ID.
     *
     * @param edge The edge ID.
     * @return The tail vertex of the given edge ID.
     */
    virtual SInt FindNodeByEdge(SInt edge) = 0;

    /*!
     * Whether this reader can read exactly the contiguous global edge range `[from_edge, to_edge)` directly
     * from disk via `ReadStrictEdgeRange()`, starting and/or ending mid-adjacency if a strict edge boundary
     * falls inside a vertex. This enables reading a strictly edge-balanced (`GraphDistribution::BALANCE_EDGES_TRUE`)
     * partition without the vertex-balanced-read + redistribution postprocessing pass. Requires the number of
     * edges reported by `ReadSize()` to be exact (i.e., the `UNKNOWN_NUM_EDGES` deficit must not be set).
     */
    virtual bool CanReadStrictEdgeRange() const {
        return false;
    }

    /*!
     * Whether a `ReadStrictEdgeRange()` result may only be trusted after collectively verifying that the input
     * is globally tail-sorted. This is `true` for edge-list readers, whose on-disk edge order is arbitrary
     * (a contiguous edge range only maps to a coherent vertex partition if the edges are sorted by tail), and
     * `false` for readers whose on-disk order is tail-grouped by construction (ParHIP, METIS). If the sortedness
     * check fails, the caller falls back to the redistribution path.
     */
    virtual bool StrictEdgeRangeRequiresSortednessCheck() const {
        return false;
    }

    /*!
     * Reads exactly the contiguous global edge range `[from_edge, to_edge)` directly from disk. The returned
     * graph holds this PE's local edges (edge list) or the CSR slice spanning the physically-present tail
     * vertices; the first and/or last vertex may be partial (its adjacency shared with an adjacent PE). The
     * returned `vertex_range` is provisional (the physically-present tail vertices `[first_tail, last_tail + 1)`);
     * the caller replaces it with the gap-free partition and fills in the split-vertex metadata. Only called
     * when `CanReadStrictEdgeRange()` returns `true`.
     *
     * @param from_edge The first global edge ID to read (inclusive).
     * @param to_edge The last global edge ID to read (exclusive).
     * @param representation Whether to return the slice as an edge list or in CSR.
     * @return The read graph slice.
     */
    virtual Graph ReadStrictEdgeRange(
        [[maybe_unused]] SInt from_edge, [[maybe_unused]] SInt to_edge,
        [[maybe_unused]] GraphRepresentation representation) {
        throw std::runtime_error("this reader does not support direct strict edge-range reads");
    }

    virtual int Deficits() const {
        return ReaderDeficits::NONE;
    }

    bool HasDeficit(const ReaderDeficits deficit) const {
        return Deficits() & deficit;
    }
};

struct GraphInfo {
    GraphInfo() = default;
    GraphInfo(const Graph& graph, MPI_Comm comm);
    GraphInfo(const GraphInfo& local, MPI_Comm comm);

    SInt local_n            = 0;
    SInt local_m            = 0;
    SInt global_n           = 0;
    SInt global_m           = 0;
    SInt offset_n           = 0;
    SInt offset_m           = 0;
    bool has_vertex_weights = false;
    bool has_edge_weights   = false;
};

class GraphWriter {
public:
    GraphWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size);

    virtual ~GraphWriter() = default;

    // Note: may not  assume that previous passes were called on the same GraphWriter object
    /*!
     * Write the graph to the given file. This can happen interleaved with the output from other PEs by calling this
     * function repeatedly in multiple passes, as indicated by the `pass` parameter. Once this function returns `false`,
     * no more passes will be made.
     *
     * More precisely, output works as follows: in sequential mode, the procedure first calls Write(pass = 0) on PE 0,
     * then on PE 1, ..., on PE <nproc>. If the functions returns `true`, the process is repeated with `pass = 1`, and
     * so on, until it returns `false` (the return value is expected to be the same on all PEs during a pass).
     *
     * The writer *MAY NOT* assume that all passes are called on the same writer object.
     *
     * @param pass The current pass.
     * @param filename The name of the file to write to.
     * @return Whether IO requires another pass after this one.
     */
    virtual bool Write(int pass, const std::string& filename) = 0;

protected:
    //! Indicates that the writer requires 2D or 3D vertex coordinates (to be called from `Write()`).
    void RequiresCoordinates() const;

    //! Indicates that the writer requires 2D vertex coordinates (to be called from `Write()`).
    void Requires2DCoordinates() const;

    //! Indicates that the writer requires 3D vertex coordinates (to be called from `Write()`).
    void Requires3DCoordinates() const;

    //! Indicates that the writer will ignore vertex coordinates if they are present (to be called from `Write()`).
    void IgnoresVertexWeights() const;

    //! Indicates that the writer will ignore edge coordinates if they are present (to be called from `Write()`).
    void IgnoresEdgeWeights() const;

    const OutputGraphConfig& config_;

    GraphInfo info_;
    Graph&    graph_;
    PEID      rank_;
    PEID      size_;
};

/*!
 * Default graph writer which performs three passes: writing a file header (only on the root PE in sequential output
 * mode), the contents and a footer (only on the root PE in sequential mode).
 */
class StandardGraphWriter : public GraphWriter {
public:
    StandardGraphWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size);

    bool Write(int pass, const std::string& filename) final;

protected:
    virtual void WriteHeader(const std::string& filename, SInt n, SInt m) = 0;

    virtual bool WriteBody(const std::string& filename) = 0;

    virtual void WriteFooter(const std::string& filename);
};

/*!
 * Factory responsible for creating reader and writer objects for a file format.
 */
class FileFormatFactory {
public:
    virtual ~FileFormatFactory() = default;

    /*!
     * Unique file extension linked to this file format. This is used when determining the file format from a file
     * extension and vice-versa.
     *
     * @return The default file extension for this file format.
     */
    virtual std::vector<std::string> DefaultExtensions() const = 0;

    //! May return `nullptr` if reading the file format is not supported.
    virtual std::unique_ptr<GraphReader> CreateReader(
        [[maybe_unused]] const InputGraphConfig& config, [[maybe_unused]] PEID rank, [[maybe_unused]] PEID size) const {
        return nullptr;
    }

    //! May return `nullptr` if writing the file format is not supported.
    virtual std::unique_ptr<GraphWriter> CreateWriter(
        [[maybe_unused]] const OutputGraphConfig& config, [[maybe_unused]] Graph& graph,
        [[maybe_unused]] GraphInfo info, [[maybe_unused]] PEID rank, [[maybe_unused]] PEID size) const {
        return nullptr;
    }
};

//
// Noop format = empty output, no input
//
class NoopWriter : public GraphWriter {
public:
    NoopWriter(const OutputGraphConfig& config, Graph& graph, const GraphInfo info, PEID rank, PEID size)
        : GraphWriter(config, graph, info, rank, size) {}

    inline bool Write(int, const std::string&) final {
        return false;
    }
};

class NoopFactory : public FileFormatFactory {
public:
    std::vector<std::string> DefaultExtensions() const final {
        return {};
    }

    std::unique_ptr<GraphWriter> CreateWriter(
        const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank,
        const PEID size) const final {
        return std::make_unique<NoopWriter>(config, graph, info, rank, size);
    }
};
} // namespace kagen
