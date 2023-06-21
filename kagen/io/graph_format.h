#pragma once

#include <memory>
#include <stdexcept>
#include <string>

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/kagen.h"

namespace kagen {
enum ReaderDeficits {
    NONE                    = 0,
    REQUIRES_REDISTRIBUTION = 1,
    EDGE_LIST_ONLY          = 2,
    CSR_ONLY                = 4,
    UNKNOWN_NUM_VERTICES    = 8,
    UNKNOWN_NUM_EDGES       = 16,
};

class GraphReader {
public:
    using GraphSize = std::pair<SInt, SInt>;

    virtual ~GraphReader() = default;

    virtual GraphSize ReadSize() = 0;

    virtual Graph Read(SInt from_vertex, SInt to_vertex, SInt to_edge, GraphRepresentation representation) = 0;

    virtual SInt FindNodeByEdge(SInt edge) = 0;

    virtual int Deficits() const {
        return ReaderDeficits::NONE;
    }
};

struct GraphInfo {
    GraphInfo() = default;

    GraphInfo(const Graph& graph, MPI_Comm comm);

    SInt global_n           = 0;
    SInt global_m           = 0;
    bool has_vertex_weights = false;
    bool has_edge_weights   = false;
};

class GraphWriter {
public:
    GraphWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size);

    virtual ~GraphWriter() = default;

    // Note: may not  assume that previous passes were called on the same GraphWriter object
    virtual bool Write(int pass, const std::string& filename) = 0;

protected:
    void RequiresCoordinates() const;
    void Requires2DCoordinates() const;
    void Requires3DCoordinates() const;
    void IgnoresVertexWeights() const;
    void IgnoresEdgeWeights() const;

    const OutputGraphConfig& config_;

    GraphInfo info_;
    Graph&    graph_;
    PEID      rank_;
    PEID      size_;
};

class StandardGraphWriter : public GraphWriter {
public:
    StandardGraphWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size);

    bool Write(int pass, const std::string& filename) final;

protected:
    virtual void WriteHeader(const std::string& filename, SInt n, SInt m) = 0;

    virtual bool WriteBody(const std::string& filename) = 0;

    virtual void WriteFooter(const std::string& filename);
};

class FileFormatFactory {
public:
    virtual ~FileFormatFactory() = default;

    virtual std::string DefaultExtension() const = 0;

    virtual std::unique_ptr<GraphReader> CreateReader(
        [[maybe_unused]] const InputGraphConfig& config, [[maybe_unused]] PEID rank, [[maybe_unused]] PEID size) const {
        return nullptr;
    }

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
    std::string DefaultExtension() const final {
        return "";
    }

    std::unique_ptr<GraphWriter> CreateWriter(
        const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank,
        const PEID size) const final {
        return std::make_unique<NoopWriter>(config, graph, info, rank, size);
    }
};
} // namespace kagen
