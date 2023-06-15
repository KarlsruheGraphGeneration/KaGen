#pragma once

#include <memory>
#include <stdexcept>
#include <string>

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/kagen.h"

namespace kagen {
class IOError : public std::exception {
public:
    IOError(std::string what) : _what(std::move(what)) {}

    const char* what() const noexcept override {
        return _what.c_str();
    }

private:
    std::string _what;
};

enum ReaderDeficits {
    NONE                     = 0,
    REQUIRES_REDISTRIBUTION  = 1,
    EDGE_LIST_ONLY           = 2,
    CSR_ONLY                 = 4,
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

class GraphWriter {
public:
    GraphWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm);

    virtual ~GraphWriter() = default;

    virtual void Write(bool report_progress) = 0;

protected:
    bool HasVertexWeights() const;
    bool HasEdgeWeights() const;

    const OutputGraphConfig& config_;

    EdgeList&      edges_;
    VertexRange&   vertex_range_;
    Coordinates&   coordinates_;
    VertexWeights& vertex_weights_;
    EdgeWeights&   edge_weights_;

    MPI_Comm comm_;
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
        [[maybe_unused]] MPI_Comm comm) const {
        return nullptr;
    }
};

//
// Noop format = empty output, no input
//
class NoopWriter : public GraphWriter {
public:
    NoopWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) : GraphWriter(config, graph, comm) {}

    void Write(bool) final {}
};

class NoopFactory : public FileFormatFactory {
public:
    std::string DefaultExtension() const final {
        return "";
    }

    std::unique_ptr<GraphWriter>
    CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const final {
        return std::make_unique<NoopWriter>(config, graph, comm);
    }
};
} // namespace kagen
