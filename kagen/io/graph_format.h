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

class GraphReader {
public:
    using GraphSize = std::pair<SInt, SInt>;

    virtual ~GraphReader() = default;

    virtual GraphSize ReadSize() = 0;

    virtual Graph Read(SInt from, SInt to, GraphRepresentation representation) = 0;

    virtual SInt FindNodeByEdge(SInt edge) = 0;
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
    virtual std::string DefaultExtension() const = 0;

    virtual std::unique_ptr<GraphReader> CreateReader(const InputGraphConfig& config) const = 0;

    virtual std::unique_ptr<GraphWriter>
    CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const = 0;
};
} // namespace kagen
