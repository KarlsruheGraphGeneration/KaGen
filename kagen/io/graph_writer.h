#pragma once

#include <string>

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/definitions.h"

namespace kagen {
class GraphWriter {
public:
    GraphWriter(EdgeList& edges, VertexRange vertex_range, Coordinates& coordinates, MPI_Comm comm);
    virtual ~GraphWriter();

    virtual std::string DefaultExtension() const = 0;

    virtual void Write(const PGeneratorConfig& config) = 0;

protected:
    EdgeList&    edges_;
    VertexRange  vertex_range_;
    Coordinates& coordinates_;
    MPI_Comm     comm_;
};

class SequentialGraphWriter : public GraphWriter {
protected:
    enum Requirement {
        NONE           = 0,
        SORTED_EDGES   = 1 << 1,
        COORDINATES    = 1 << 2,
        COORDINATES_2D = 1 << 4,
        COORDINATES_3D = 1 << 8,
    };

public:
    SequentialGraphWriter(EdgeList& edges, VertexRange vertex_range, Coordinates& coordinates, MPI_Comm comm);

    virtual void Write(const PGeneratorConfig& config) override;

protected:
    virtual void AppendHeaderTo(const std::string& filename, SInt n, SInt m) = 0;

    virtual void AppendTo(const std::string& filename) = 0;

    virtual void AppendFooterTo(const std::string& filename);

    virtual Requirement Requirements() const;

private:
    static void CreateFile(const std::string& filename);
};

class NoopWriter : public GraphWriter {
public:
    NoopWriter(EdgeList& edges, VertexRange vertex_range, Coordinates& coordinates, MPI_Comm comm);

    std::string DefaultExtension() const final;

    void Write(const PGeneratorConfig& config) final;
};
} // namespace kagen
