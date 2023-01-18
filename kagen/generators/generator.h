#pragma once

#include <cstdint>
#include <exception>
#include <memory>

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/definitions.h"

namespace kagen {
class Generator {
public:
    virtual ~Generator();

    void Generate(GraphRepresentation representation);

    void Finalize(MPI_Comm comm);

    SInt GetNumberOfEdges() const;

    Graph Take();

protected:
    virtual void GenerateEdgeList() = 0;

    virtual void GenerateCSR() = 0;

    virtual void FinalizeEdgeList(MPI_Comm comm);

    virtual void FinalizeCSR(MPI_Comm comm);

    void SetVertexRange(VertexRange vetrex_range);

    inline void PushCoordinate(const HPFloat x, const HPFloat y) {
        coordinates_.first.emplace_back(x, y);
    }

    inline void PushCoordinate(const HPFloat x, const HPFloat y, const HPFloat z) {
        coordinates_.second.emplace_back(x, y, z);
    }

    inline void PushVertexWeight(const SSInt weight) {
        vertex_weights_.push_back(weight);
    }

    inline void PushEdge(const SInt from, const SInt to) {
        edges_.emplace_back(from, to);
    }

    inline void PushEdgeWeight(const SSInt weight) {
        edge_weights_.push_back(weight);
    }

    inline void SetVertexRange(const SInt first_vertex, const SInt first_invalid_vertex) {
        vertex_range_ = std::make_pair(first_vertex, first_invalid_vertex);
    }

    void FilterDuplicateEdges();

    VertexRange   vertex_range_;
    EdgeList      edges_;
    XadjArray     xadj_;
    AdjncyArray   adjncy_;
    Coordinates   coordinates_;
    VertexWeights vertex_weights_;
    EdgeWeights   edge_weights_;

private:
    void Reset();

    GraphRepresentation representation_;
};

class ConfigurationError : public std::exception {
public:
    ConfigurationError(std::string what) : _what(std::move(what)) {}

    const char* what() const noexcept override {
        return _what.c_str();
    }

private:
    std::string _what;
};

class EdgeListOnlyGenerator : virtual Generator {
public:
    void GenerateCSR() final;
    void FinalizeCSR(MPI_Comm comm) final;
};

class CSROnlyGenerator : virtual Generator {
public:
    void GenerateEdgeList() final;
    void FinalizeEdgeList(MPI_Comm comm) final;
};

class GeneratorFactory {
public:
    virtual ~GeneratorFactory();

    virtual PGeneratorConfig NormalizeParameters(PGeneratorConfig config, PEID rank, PEID size, bool output) const;

    virtual std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const = 0;

protected:
    void EnsurePowerOfTwoCommunicatorSize(PGeneratorConfig& config, PEID size) const;
    void EnsureSquareChunkSize(PGeneratorConfig& config, PEID size) const;
    void EnsureCubicChunkSize(PGeneratorConfig& config, PEID size) const;
    void EnsureOneChunkPerPE(PGeneratorConfig& config, PEID size) const;
};
} // namespace kagen
