#pragma once

#include "kagen/context.h"
#include "kagen/kagen.h"

#include <mpi.h>

#include <exception>
#include <memory>

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
        graph_.coordinates.first.emplace_back(x, y);
    }

    inline void PushCoordinate(const HPFloat x, const HPFloat y, const HPFloat z) {
        graph_.coordinates.second.emplace_back(x, y, z);
    }

    inline void PushVertexWeight(const SSInt weight) {
        graph_.vertex_weights.push_back(weight);
    }

    inline void PushEdge(const SInt from, const SInt to) {
        graph_.edges.emplace_back(from, to);
    }

    inline void PushEdgeWeight(const SSInt weight) {
        graph_.edge_weights.push_back(weight);
    }

    inline void SetVertexRange(const SInt first_vertex, const SInt first_invalid_vertex) {
        graph_.vertex_range = std::make_pair(first_vertex, first_invalid_vertex);
    }

    void FilterDuplicateEdges();

    GraphRepresentation desired_representation_;
    Graph               graph_;

private:
    void Reset();
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
    void EnsureSquarePowerOfTwoChunkSize(PGeneratorConfig& config, PEID size, bool output) const;
    void EnsureCubicPowerOfTwoChunkSize(PGeneratorConfig& config, PEID size, bool output) const;
    void EnsureOneChunkPerPE(PGeneratorConfig& config, PEID size) const;
};
} // namespace kagen
