#pragma once

#include <cstdint>
#include <exception>
#include <memory>

#include "kagen/context.h"
#include "kagen/definitions.h"

namespace kagen {
enum GeneratorRequirement {
    NONE                           = 0,
    POWER_OF_TWO_COMMUNICATOR_SIZE = 1 << 0,
    SQUARE_CHUNKS                  = 1 << 1,
    CUBIC_CHUNKS                   = 1 << 2,
    ONE_CHUNK_PER_PE               = 1 << 3,
};

class Generator {
public:
    virtual ~Generator();

    virtual bool IsAlmostUndirected() const;

    void Generate();

    const EdgeList& GetEdges() const;

    EdgeList&& TakeEdges();

    VertexRange GetVertexRange() const;

    void SetVertexRange(VertexRange vetrex_range);

    const Coordinates& GetCoordinates() const;

    Coordinates&& TakeCoordinates();

protected:
    virtual void GenerateImpl() = 0;

    inline void PushCoordinate(const HPFloat x, const HPFloat y) {
        coordinates_.first.emplace_back(x, y);
    }

    inline void PushCoordinate(const HPFloat x, const HPFloat y, const HPFloat z) {
        coordinates_.second.emplace_back(x, y, z);
    }

    inline void PushEdge(const SInt from, const SInt to) {
        edges_.emplace_back(from, to);
    }

    inline void SetVertexRange(const SInt first_vertex, const SInt first_invalid_vertex) {
        vertex_range_ = std::make_pair(first_vertex, first_invalid_vertex);
    }

private:
    EdgeList    edges_;
    VertexRange vertex_range_;
    Coordinates coordinates_;
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

class GeneratorFactory {
public:
    virtual ~GeneratorFactory();

    virtual int Requirements() const;

    virtual PGeneratorConfig NormalizeParameters(PGeneratorConfig config, bool output) const;

    virtual std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const = 0;
};
} // namespace kagen
