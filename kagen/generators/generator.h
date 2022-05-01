#pragma once

#include <cstdint>

#include "kagen/definitions.h"

namespace kagen {
enum GeneratorRequirement {
    POWER_OF_TWO_COMMUNICATOR_SIZE = 1 << 0,
    SQUARE_CHUNKS                  = 1 << 1,
    CUBIC_CHUNKS                   = 1 << 2,
    ONE_CHUNK_PER_PE               = 1 << 3,
};

class Generator {
public:
    virtual ~Generator();

    virtual int Requirements() const;

    virtual bool AlmostUndirected() const;

    std::pair<EdgeList, VertexRange> Generate();

protected:
    virtual void GenerateImpl() = 0;

    inline void PushEdge(const SInt from, const SInt to) {
        edges_.emplace_back(from, to);
    }

    inline void SetVertexRange(const SInt first_vertex, const SInt first_invalid_vertex) {
        vertex_range_ = std::make_pair(first_vertex, first_invalid_vertex);
    }

private:
    EdgeList    edges_;
    VertexRange vertex_range_;
};
} // namespace kagen
