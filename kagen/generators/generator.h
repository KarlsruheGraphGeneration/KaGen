#pragma once

#include <cstdint>

#include "kagen/definitions.h"

namespace kagen {
enum class GeneratorFeature : std::uint8_t {
    NONE              = 0,
    UNDIRECTED        = 1 << 0,
    ALMOST_UNDIRECTED = 1 << 1,
    SELF_LOOPS        = 1 << 2,
};

inline GeneratorFeature operator&(const GeneratorFeature x, const GeneratorFeature y) {
    return static_cast<GeneratorFeature>(static_cast<std::uint8_t>(x) & static_cast<std::uint8_t>(y));
}

inline GeneratorFeature operator|(const GeneratorFeature x, const GeneratorFeature y) {
    return static_cast<GeneratorFeature>(static_cast<std::uint8_t>(x) | static_cast<std::uint8_t>(y));
}

enum class GeneratorRequirement : std::uint8_t {
    NONE                           = 0,
    POWER_OF_TWO_COMMUNICATOR_SIZE = 1 << 0,
    SQAURE_CHUNKS                  = 1 << 1,
    CUBIC_CHUNKS                   = 1 << 2,
};

inline GeneratorRequirement operator&(const GeneratorRequirement x, const GeneratorRequirement y) {
    return static_cast<GeneratorRequirement>(static_cast<std::uint8_t>(x) & static_cast<std::uint8_t>(y));
}

inline GeneratorRequirement operator|(const GeneratorRequirement x, const GeneratorRequirement y) {
    return static_cast<GeneratorRequirement>(static_cast<std::uint8_t>(x) | static_cast<std::uint8_t>(y));
}

class Generator {
public:
    virtual GeneratorRequirement Requirements() const = 0;
    virtual GeneratorFeature     Features() const     = 0;

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
