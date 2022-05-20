#include "kagen/generators/generator.h"
#include "kagen/definitions.h"
namespace kagen {
Generator::~Generator() = default;

int Generator::Requirements() const {
    return 0;
}

bool Generator::AlmostUndirected() const {
    return false;
}

bool Generator::InvalidVertexRangeIfEmpty() const {
    return false;
}

void Generator::Generate() {
    edges_.clear();
    coordinates_.first.clear();
    coordinates_.second.clear();
    GenerateImpl();
}

const EdgeList& Generator::GetEdges() const {
    return edges_;
}

EdgeList&& Generator::TakeEdges() {
    return std::move(edges_);
}

VertexRange Generator::GetVertexRange() const {
    return vertex_range_;
}

void Generator::SetVertexRange(const VertexRange vertex_range) {
    vertex_range_ = vertex_range;
}

const Coordinates& Generator::GetCoordinates() const {
    return coordinates_;
}

Coordinates&& Generator::TakeCoordinates() {
    return std::move(coordinates_);
}
} // namespace kagen
