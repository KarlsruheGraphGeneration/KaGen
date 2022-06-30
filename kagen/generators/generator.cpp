#include "kagen/generators/generator.h"
#include "kagen/context.h"
#include "kagen/definitions.h"

#include <algorithm>

namespace kagen {
Generator::~Generator() = default;

bool Generator::IsAlmostUndirected() const {
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

void Generator::FilterDuplicateEdges() {
    std::sort(edges_.begin(), edges_.end());
    auto it = std::unique(edges_.begin(), edges_.end());
    edges_.erase(it, edges_.end());
}

GeneratorFactory::~GeneratorFactory() = default;

int GeneratorFactory::Requirements() const {
    return GeneratorRequirement::NONE;
}

PGeneratorConfig GeneratorFactory::NormalizeParameters(const PGeneratorConfig config, bool) const {
    return config;
}
} // namespace kagen
