#include "kagen/generators/generator.h"

namespace kagen {
Generator::~Generator() = default;

int Generator::Requirements() const {
    return 0;
}

bool Generator::AlmostUndirected() const {
    return false;
}

std::pair<EdgeList, VertexRange> Generator::Generate() {
    GenerateImpl();
    return {std::move(edges_), vertex_range_};
}
} // namespace kagen
