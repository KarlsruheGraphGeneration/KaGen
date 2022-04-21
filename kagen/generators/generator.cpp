#include "kagen/generators/generator.h"

namespace kagen {
std::pair<EdgeList, VertexRange> Generator::Generate() {
    GenerateImpl();
    return {std::move(edges_), vertex_range_};
}
} // namespace kagen
