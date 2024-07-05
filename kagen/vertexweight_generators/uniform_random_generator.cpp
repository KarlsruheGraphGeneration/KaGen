#include "kagen/vertexweight_generators//uniform_random_generator.h"

#include <random>

namespace kagen {

UniformRandomVertexWeightGenerator::UniformRandomVertexWeightGenerator(VertexWeightConfig config, MPI_Comm comm)
    : config_(config),
      comm_{comm} {
    if (config_.weight_range_begin >= config_.weight_range_end) {
        throw std::runtime_error("Weight causes undefined behavior, need weight_range_begin > weight_range_end.");
    }
}

namespace {
VertexWeights GenerateRandomWeights(const VertexWeightConfig& config, VertexRange vertex_range, MPI_Comm comm) {
    PEID rank;
    MPI_Comm_rank(comm, &rank);
    std::mt19937                         gen((rank + 42) * 3);
    std::uniform_int_distribution<SSInt> weight_dist(config.weight_range_begin, config.weight_range_end - 1);
    const size_t                         num_vertices = vertex_range.second - vertex_range.first;
    VertexWeights                        vertex_weights(num_vertices);
    std::generate(vertex_weights.begin(), vertex_weights.end(), [&]() { return weight_dist(gen); });
    return vertex_weights;
}
} // namespace

VertexWeights
UniformRandomVertexWeightGenerator::GenerateVertexWeights(const VertexRange& vertex_range, const Edgelist&) {
    return GenerateRandomWeights(config_, vertex_range, comm_);
}
VertexWeights UniformRandomVertexWeightGenerator::GenerateVertexWeights(
    const VertexRange& vertex_range, const XadjArray&, const AdjncyArray&) {
    return GenerateRandomWeights(config_, vertex_range, comm_);
}
} // namespace kagen
