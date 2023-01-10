#include "kagen/generators/static/static_graph.h"

#include "kagen/generators/static/binary_parhip.h"
#include "kagen/generators/static/graph_reader.h"
#include "kagen/generators/static/metis.h"

namespace kagen {
using namespace staticgraph;

PGeneratorConfig StaticGraphFactory::NormalizeParameters(PGeneratorConfig config, PEID, bool) const {
    return config;
}

std::unique_ptr<Generator> StaticGraphFactory::Create(const PGeneratorConfig& config, PEID rank, PEID size) const {
    return std::make_unique<StaticGraph>(config, rank, size);
}

namespace {
std::unique_ptr<GraphReader> CreateReader(const PGeneratorConfig& config) {
    switch (config.static_graph.format) {
        case StaticGraphFormat::METIS:
            return std::make_unique<MetisReader>(config.static_graph.filename);

        case StaticGraphFormat::BINARY_PARHIP:
            break;
    }

    throw std::invalid_argument("bad format");
}

std::pair<SInt, SInt> ComputeRange(const SInt n, const PEID size, const PEID rank) {
    const SInt chunk = n / size;
    const SInt rem   = n % size;
    const SInt from  = rank * chunk + std::min<SInt>(rank, rem);
    const SInt to    = std::min<SInt>(from + ((static_cast<SInt>(rank) < rem) ? chunk + 1 : chunk), n);
    return {from, to};
}
} // namespace

StaticGraph::StaticGraph(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size) {}

void StaticGraph::GenerateImpl() {
    auto reader       = CreateReader(config_);
    const auto [n, m] = reader->ReadSize();

    SInt from      = 0;
    SInt to        = std::numeric_limits<SInt>::max();
    SInt num_edges = std::numeric_limits<SInt>::max();

    switch (config_.static_graph.distribution) {
        case StaticGraphDistribution::BALANCE_NODES:
            std::tie(from, to) = ComputeRange(n, size_, rank_);
            break;

        case StaticGraphDistribution::BALANCE_EDGES: {
            const auto edge_range = ComputeRange(m, size_, rank_);
            from                  = reader->FindNodeByEdge(edge_range.first);
            num_edges             = edge_range.second - edge_range.first;
            break;
        }
    }

    auto graph      = reader->Read(from, to, num_edges);
    edges_          = std::move(graph.edges);
    edge_weights_   = std::move(graph.edge_weights);
    vertex_weights_ = std::move(graph.vertex_weights);
    coordinates_    = std::move(graph.coordinates);
    SetVertexRange(from, to);
}
} // namespace kagen
