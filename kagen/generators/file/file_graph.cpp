#include "kagen/generators/file/file_graph.h"

#include "kagen/io/binary_parhip.h"
#include "kagen/io/graph_reader.h"
#include "kagen/io/metis.h"

namespace kagen {
std::unique_ptr<Generator> FileGraphFactory::Create(const PGeneratorConfig& config, PEID rank, PEID size) const {
    return std::make_unique<FileGraphGenerator>(config, rank, size);
}

namespace {
std::unique_ptr<GraphReader> CreateReader(const PGeneratorConfig& config) {
    switch (config.static_graph.format) {
        case InputFormat::METIS:
            return std::make_unique<MetisReader>(config.static_graph.filename);
        case InputFormat::PARHIP:
            return std::make_unique<BinaryParhipReader>(config.static_graph.filename);
    }

    __builtin_unreachable();
}

std::pair<SInt, SInt> ComputeRange(const SInt n, const PEID size, const PEID rank) {
    const SInt chunk = n / size;
    const SInt rem   = n % size;
    const SInt from  = rank * chunk + std::min<SInt>(rank, rem);
    const SInt to    = std::min<SInt>(from + ((static_cast<SInt>(rank) < rem) ? chunk + 1 : chunk), n);
    return {from, to};
}
} // namespace

FileGraphGenerator::FileGraphGenerator(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size) {}

void FileGraphGenerator::GenerateEdgeList() {
    GenerateImpl(GraphRepresentation::EDGE_LIST);
}

void FileGraphGenerator::GenerateCSR() {
    GenerateImpl(GraphRepresentation::CSR);
}

void FileGraphGenerator::GenerateImpl(const GraphRepresentation representation) {
    auto reader       = CreateReader(config_);
    const auto [n, m] = reader->ReadSize();

    SInt from    = 0;
    SInt to_node = std::numeric_limits<SInt>::max();
    SInt to_edge = std::numeric_limits<SInt>::max();

    switch (config_.static_graph.distribution) {
        case GraphDistribution::BALANCE_VERTICES:
            std::tie(from, to_node) = ComputeRange(n, size_, rank_);
            break;

        case GraphDistribution::BALANCE_EDGES: {
            const auto edge_range = ComputeRange(m, size_, rank_);
            from                  = reader->FindNodeByEdge(edge_range.first);
            to_edge               = edge_range.second;
            break;
        }
    }

    auto graph      = reader->Read(from, to_node, to_edge, representation);
    vertex_range_   = graph.vertex_range;
    xadj_           = std::move(graph.xadj);
    adjncy_         = std::move(graph.adjncy);
    edges_          = std::move(graph.edges);
    edge_weights_   = std::move(graph.edge_weights);
    vertex_weights_ = std::move(graph.vertex_weights);
    coordinates_    = std::move(graph.coordinates);
}
} // namespace kagen
