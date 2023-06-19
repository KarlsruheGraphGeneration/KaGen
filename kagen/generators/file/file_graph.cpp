#include "kagen/generators/file/file_graph.h"

#include "kagen/io.h"
#include "kagen/tools/converter.h"
#include "kagen/tools/postprocessor.h"
#include "kagen/tools/utils.h"

namespace kagen {
std::unique_ptr<Generator> FileGraphFactory::Create(const PGeneratorConfig& config, PEID rank, PEID size) const {
    return std::make_unique<FileGraphGenerator>(config, rank, size);
}

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

bool FileGraphGenerator::CheckDeficit(const ReaderDeficits deficit) const {
    return (deficits_ & deficit) == 1;
}

void FileGraphGenerator::GenerateImpl(const GraphRepresentation representation) {
    auto reader       = CreateGraphReader(config_.input_graph.format, config_.input_graph, rank_, size_);
    const auto [n, m] = reader->ReadSize();

    deficits_ = reader->Deficits();
    if (CheckDeficit(ReaderDeficits::REQUIRES_REDISTRIBUTION)
        && config_.input_graph.distribution == GraphDistribution::BALANCE_EDGES) {
        throw std::invalid_argument("not implemented");
    }

    // If we need postprocessing, always generate an edge list because postprocessing is not implemented for CSR
    actual_representation_ =
        CheckDeficit(ReaderDeficits::REQUIRES_REDISTRIBUTION) ? GraphRepresentation::EDGE_LIST : representation;

    SInt from    = 0;
    SInt to_node = std::numeric_limits<SInt>::max();
    SInt to_edge = std::numeric_limits<SInt>::max();

    switch (config_.input_graph.distribution) {
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

    auto graph      = reader->Read(from, to_node, to_edge, actual_representation_);
    vertex_range_   = graph.vertex_range;
    xadj_           = std::move(graph.xadj);
    adjncy_         = std::move(graph.adjncy);
    edges_          = std::move(graph.edges);
    edge_weights_   = std::move(graph.edge_weights);
    vertex_weights_ = std::move(graph.vertex_weights);
    coordinates_    = std::move(graph.coordinates);
}

void FileGraphGenerator::FinalizeEdgeList(MPI_Comm comm) {
    if (actual_representation_ == GraphRepresentation::CSR) {
        FinalizeCSR(comm);

        if (Output()) {
            std::cout << "converting to edge list ... " << std::flush;
        }

        edges_ = BuildEdgeListFromCSR(vertex_range_, xadj_, adjncy_);
        {
            [[maybe_unused]] auto free_xadj   = std::move(xadj_);
            [[maybe_unused]] auto free_adjncy = std::move(adjncy_);
        }
    }

    if (CheckDeficit(ReaderDeficits::REQUIRES_REDISTRIBUTION)) {
        if (Output()) {
            std::cout << "redistributing edges ... " << std::flush;
        }

        // Find the number of vertices in the graph
        SInt n = 0;
        for (const auto& [u, v]: edges_) {
            n = std::max(n, std::max(u, v));
        }
        MPI_Allreduce(MPI_IN_PLACE, &n, 1, KAGEN_MPI_SINT, MPI_MAX, comm);
        ++n;

        std::tie(vertex_range_.first, vertex_range_.second) = ComputeRange(n, size_, rank_);
        AddReverseEdgesAndRedistribute(edges_, vertex_range_, false, comm);
    }
}

void FileGraphGenerator::FinalizeCSR(MPI_Comm comm) {
    if (actual_representation_ == GraphRepresentation::EDGE_LIST) {
        FinalizeEdgeList(comm);

        if (Output()) {
            std::cout << "converting to CSR ... " << std::flush;
        }

        std::tie(xadj_, adjncy_) = BuildCSRFromEdgeList(vertex_range_, edges_, edge_weights_);
        { [[maybe_unused]] auto free_edges = std::move(edges_); }
    }
}

bool FileGraphGenerator::Output() const {
    return !config_.quiet && rank_ == 0;
}
} // namespace kagen
