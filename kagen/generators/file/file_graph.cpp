#include "kagen/generators/file/file_graph.h"

#include "kagen/io.h"
#include "kagen/tools/converter.h"
#include "kagen/tools/postprocessor.h"

namespace kagen {
std::unique_ptr<Generator> FileGraphFactory::Create(const PGeneratorConfig& config, PEID rank, PEID size) const {
    return std::make_unique<FileGraphGenerator>(config, rank, size);
}

namespace {
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

bool FileGraphGenerator::CheckDeficit(const ReaderDeficits deficit) const {
    return (deficits_ & deficit) == 1;
}

void FileGraphGenerator::GenerateImpl(const GraphRepresentation representation) {
    auto reader       = CreateGraphReader(config_.input_graph.format, config_.input_graph, rank_, size_);
    const auto [n, m] = reader->ReadSize();

    deficits_ = reader->Deficits();
    if (CheckDeficit(ReaderDeficits::REQUIRES_REDISTRIBUTION)
        && config_.input_graph.distribution == GraphDistribution::BALANCE_EDGES) {
        throw std::invalid_argument("edge balanced IO is not implemented for graph readers requiring redistribution");
    }

    // If we need postprocessing, always generate an edge list because postprocessing is not implemented for CSR
    const GraphRepresentation actual_representation =
        RequiresPostprocessing() ? GraphRepresentation::EDGE_LIST : representation;

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

    auto graph      = reader->Read(from, to_node, to_edge, actual_representation);
    vertex_range_   = graph.vertex_range;
    xadj_           = std::move(graph.xadj);
    adjncy_         = std::move(graph.adjncy);
    edges_          = std::move(graph.edges);
    edge_weights_   = std::move(graph.edge_weights);
    vertex_weights_ = std::move(graph.vertex_weights);
    coordinates_    = std::move(graph.coordinates);
}

void FileGraphGenerator::FinalizeEdgeList(MPI_Comm comm) {
    if (CheckDeficit(ReaderDeficits::REQUIRES_REDISTRIBUTION)) {
        if (CheckDeficit(ReaderDeficits::CSR_ONLY)) {
            throw std::invalid_argument("unsupported");
        }

        // Find the number of vertices in the graph
        SInt n = 0;
        for (const auto& [u, v]: edges_) {
            n = std::max(n, std::max(u, v));
        }
        MPI_Allreduce(MPI_IN_PLACE, &n, 1, KAGEN_MPI_SINT, MPI_MAX, comm);
        ++n;

        // Redistribute the graph by assigning the same number of vertices to each PE
        vertex_range_ = RedistributeEdges(edges_, edges_, true, true, n, comm);
    }

    if (config_.input_graph.add_reverse_edges) {
        std::cout << "adding reverse edges ..." << std::flush;
        AddReverseEdges(edges_, vertex_range_, comm);
    }

    if (CheckDeficit(ReaderDeficits::CSR_ONLY)) {
        FinalizeCSR(comm);
        edges_ = BuildEdgeListFromCSR(vertex_range_, xadj_, adjncy_);
        {
            XadjArray tmp;
            std::swap(xadj_, tmp);
        }
        {
            AdjncyArray tmp;
            std::swap(adjncy_, tmp);
        }
    }
}

void FileGraphGenerator::FinalizeCSR(MPI_Comm comm) {
    if (CheckDeficit(ReaderDeficits::EDGE_LIST_ONLY) || RequiresPostprocessing()) {
        FinalizeEdgeList(comm);
        std::tie(xadj_, adjncy_) = BuildCSRFromEdgeList(vertex_range_, edges_, edge_weights_);
        {
            EdgeList tmp;
            std::swap(edges_, tmp);
        }
    }
}

bool FileGraphGenerator::RequiresPostprocessing() const {
    return config_.input_graph.add_reverse_edges;
}
} // namespace kagen
