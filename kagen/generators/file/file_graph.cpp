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
    return deficits_ & deficit;
}

void FileGraphGenerator::GenerateImpl(const GraphRepresentation representation) {
    auto reader = CreateGraphReader(config_.input_graph.format, config_.input_graph, rank_, size_);
    deficits_   = reader->Deficits();
    graph_      = ReadGraph(*reader, representation, config_.input_graph, rank_, size_);
}

void FileGraphGenerator::FinalizeEdgeList(MPI_Comm comm) {
    if (graph_.representation == GraphRepresentation::CSR) {
        FinalizeCSR(comm);

        if (Output()) {
            std::cout << "converting to edge list ... " << std::flush;
        }

        graph_.edges = BuildEdgeListFromCSR(graph_.vertex_range, graph_.xadj, graph_.adjncy);
        graph_.FreeCSR();
    }

    graph_ = FinalizeReadGraph(deficits_, std::move(graph_), Output(), comm);
}

void FileGraphGenerator::FinalizeCSR(MPI_Comm comm) {
    if (graph_.representation == GraphRepresentation::EDGE_LIST) {
        FinalizeEdgeList(comm);

        if (Output()) {
            std::cout << "converting to CSR ... " << std::flush;
        }

        std::tie(graph_.xadj, graph_.adjncy) =
            BuildCSRFromEdgeList(graph_.vertex_range, graph_.edges, graph_.edge_weights);
        graph_.FreeEdgelist();
    }

    graph_ = FinalizeReadGraph(deficits_, std::move(graph_), Output(), comm);
}

bool FileGraphGenerator::Output() const {
    return !config_.quiet && rank_ == 0;
}
} // namespace kagen
