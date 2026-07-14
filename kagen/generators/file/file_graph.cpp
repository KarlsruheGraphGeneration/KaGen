#include "kagen/generators/file/file_graph.h"

#include "kagen/io.h"
#include "kagen/tools/converter.h"

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

void FileGraphGenerator::GenerateImpl(const GraphRepresentation representation) {
    auto reader = CreateGraphReader(config_.input_graph.format, config_.input_graph, rank_, size_);
    fragment_   = ReadGraphFragment(*reader, representation, config_.input_graph, rank_, size_);

    if (config_.input_graph.drop_vertex_weights) {
        fragment_.graph.vertex_weights.clear();
    }
    if (config_.input_graph.drop_edge_weights) {
        fragment_.graph.edge_weights.clear();
    }
}

void FileGraphGenerator::FinalizeEdgeList(MPI_Comm comm) {
    if (fragment_.graph.representation == GraphRepresentation::CSR) {
        FinalizeCSR(comm);

        if (Output()) {
            std::cout << "converting to edge list ... " << std::flush;
        }

        graph_.edges          = BuildEdgeListFromCSR(graph_.vertex_range, graph_.xadj, graph_.adjncy);
        graph_.representation = GraphRepresentation::EDGE_LIST;
        graph_.FreeCSR();
    } else {
        graph_ = FinalizeGraphFragment(std::move(fragment_), config_.input_graph, Output(), comm);
    }
}

void FileGraphGenerator::FinalizeCSR(MPI_Comm comm) {
    if (fragment_.graph.representation == GraphRepresentation::EDGE_LIST) {
        FinalizeEdgeList(comm);

        if (graph_.has_split_vertices) {
            // BuildCSRFromEdgeList(vertex_range, edges, ...) assumes every edge's tail lies within vertex_range,
            // which does not hold for a split/boundary vertex's edges (from
            // --distribution=balance-edges-strict) -- would corrupt or crash building the resulting
            // xadj/adjncy arrays. has_split_vertices is already collectively consistent across all PEs
            // (RedistributeEdgesTrueBalance Allreduces it), so every PE reaches this same conclusion and
            // throws together.
            throw ConfigurationError(
                "CSR representation requires single-PE vertex ownership and is not supported together with a "
                "graph that has split vertices (from --distribution=balance-edges-strict); use the edge list "
                "representation instead");
        }

        if (Output()) {
            std::cout << "converting to CSR ... " << std::flush;
        }

        std::tie(graph_.xadj, graph_.adjncy) =
            BuildCSRFromEdgeList(graph_.vertex_range, graph_.edges, graph_.edge_weights);
        graph_.representation = GraphRepresentation::CSR;
        graph_.FreeEdgelist();
    } else {
        graph_ = FinalizeGraphFragment(std::move(fragment_), config_.input_graph, Output(), comm);
    }
}

bool FileGraphGenerator::Output() const {
    return !config_.quiet && rank_ == 0;
}
} // namespace kagen
