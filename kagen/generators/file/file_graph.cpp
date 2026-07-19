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

        if (Output()) {
            std::cout << "converting to CSR ... " << std::flush;
        }

        // Split graphs are supported: extend the CSR row space down by one on a PE that holds a replica of its
        // first vertex (left-partial, local_offset == 0), whose edges sit just below the gap-free vertex_range
        // and would otherwise underflow BuildCSRFromEdgeList's `from - vertex_range.first`. This yields the
        // overlapping physically-present row space the strict CSR reader also produces; see
        // EdgeListOnlyGenerator::FinalizeCSR for the full rationale.
        VertexRange csr_range = graph_.vertex_range;
        for (const auto& partial: graph_.partial_vertices) {
            if (partial.local_offset == 0) {
                csr_range.first = partial.vertex;
            }
        }
        graph_.vertex_range = csr_range;
        std::tie(graph_.xadj, graph_.adjncy) =
            BuildCSRFromEdgeList(csr_range, graph_.edges, graph_.edge_weights);
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
