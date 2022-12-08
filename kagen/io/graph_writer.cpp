#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "kagen/context.h"
#include "kagen/io/graph_writer.h"
#include "kagen/tools/statistics.h"

namespace kagen {
GraphWriter::GraphWriter(Graph& graph, MPI_Comm comm)
    : edges_(graph.edges),
      vertex_range_(graph.vertex_range),
      coordinates_(graph.coordinates),
      vertex_weights_(graph.vertex_weights),
      edge_weights_(graph.edge_weights),
      comm_(comm) {
    has_vertex_weights_ = vertex_weights_.size() == vertex_range_.second - vertex_range_.first;
    has_edge_weights_   = edge_weights_.size() == edges_.size();
    MPI_Allreduce(MPI_IN_PLACE, &has_vertex_weights_, 1, MPI_CXX_BOOL, MPI_LAND, comm_);
    MPI_Allreduce(MPI_IN_PLACE, &has_edge_weights_, 1, MPI_CXX_BOOL, MPI_LAND, comm_);
}

GraphWriter::~GraphWriter() = default;

bool GraphWriter::HasVertexWeights() const {
    return has_vertex_weights_;
}

bool GraphWriter::HasEdgeWeights() const {
    return has_edge_weights_;
}

SequentialGraphWriter::SequentialGraphWriter(Graph& graph, MPI_Comm comm) : GraphWriter(graph, comm) {}

void SequentialGraphWriter::Write(const PGeneratorConfig& config) {
    PEID rank, size;
    MPI_Comm_rank(comm_, &rank);
    MPI_Comm_size(comm_, &size);
    const bool output = !config.quiet && rank == ROOT;

    const std::string base_filename = config.output_file + "." + DefaultExtension();
    const std::string filename = config.output_single_file ? base_filename : base_filename + "." + std::to_string(rank);

    const bool requires_sorted_edges      = Requirements() & Requirement::SORTED_EDGES;
    const bool requires_coordinates       = Requirements() & Requirement::COORDINATES;
    const bool requires_coordinates2d     = Requirements() & Requirement::COORDINATES_2D;
    const bool requires_coordinates3d     = Requirements() & Requirement::COORDINATES_3D;
    const bool supports_no_vertex_weights = Requirement() & Requirement::NO_VERTEX_WEIGHTS;
    const bool supports_no_edge_weights   = Requirement() & Requirement::NO_EDGE_WEIGHTS;
    const bool has_coordinates2d          = coordinates_.first.size() == vertex_range_.second - vertex_range_.first;
    const bool has_coordinates3d          = coordinates_.second.size() == vertex_range_.second - vertex_range_.first;

    // Check if edges have to be sorted
    if (requires_sorted_edges) {
        if (!std::is_sorted(edges_.begin(), edges_.end())) {
            std::sort(edges_.begin(), edges_.end());
        }
    }

    // Check if coordinates are required
    const bool local_coordinates_ok = (!requires_coordinates2d || has_coordinates2d)
                                      && (!requires_coordinates3d || has_coordinates3d)
                                      && (!requires_coordinates || has_coordinates2d || has_coordinates3d);
    bool global_coordinates_ok;
    MPI_Allreduce(&local_coordinates_ok, &global_coordinates_ok, 1, MPI_CXX_BOOL, MPI_LAND, comm_);
    if (!global_coordinates_ok) {
        if (rank == ROOT) {
            std::cerr << "Output format requires coordinates, but the graph was generated without coordinates\n";
            std::cerr << "This may happen because coordinates were disabled or the graph format does not support "
                         "coordinates\n";
        }
        std::exit(1);
    }

    // Warn if we have vertex or edge weights and the output format does not support them
    if (supports_no_vertex_weights && HasVertexWeights() && output) {
        std::cout
            << "Warning: graph was generated with vertex weights, but the output format does not support vertex weights"
            << std::endl;
    }
    if (supports_no_edge_weights && HasEdgeWeights() && output) {
        std::cout
            << "Warning: graph was generated with edge weights, but the output format does not support edge weights"
            << std::endl;
    }

    // Everything OK, write graph
    CreateFile(filename);

    if (config.output_single_file) {
        const SInt n = FindNumberOfGlobalNodes(vertex_range_, comm_);
        const SInt m = FindNumberOfGlobalEdges(edges_, comm_);

        if (output) {
            std::cout << "Writing graph to " << filename << " in " << config.output_format << " format" << std::endl;
        }

        if (rank == ROOT) {
            AppendHeaderTo(filename, n, m);
        }

        for (PEID pe = 0; pe < size; ++pe) {
            if (output) {
                std::cout << "  Writing subgraph of PE " << pe << " ... " << std::flush;
            }
            if (rank == pe) {
                AppendTo(filename);
            }
            if (output) {
                std::cout << "done" << std::endl;
            }
            MPI_Barrier(comm_);
        }

        if (rank == ROOT) {
            AppendFooterTo(filename);
        }
    } else {
        const SInt n                   = vertex_range_.second - vertex_range_.first;
        const SInt m                   = edges_.size();
        const bool write_header_footer = (rank == ROOT && config.output_header == OutputHeader::ROOT)
                                         || config.output_header == OutputHeader::ALWAYS;

        if (output) {
            std::cout << "Writing graph to [" << base_filename << ".0";
            if (size > 2) {
                std::cout << ", ...";
            }
            if (size > 1) {
                std::cout << ", " << base_filename << "." << size - 1;
            }
            std::cout << "] in " << config.output_format << " format" << std::endl;
        }

        if (write_header_footer) {
            AppendHeaderTo(filename, n, m);
        }
        AppendTo(filename);
        if (write_header_footer) {
            AppendFooterTo(filename);
        }
    }
}

void SequentialGraphWriter::AppendFooterTo(const std::string&) {}

int SequentialGraphWriter::Requirements() const {
    return Requirement::NONE;
}

void SequentialGraphWriter::CreateFile(const std::string& filename) {
    std::ofstream ofs(filename);
}

NoopWriter::NoopWriter(Graph& graph, MPI_Comm comm) : GraphWriter(graph, comm) {}

std::string NoopWriter::DefaultExtension() const {
    return "";
}

void NoopWriter::Write(const PGeneratorConfig&) {}
} // namespace kagen
