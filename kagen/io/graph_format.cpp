#include "kagen/io/graph_format.h"

#include <fstream>

namespace kagen {
GraphInfo::GraphInfo(const Graph& graph, MPI_Comm comm) {
    // @todo
}

GraphWriter::GraphWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size)
    : config_(config),
      info_(info),
      graph_(graph),
      rank_(rank),
      size_(size) {}

void GraphWriter::SortEdges() {
    auto cmp_from = [](const auto& lhs, const auto& rhs) {
        return std::get<0>(lhs) < std::get<0>(rhs);
    };

    if (!std::is_sorted(graph_.edges.begin(), graph_.edges.end(), cmp_from)) {
        if (!graph_.edge_weights.empty()) {
            std::vector<SSInt> indices(graph_.edges.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&](const auto& lhs, const auto& rhs) {
                return cmp_from(graph_.edges[lhs], graph_.edges[rhs]);
            });
            for (std::size_t e = 0; e < graph_.edges.size(); ++e) {
                indices[e] = graph_.edge_weights[indices[e]];
            }
            std::swap(graph_.edge_weights, indices);
        }

        std::sort(graph_.edges.begin(), graph_.edges.end(), cmp_from);
    }
}

void GraphWriter::RequiresCoordinates() const {
    const SInt local_n = graph_.vertex_range.second - graph_.vertex_range.first;
    if (graph_.coordinates.first.size() != local_n && graph_.coordinates.second.size() != local_n) {
        throw IOError("output format requires coordinates, but the graph was generated without coordinates");
    }
}

void GraphWriter::Requires2DCoordinates() const {
    if (graph_.coordinates.first.size() != graph_.vertex_range.second - graph_.vertex_range.first) {
        throw IOError("output format requires 2D coordinates, but the graph was generated without 2D coordinates");
    }
}

void GraphWriter::Requires3DCoordinates() const {
    if (graph_.coordinates.second.size() != graph_.vertex_range.second - graph_.vertex_range.first) {
        throw IOError("output format requires 3D coordinates, but the graph was generated without 3D coordinates");
    }
}

void GraphWriter::IgnoresVertexWeights() const {
    if (rank_ == ROOT && info_.has_vertex_weights) {
        std::cerr
            << "Warning: output format does not support vertex weights, but the graph was generated with vertex weights"
            << std::endl;
    }
}

void GraphWriter::IgnoresEdgeWeights() const {
    if (rank_ == ROOT && info_.has_edge_weights) {
        std::cerr
            << "Warning: output format does not support edge weights, but the graph was generated with edge weights"
            << std::endl;
    }
}

void WriteGraph(GraphWriter& writer, const OutputGraphConfig& config, const bool output, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const std::string filename = config.distributed ? config.filename + "." + std::to_string(rank) : config.filename;
    { std::ofstream out(filename); }

    if (!config.distributed) {
        if (output) {
            std::cout << "Writing graph to " << filename << " ..." << std::endl;
        }

        bool continue_with_next_pass = true;
        for (int pass = 0; continue_with_next_pass; ++pass) {
            for (PEID pe = 0; pe < size; ++pe) {
                if (output) {
                    std::cout << "  Writing subgraph of PE " << pe << " (" << pass << ") ... " << std::flush;
                }
                if (rank == pe) {
                    continue_with_next_pass = writer.Write(pass, filename);
                }
                MPI_Barrier(comm);
                if (output) {
                    std::cout << "OK" << std::endl;
                }
            }
        }
    } else {
        if (output) {
            std::cout << "Writing graph to [" << filename << ".0";
            if (size > 2) {
                std::cout << ", ...";
            }
            if (size > 1) {
                std::cout << ", " << filename << "." << size - 1;
            }
            std::cout << "] ... " << std::flush;
        }

        bool continue_with_next_pass = true;
        for (int pass = 0; continue_with_next_pass; ++pass) {
            continue_with_next_pass = writer.Write(pass, filename);
        }

        if (output) {
            std::cout << "OK" << std::endl;
        }
    }
}

StandardGraphWriter::StandardGraphWriter(
    const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size)
    : GraphWriter(config, graph, info, rank, size) {}

bool StandardGraphWriter::Write(const int pass, const std::string& filename) {
    const bool write_header_footer = [&] {
        if (config_.distributed) {
            return (rank_ == ROOT && config_.header == OutputHeader::ROOT) || config_.header == OutputHeader::ALWAYS;
        } else {
            return rank_ == ROOT;
        }
    }();

    if (pass == 0) {
        if (write_header_footer) {
            const SInt n =
                config_.distributed ? graph_.vertex_range.second - graph_.vertex_range.first : info_.global_n;
            const SInt m = config_.distributed ? graph_.edges.size() : info_.global_m;
            WriteHeader(filename, n, m);
        }
        WriteBody(filename);
        return WriteBody(filename);
    }
    if (write_header_footer) {
        WriteFooter(filename);
    }
    return false;
}

void StandardGraphWriter::WriteFooter(const std::string&) {}
} // namespace kagen
