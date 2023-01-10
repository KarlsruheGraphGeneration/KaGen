#include "kagen/io/binary_parhip.h"

#include <algorithm>
#include <fstream>

#include <mpi.h>

#include "kagen/io/buffered_writer.h"
#include "kagen/io/graph_writer.h"
#include "kagen/tools/statistics.h"

namespace kagen {
using ParHipID     = unsigned long long;
using ParHipWeight = SSInt;

BinaryParHipWriter::BinaryParHipWriter(Graph& graph, MPI_Comm comm) : GraphWriter(graph, comm) {}

std::string BinaryParHipWriter::DefaultExtension() const {
    return "bgf";
}

void BinaryParHipWriter::Write(const PGeneratorConfig& config) {
    PEID rank = 0, size = 0;
    MPI_Comm_rank(comm_, &rank);
    MPI_Comm_size(comm_, &size);

    const bool        output   = !config.quiet && rank == ROOT;
    const std::string filename = config.output_file + "." + DefaultExtension();

    if (!config.output_single_file && output) {
        std::cout
            << "Warning: this file format does not support distributed output; writing the graph to a single file."
            << std::endl;
    }

    // Edges must be sorted in order to convert them to the CSR format
    if (!std::is_sorted(edges_.begin(), edges_.end())) {
        std::sort(edges_.begin(), edges_.end());
    }

    const ParHipID num_global_vertices = FindNumberOfGlobalNodes(vertex_range_, comm_);
    const ParHipID num_global_edges    = FindNumberOfGlobalEdges(edges_, comm_);

    // Header
    if (rank == ROOT && config.output_header != OutputHeader::NEVER) {
        // To be compatible with the original format used by ParHiP, we use the following versions:
        // 3 = no vertex or edge weights = compatible with ParHiP
        // 2 = no vertex weights, but edge weights
        // 1 = vertex weights, but no edge weights
        // 0 = vertex weights and edge weights
        // I.e., the negated "format" code used by the Metis format
        const ParHipID vertex_weights_bit = (static_cast<SInt>(HasVertexWeights()) ^ 1) << 1;
        const ParHipID edge_weights_bit   = static_cast<SInt>(HasEdgeWeights()) ^ 1;
        const ParHipID version            = vertex_weights_bit | edge_weights_bit;

        std::ofstream out(filename, std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
        out.write(reinterpret_cast<const char*>(&version), sizeof(ParHipID));
        out.write(reinterpret_cast<const char*>(&num_global_vertices), sizeof(ParHipID));
        out.write(reinterpret_cast<const char*>(&num_global_edges), sizeof(ParHipID));
    }

    // Generate the CSR data
    const SInt            num_local_vertices = vertex_range_.second - vertex_range_.first;
    std::vector<ParHipID> offset(num_local_vertices + 1);
    std::vector<ParHipID> edges(edges_.size());

    SInt cur_offset = (3 + num_global_vertices + 1) * sizeof(ParHipID); // 3 = header size
    SInt cur_edge   = 0;
    SInt cur_vertex = 0;
    for (SInt from = vertex_range_.first; from < vertex_range_.second; ++from) {
        offset[cur_vertex]         = cur_offset;
        const SInt cur_edge_before = cur_edge;

        while (cur_edge < edges_.size() && std::get<0>(edges_[cur_edge]) == from) {
            edges[cur_edge] = std::get<1>(edges_[cur_edge]);
            ++cur_edge;
        }

        const SInt degree = cur_edge - cur_edge_before;
        cur_offset += degree * sizeof(ParHipID);
        ++cur_vertex;
    }
    offset[cur_vertex] = cur_offset;

    // Write graph to binary file
    auto output_loop = [&](auto&& action) {
        for (PEID pe = 0; pe < size; ++pe) {
            if (pe == rank) {
                std::ofstream out(filename, std::ios_base::out | std::ios_base::binary | std::ios_base::app);
                action(out);
            }
            MPI_Barrier(comm_);
        }
    };

    // Write offset array
    output_loop([&offset](std::ofstream& out) {
        out.write(reinterpret_cast<const char*>(offset.data()), offset.size() * sizeof(ParHipID));
    });

    // Write edges array
    output_loop([&edges](std::ofstream& out) {
        out.write(reinterpret_cast<const char*>(edges.data()), edges.size() * sizeof(ParHipID));
    });

    // Write  weights
    static_assert(std::is_same_v<SSInt, ParHipWeight>); // @todo create copy if the data types do not match

    if (HasVertexWeights()) {
        output_loop([this](std::ofstream& out) {
            out.write(
                reinterpret_cast<const char*>(vertex_weights_.data()), vertex_weights_.size() * sizeof(ParHipWeight));
        });
    }
    if (HasEdgeWeights()) {
        output_loop([this](std::ofstream& out) {
            out.write(reinterpret_cast<const char*>(edge_weights_.data()), edge_weights_.size() * sizeof(ParHipWeight));
        });
    }
}
} // namespace kagen
