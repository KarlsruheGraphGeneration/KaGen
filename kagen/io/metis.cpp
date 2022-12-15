#include "kagen/io/metis.h"

#include "kagen/io/buffered_writer.h"
#include "kagen/io/graph_writer.h"

namespace kagen {
MetisWriter::MetisWriter(Graph& graph, MPI_Comm comm) : SequentialGraphWriter(graph, comm) {}

std::string MetisWriter::DefaultExtension() const {
    return "graph";
}

int MetisWriter::Requirements() const {
    return SequentialGraphWriter::Requirement::SORTED_EDGES;
}

void MetisWriter::AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) {
    BufferedTextOutput<> out(tag::append, filename);
    out.WriteInt(n).WriteChar(' ').WriteInt(m / 2);

    const bool has_vertex_weights = HasVertexWeights();
    const bool has_edge_weights   = HasEdgeWeights();
    if (has_vertex_weights || has_edge_weights) {
        out.WriteChar(' ');
        if (has_vertex_weights) {
            out.WriteInt(has_vertex_weights).WriteInt(has_edge_weights);
        } else {
            out.WriteInt(has_edge_weights);
        }
    }
    out.WriteChar('\n').Flush();
}

void MetisWriter::AppendTo(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);

    const bool has_vertex_weights = HasVertexWeights();
    const bool has_edge_weights   = HasEdgeWeights();

    SInt cur_edge = 0;
    for (SInt from = vertex_range_.first; from < vertex_range_.second; ++from) {
        if (has_vertex_weights) {
            out.WriteInt(vertex_weights_[from - vertex_range_.first]).WriteChar(' ').Flush();
        }

        while (cur_edge < edges_.size() && std::get<0>(edges_[cur_edge]) == from) {
            const SInt to = std::get<1>(edges_[cur_edge]) + 1;
            out.WriteInt(to).WriteChar(' ');
            if (has_edge_weights) {
                out.WriteInt(edge_weights_[cur_edge]).WriteChar(' ');
            }
            out.Flush();
            ++cur_edge;
        }
        out.WriteChar('\n').Flush();
    }
}
} // namespace kagen
