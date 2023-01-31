#include "kagen/io/hmetis.h"

#include "kagen/io/buffered_writer.h"

namespace kagen {
HMetisWriter::HMetisWriter(Graph& graph, MPI_Comm comm) : SequentialGraphWriter(graph, comm) {}

std::string HMetisWriter::DefaultExtension() const {
    return "hgr";
}

void HMetisWriter::AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) {
    BufferedTextOutput<> out(tag::append, filename);
    out.WriteInt(m / 2).WriteChar(' ').WriteInt(n);
    if (HasVertexWeights() || HasEdgeWeights()) {
        out.WriteChar(' ');
        if (HasVertexWeights()) {
            out.WriteChar('1');
        }
        if (HasEdgeWeights()) {
            out.WriteChar('1');
        } else {
            out.WriteChar('0');
        }
    }
    out.WriteChar('\n').Flush();
}

void HMetisWriter::AppendTo(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);

    for (SInt e = 0; e < edges_.size(); ++e) {
        const auto& [from, to] = edges_[e];

        if (HasEdgeWeights()) {
            out.WriteInt(static_cast<SInt>(edge_weights_[e])).WriteChar(' ');
        }

        // Edges should only occur in one direction
        if (from < to) {
            out.WriteInt(from + 1).WriteChar(' ').WriteInt(to + 1).WriteChar('\n').Flush();
        }
    }
}

void HMetisWriter::AppendFooterTo(const std::string& filename) {
    if (!HasVertexWeights()) {
        return;
    }

    BufferedTextOutput<> out(tag::append, filename);
    for (const SSInt& weight: vertex_weights_) {
        out.WriteInt(weight).WriteChar('\n').Flush();
    }
}
} // namespace kagen
