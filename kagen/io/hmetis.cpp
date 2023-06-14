#include "kagen/io/hmetis.h"

#include "kagen/io/buffered_writer.h"
#include "kagen/io/seq_graph_writer.h"

namespace kagen {
HmetisWriter::HmetisWriter(const bool directed, const OutputGraphConfig& config, Graph& graph, MPI_Comm comm)
    : SequentialGraphWriter(config, graph, comm),
      directed_(directed) {}

void HmetisWriter::AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) {
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

void HmetisWriter::AppendTo(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);

    for (SInt e = 0; e < edges_.size(); ++e) {
        const auto& [from, to] = edges_[e];

        if (HasEdgeWeights()) {
            out.WriteInt(static_cast<SInt>(edge_weights_[e])).WriteChar(' ');
        }

        // Edges should only occur in one direction
        if (directed_ || from < to) {
            out.WriteInt(from + 1).WriteChar(' ').WriteInt(to + 1).WriteChar('\n').Flush();
        }
    }
}

void HmetisWriter::AppendFooterTo(const std::string& filename) {
    if (!HasVertexWeights()) {
        return;
    }

    BufferedTextOutput<> out(tag::append, filename);
    for (const SSInt& weight: vertex_weights_) {
        out.WriteInt(weight).WriteChar('\n').Flush();
    }
}

std::unique_ptr<GraphWriter>
HmetisFactory::CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const {
    return std::make_unique<HmetisWriter>(false, config, graph, comm);
}

std::unique_ptr<GraphWriter>
DirectedHmetisFactory::CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const {
    return std::make_unique<HmetisWriter>(true, config, graph, comm);
}
} // namespace kagen
