#include "kagen/io/metis.h"

#include "kagen/io/buffered_writer.h"

namespace kagen {
MetisWriter::MetisWriter(EdgeList& edges, const VertexRange vertex_range, Coordinates& coordinates, MPI_Comm comm)
    : SequentialGraphWriter(edges, vertex_range, coordinates, comm) {}

std::string MetisWriter::DefaultExtension() const {
    return "metis";
}

void MetisWriter::AppendHeaderTo(const std::string &filename, const SInt n, const SInt m) {
    BufferedTextOutput<> out(tag::append, filename);
    out.WriteInt(n).WriteChar(' ').WriteInt(m / 2).WriteChar('\n').Flush();
}

void MetisWriter::AppendTo(const std::string &filename) {
    BufferedTextOutput<> out(tag::append, filename);

    SInt cur_edge = 0;
    for (SInt from = vertex_range_.first; from < vertex_range_.second; ++from) {
        while (cur_edge < edges_.size() && std::get<0>(edges_[cur_edge]) == from) {
            const SInt to = std::get<1>(edges_[cur_edge]) + 1;
            out.WriteInt(to).WriteChar(' ').Flush();
            ++cur_edge;
        }
        out.WriteChar('\n').Flush();
    }
}
} // namespace kagen
