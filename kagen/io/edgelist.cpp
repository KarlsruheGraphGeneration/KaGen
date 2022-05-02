#include "kagen/io/edgelist.h"

#include "kagen/definitions.h"
#include "kagen/io/buffered_writer.h"
#include "kagen/io/graph_writer.h"

namespace kagen {
EdgeListWriter::EdgeListWriter(EdgeList& edges, const VertexRange vertex_range, Coordinates& coordinates, MPI_Comm comm)
    : SequentialGraphWriter(edges, vertex_range, coordinates, comm) {}

std::string EdgeListWriter::DefaultExtension() const {
    return "edgelist";
}

void EdgeListWriter::AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) {
    BufferedTextOutput<> out(tag::append, filename);
    out.WriteString("p ").WriteInt(n).WriteChar(' ').WriteInt(m).WriteChar('\n').Flush();
}

void EdgeListWriter::AppendTo(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);
    for (const auto& [from, to]: edges_) {
        out.WriteString("e ").WriteInt(from + 1).WriteChar(' ').WriteInt(to + 1).WriteChar('\n').Flush();
    }
}

BinaryEdgeListWriter::BinaryEdgeListWriter(
    EdgeList& edges, const VertexRange vertex_range, Coordinates& coordinates, MPI_Comm comm)
    : SequentialGraphWriter(edges, vertex_range, coordinates, comm) {}

std::string BinaryEdgeListWriter::DefaultExtension() const {
    return "binaryedgelist";
}

void BinaryEdgeListWriter::AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) {
    auto *fout = std::fopen(filename.c_str(), "ab");
    std::fwrite(&n, sizeof(SInt), 1, fout);
    std::fwrite(&m, sizeof(SInt), 1, fout);
    std::fclose(fout);
}

void BinaryEdgeListWriter::AppendTo(const std::string& filename) {
    auto *fout = std::fopen(filename.c_str(), "ab");
    for (const auto &[from, to] : edges_) {
        const SInt edges[2] = {from, to};
        std::fwrite(edges, sizeof(SInt), 2, fout);
    }
    fclose(fout);
}

} // namespace kagen