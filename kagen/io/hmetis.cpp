#include "kagen/io/hmetis.h"

#include "kagen/io/buffered_writer.h"

namespace kagen {
HMetisWriter::HMetisWriter(Graph& graph, MPI_Comm comm) : SequentialGraphWriter(graph, comm) {}

std::string HMetisWriter::DefaultExtension() const {
    return "hgr";
}

void HMetisWriter::AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) {
    BufferedTextOutput<> out(tag::append, filename);
    out.WriteInt(m).WriteChar(' ').WriteInt(n).WriteChar('\n').Flush();
}

void HMetisWriter::AppendTo(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);

    for (const auto& [from, to]: edges_) {
        out.WriteInt(from + 1).WriteChar(' ').WriteInt(to + 1).WriteChar('\n').Flush();
    }
}
} // namespace kagen
