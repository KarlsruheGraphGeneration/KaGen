#include "kagen/io/edgelist.h"

#include "kagen/definitions.h"
#include "kagen/io/buffered_writer.h"
#include "kagen/io/graph_format.h"

namespace kagen {
EdgelistWriter::EdgelistWriter(
    const bool header, const bool directed, const OutputGraphConfig& config, Graph& graph, MPI_Comm comm)
    : SequentialGraphWriter(config, graph, comm),
      header_(header),
      directed_(directed) {}

int EdgelistWriter::Requirements() const {
    return Requirement::NO_VERTEX_WEIGHTS | Requirement::NO_EDGE_WEIGHTS;
}

void EdgelistWriter::AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) {
    if (!header_) {
        return;
    }

    BufferedTextOutput<> out(tag::append, filename);
    out.WriteString("p ").WriteInt(n).WriteChar(' ').WriteInt(m).WriteChar('\n').Flush();
}

void EdgelistWriter::AppendTo(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);
    for (const auto& [from, to]: edges_) {
        if (!directed_ && from > to) {
            continue;
        }
        out.WriteString("e ").WriteInt(from + 1).WriteChar(' ').WriteInt(to + 1).WriteChar('\n').Flush();
    }
}

std::unique_ptr<GraphReader> EdgelistFactory::CreateReader(const InputGraphConfig&) const {
    return nullptr;
}

std::unique_ptr<GraphWriter>
EdgelistFactory::CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const {
    return std::make_unique<EdgelistWriter>(true, true, config, graph, comm);
}

std::unique_ptr<GraphReader> UndirectedEdgelistFactory::CreateReader(const InputGraphConfig&) const {
    return nullptr;
}

std::unique_ptr<GraphWriter>
UndirectedEdgelistFactory::CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const {
    return std::make_unique<EdgelistWriter>(true, false, config, graph, comm);
}

BinaryEdgelistWriter::BinaryEdgelistWriter(
    const bool header, const bool directed, const int datatype_size, const OutputGraphConfig& config, Graph& graph,
    MPI_Comm comm)
    : SequentialGraphWriter(config, graph, comm),
      header_(header),
      directed_(directed),
      datatype_size_(datatype_size) {}

int BinaryEdgelistWriter::Requirements() const {
    return Requirement::NO_VERTEX_WEIGHTS | Requirement::NO_EDGE_WEIGHTS;
}

void BinaryEdgelistWriter::AppendHeaderTo(const std::string& filename, const SInt n, const SInt m) {
    if (!header_) {
        return;
    }

    auto* fout = std::fopen(filename.c_str(), "ab");
    std::fwrite(&n, sizeof(SInt), 1, fout);
    std::fwrite(&m, sizeof(SInt), 1, fout);
    std::fclose(fout);
}

void BinaryEdgelistWriter::AppendTo(const std::string& filename) {
    auto* fout = std::fopen(filename.c_str(), "ab");
    for (const auto& [from, to]: edges_) {
        if (!directed_ && from > to) {
            continue;
        }
        if (datatype_size_ == 32) {
            const std::uint32_t edges[2] = {static_cast<std::uint32_t>(from), static_cast<std::uint32_t>(to)};
            std::fwrite(edges, sizeof(std::uint32_t), 2, fout);
        } else { // Default: 64 bit
            const std::uint64_t edges[2] = {static_cast<std::uint64_t>(from), static_cast<std::uint64_t>(to)};
            std::fwrite(edges, sizeof(std::uint64_t), 2, fout);
        }
    }
    fclose(fout);
}

std::unique_ptr<GraphReader> BinaryEdgelistFactory::CreateReader(const InputGraphConfig&) const {
    return nullptr;
}

std::unique_ptr<GraphWriter>
BinaryEdgelistFactory::CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const {
    return std::make_unique<BinaryEdgelistWriter>(true, true, config.width, config, graph, comm);
}

std::unique_ptr<GraphReader> UndirectedBinaryEdgelistFactory::CreateReader(const InputGraphConfig&) const {
    return nullptr;
}

std::unique_ptr<GraphWriter>
UndirectedBinaryEdgelistFactory::CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const {
    return std::make_unique<BinaryEdgelistWriter>(true, false, config.width, config, graph, comm);
}

std::unique_ptr<GraphReader> XtrapulpFactory::CreateReader(const InputGraphConfig&) const {
    return nullptr;
}

std::unique_ptr<GraphWriter>
XtrapulpFactory::CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const {
    return std::make_unique<BinaryEdgelistWriter>(false, false, config.width, config, graph, comm);
}
} // namespace kagen
