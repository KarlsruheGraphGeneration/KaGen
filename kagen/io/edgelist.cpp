#include "kagen/io/edgelist.h"

#include "kagen/definitions.h"
#include "kagen/io/buffered_writer.h"
#include "kagen/io/graph_format.h"
#include "kagen/tools/utils.h"

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

std::unique_ptr<GraphWriter>
EdgelistFactory::CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const {
    return std::make_unique<EdgelistWriter>(true, true, config, graph, comm);
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

std::unique_ptr<GraphWriter>
BinaryEdgelistFactory::CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const {
    return std::make_unique<BinaryEdgelistWriter>(true, true, config.width, config, graph, comm);
}

std::unique_ptr<GraphWriter>
UndirectedBinaryEdgelistFactory::CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const {
    return std::make_unique<BinaryEdgelistWriter>(true, false, config.width, config, graph, comm);
}

//
// Xtrapulp
//

std::unique_ptr<GraphWriter>
XtrapulpFactory::CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const {
    return std::make_unique<BinaryEdgelistWriter>(false, false, config.width, config, graph, comm);
}

//
// PlainEdgeList
//

PlainEdgelistWriter::PlainEdgelistWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm)
    : SequentialGraphWriter(config, graph, comm) {}

void PlainEdgelistWriter::AppendHeaderTo(const std::string&, SInt, SInt) {}

void PlainEdgelistWriter::AppendTo(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);
    for (const auto& [from, to]: edges_) {
        out.WriteInt(from).WriteChar('\t').WriteInt(to).WriteChar('\n').Flush();
    }
}

PlainEdgelistReader::PlainEdgelistReader(const std::string& filename, PEID rank, PEID size)
    : toker_(filename),
      rank_(rank),
      size_(size) {}

std::pair<SInt, SInt> PlainEdgelistReader::ReadSize() {
    // We don't know the graph size yet; using size_ as a fake graph size will assign one vertex to each PE, which we
    // can then ignore
    return {static_cast<SInt>(size_), static_cast<SInt>(size_)};
}

Graph PlainEdgelistReader::Read(SInt, SInt, SInt, GraphRepresentation) {
    // Start reading in the next line after from, read until the next line end after to
    const std::size_t length = toker_.Length();
    const auto [from, to]    = ComputeRange(length, size_, rank_);

    if (from > 0) {
        toker_.Seek(from - 1);
        while (toker_.ValidPosition() && toker_.Current() != '\n') {
            toker_.Advance();
        }
        toker_.Advance();
    }

    Graph graph;
    while (toker_.ValidPosition() && toker_.Position() < to) {
        const SInt u = toker_.ScanUnsigned();
        const SInt v = toker_.ScanUnsigned();
        graph.edges.emplace_back(u, v);

        if (toker_.ValidPosition() && !toker_.ConsumeChar('\n')) {
            throw IOError("unexpected char in edge list");
        }
    }
    return graph;
}

SInt PlainEdgelistReader::FindNodeByEdge(SInt) {
    return 0;
}

int PlainEdgelistReader::Deficits() const {
    return ReaderDeficits::REQUIRES_REDISTRIBUTION | ReaderDeficits::EDGE_LIST_ONLY;
}

std::unique_ptr<GraphReader>
PlainEdgelistFactory::CreateReader(const InputGraphConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<PlainEdgelistReader>(config.filename, rank, size);
}

std::unique_ptr<GraphWriter>
PlainEdgelistFactory::CreateWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm) const {
    return std::make_unique<PlainEdgelistWriter>(config, graph, comm);
}
} // namespace kagen
