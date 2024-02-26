#include "kagen/io/edgelist.h"

#include "kagen/io.h"
#include "kagen/io/buffered_writer.h"
#include "kagen/io/graph_format.h"
#include "kagen/kagen.h"

namespace kagen {
EdgelistWriter::EdgelistWriter(
    const bool header, const bool directed, const OutputGraphConfig& config, Graph& graph, const GraphInfo info,
    const PEID rank, const PEID size)
    : StandardGraphWriter(config, graph, info, rank, size),
      header_(header),
      directed_(directed) {}

void EdgelistWriter::WriteHeader(const std::string& filename, const SInt n, const SInt m) {
    IgnoresVertexWeights();
    IgnoresEdgeWeights();

    if (!header_) {
        return;
    }

    BufferedTextOutput<> out(tag::append, filename);
    out.WriteString("p ").WriteInt(n).WriteChar(' ').WriteInt(m).WriteChar('\n').Flush();
}

bool EdgelistWriter::WriteBody(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);
    for (const auto& [from, to]: graph_.edges) {
        if (!directed_ && from > to) {
            continue;
        }
        out.WriteString("e ").WriteInt(from + 1).WriteChar(' ').WriteInt(to + 1).WriteChar('\n').Flush();
    }

    return false;
}

std::unique_ptr<GraphWriter> EdgelistFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<EdgelistWriter>(true, true, config, graph, info, rank, size);
}

std::unique_ptr<GraphWriter> UndirectedEdgelistFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<EdgelistWriter>(true, false, config, graph, info, rank, size);
}

BinaryEdgelistWriter::BinaryEdgelistWriter(
    const bool header, const bool directed, const int datatype_size, const OutputGraphConfig& config, Graph& graph,
    const GraphInfo info, const PEID rank, const PEID size)
    : StandardGraphWriter(config, graph, info, rank, size),
      header_(header),
      directed_(directed),
      datatype_size_(datatype_size) {}

void BinaryEdgelistWriter::WriteHeader(const std::string& filename, const SInt n, const SInt m) {
    IgnoresVertexWeights();
    IgnoresEdgeWeights();

    if (!header_) {
        return;
    }

    auto* fout = std::fopen(filename.c_str(), "ab");
    std::fwrite(&n, sizeof(SInt), 1, fout);
    std::fwrite(&m, sizeof(SInt), 1, fout);
    std::fclose(fout);
}

bool BinaryEdgelistWriter::WriteBody(const std::string& filename) {
    auto* fout = std::fopen(filename.c_str(), "ab");
    for (const auto& [from, to]: graph_.edges) {
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

    return false;
}

std::unique_ptr<GraphWriter> BinaryEdgelistFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<BinaryEdgelistWriter>(true, true, config.width, config, graph, info, rank, size);
}

std::unique_ptr<GraphWriter> UndirectedBinaryEdgelistFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<BinaryEdgelistWriter>(true, false, config.width, config, graph, info, rank, size);
}

//
// Xtrapulp
//

std::unique_ptr<GraphWriter> XtrapulpFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<BinaryEdgelistWriter>(false, false, config.width, config, graph, info, rank, size);
}

//
// PlainEdgeList
//

PlainEdgelistWriter::PlainEdgelistWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size)
    : StandardGraphWriter(config, graph, info, rank, size) {}

void PlainEdgelistWriter::WriteHeader(const std::string&, SInt, SInt) {}

bool PlainEdgelistWriter::WriteBody(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);
    for (const auto& [from, to]: graph_.edges) {
        out.WriteInt(from).WriteChar('\t').WriteInt(to).WriteChar('\n').Flush();
    }
    return false;
}

PlainEdgelistReader::PlainEdgelistReader(const std::string& filename) : toker_(filename) {}

std::pair<SInt, SInt> PlainEdgelistReader::ReadSize() {
    return {toker_.Length(), toker_.Length()};
}

Graph PlainEdgelistReader::Read(const SInt from, const SInt to, SInt, GraphRepresentation) {
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
    return ReaderDeficits::REQUIRES_REDISTRIBUTION | ReaderDeficits::EDGE_LIST_ONLY
           | ReaderDeficits::UNKNOWN_NUM_VERTICES | ReaderDeficits::UNKNOWN_NUM_EDGES;
}

std::unique_ptr<GraphReader> PlainEdgelistFactory::CreateReader(const InputGraphConfig& config, PEID, PEID) const {
    return std::make_unique<PlainEdgelistReader>(config.filename);
}

std::unique_ptr<GraphWriter> PlainEdgelistFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<PlainEdgelistWriter>(config, graph, info, rank, size);
}
} // namespace kagen
