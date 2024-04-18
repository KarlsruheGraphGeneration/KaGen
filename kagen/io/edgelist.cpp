#include "kagen/io/edgelist.h"

#include "kagen/io.h"
#include "kagen/io/buffered_writer.h"
#include "kagen/io/graph_format.h"
#include "kagen/kagen.h"

#include <cstring>

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

        toker_.ConsumeChar('\n');
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

//
// WeightedBinaryEdgelist
//

WeightedBinaryEdgelistWriter::WeightedBinaryEdgelistWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size)
    : StandardGraphWriter(config, graph, info, rank, size) {}

void WeightedBinaryEdgelistWriter::WriteHeader(const std::string&, SInt, SInt) {}

bool WeightedBinaryEdgelistWriter::WriteBody(const std::string& filename) {
    std::ofstream out(filename, std::ios_base::binary | std::ios_base::app);
    for (std::size_t e = 0; e < graph_.edges.size(); ++e) {
        const auto& [from, to] = graph_.edges[e];
        const SInt weight      = info_.has_edge_weights ? graph_.edge_weights[e] : 1;

        // @todo replace with something less ub'ish
        out.write(reinterpret_cast<const char*>(&from), config_.vtx_width / 8);
        out.write(reinterpret_cast<const char*>(&to), config_.vtx_width / 8);
        out.write(reinterpret_cast<const char*>(&weight), config_.adjwgt_width / 8);
    }

    return false;
}

WeightedBinaryEdgelistReader::WeightedBinaryEdgelistReader(
    const std::string& filename, const SInt vtx_width, const SInt adjwgt_width)
    : in_(filename, std::ios::binary),
      vtx_width_(vtx_width / 8),
      adjwgt_width_(adjwgt_width / 8) {
    if (vtx_width == 1 || (vtx_width % 8 != 0)) {
        throw IOError("vertex width must be non-zero and dividable by 8");
    }
    if (adjwgt_width == 0 || (adjwgt_width % 8 != 0)) {
        throw IOError("edge weight width must be non-zero and dividable by 8");
    }
}

SInt WeightedBinaryEdgelistReader::StepSize() const {
    return 2 * vtx_width_ + adjwgt_width_;
}

std::pair<SInt, SInt> WeightedBinaryEdgelistReader::ReadSize() {
    in_.seekg(0, std::ios_base::end);
    const std::size_t size = in_.tellg() / StepSize();
    return {size, size};
}

Graph WeightedBinaryEdgelistReader::Read(const SInt from, const SInt to, SInt, GraphRepresentation) {
    // This is fine since ReadSize() reports the size of the graph as number of edges in the graph
    in_.seekg(from * StepSize(), std::ios_base::beg);

    std::vector<char> data((to - from) * StepSize());
    in_.read(data.data(), data.size());
    if (!in_) {
        throw IOError("could not read enough bytes");
    }

    Graph graph;
    for (std::size_t i = 0; i < data.size();) {
        SInt from = 0;
        std::memcpy(&from, data.data() + i, vtx_width_);
        i += vtx_width_;

        SInt to = 0;
        std::memcpy(&to, data.data() + i, vtx_width_);
        i += vtx_width_;

        SSInt weight = 0;
        std::memcpy(&weight, data.data() + i, adjwgt_width_);
        i += adjwgt_width_;

        graph.edges.emplace_back(from, to);
        graph.edge_weights.emplace_back(weight);
    }

    return graph;
}

SInt WeightedBinaryEdgelistReader::FindNodeByEdge(SInt) {
    return 0;
}

int WeightedBinaryEdgelistReader::Deficits() const {
    return ReaderDeficits::REQUIRES_REDISTRIBUTION | ReaderDeficits::EDGE_LIST_ONLY
           | ReaderDeficits::UNKNOWN_NUM_VERTICES;
}

std::unique_ptr<GraphReader>
WeightedBinaryEdgelistFactory::CreateReader(const InputGraphConfig& config, PEID, PEID) const {
    return std::make_unique<WeightedBinaryEdgelistReader>(config.filename, config.vtx_width, config.adjwgt_width);
}

std::unique_ptr<GraphWriter> WeightedBinaryEdgelistFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<WeightedBinaryEdgelistWriter>(config, graph, info, rank, size);
}
} // namespace kagen
