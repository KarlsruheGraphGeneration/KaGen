#include "kagen/io/metis.h"

#include "kagen/io/buffered_writer.h"

namespace kagen {
MetisWriter::MetisWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size)
    : StandardGraphWriter(config, graph, info, rank, size) {}

void MetisWriter::WriteHeader(const std::string& filename, const SInt n, const SInt m) {
    BufferedTextOutput<> out(tag::append, filename);
    out.WriteInt(n).WriteChar(' ').WriteInt(m / 2);

    const bool has_vertex_weights = info_.has_vertex_weights;
    const bool has_edge_weights   = info_.has_edge_weights;
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

bool MetisWriter::WriteBody(const std::string& filename) {
    graph_.SortEdgelist();

    BufferedTextOutput<> out(tag::append, filename);

    SInt cur_edge = 0;
    for (SInt from = graph_.vertex_range.first; from < graph_.vertex_range.second; ++from) {
        if (info_.has_vertex_weights) {
            out.WriteInt(graph_.vertex_weights[from - graph_.vertex_range.first]).WriteChar(' ').Flush();
        }

        while (cur_edge < graph_.edges.size() && std::get<0>(graph_.edges[cur_edge]) == from) {
            const SInt to = std::get<1>(graph_.edges[cur_edge]) + 1;
            out.WriteInt(to).WriteChar(' ');
            if (info_.has_edge_weights) {
                out.WriteInt(graph_.edge_weights[cur_edge]).WriteChar(' ');
            }
            out.Flush();
            ++cur_edge;
        }
        out.WriteChar('\n').Flush();
    }

    return false;
}

namespace {
struct Format {
    std::uint64_t number_of_nodes  = 0;
    std::uint64_t number_of_edges  = 0;
    bool          has_node_weights = false;
    bool          has_edge_weights = false;
};

inline Format ParseHeader(MappedFileToker& toker) {
    toker.SkipSpaces();
    while (toker.TestChar('%')) {
        toker.SkipLine();
        toker.SkipSpaces();
    }

    const std::uint64_t number_of_nodes = toker.ScanUnsigned();
    const std::uint64_t number_of_edges = toker.ScanUnsigned() * 2;
    const std::uint64_t format          = (toker.Current() != '\n') ? toker.ScanUnsigned() : 0;
    if (!toker.ConsumeChar('\n')) {
        throw IOError("unexpected char");
    }

    if (format != 0 && format != 1 && format != 10 && format != 11 && format && format != 100 && format != 110
        && format != 101 && format != 111) {
        throw IOError("invalid or unsupported METIS graph format");
    }

    const bool has_node_weights = (format % 100) / 10; // == x1x
    const bool has_edge_weights = format % 10;         // == xx1

    return {number_of_nodes, number_of_edges, has_node_weights, has_edge_weights};
}

template <typename NodeCB, typename EdgeCB>
std::size_t ParseBody(
    MappedFileToker& toker, NodeCB&& node_cb, EdgeCB&& edge_cb, const SInt num_nodes, const bool has_node_weights,
    const bool has_edge_weights) {
    static_assert(std::is_invocable_v<NodeCB, std::uint64_t>);
    static_assert(std::is_invocable_v<EdgeCB, std::uint64_t, std::uint64_t>);

    bool exited_preemptively = false;
    for (std::uint64_t u = 0; u < num_nodes; ++u) {
        toker.SkipSpaces();
        while (toker.TestChar('%')) {
            toker.SkipLine();
            toker.SkipSpaces();
        }

        std::uint64_t node_weight = 1;
        toker.Mark();
        if (has_node_weights) {
            node_weight = toker.ScanUnsigned();
        }
        if (!node_cb(node_weight)) {
            exited_preemptively = true;
            break;
        }

        while (std::isdigit(toker.Current())) {
            const std::uint64_t v           = toker.ScanUnsigned() - 1;
            std::uint64_t       edge_weight = 1;
            if (has_edge_weights) {
                edge_weight = toker.ScanUnsigned();
            }
            edge_cb(edge_weight, v);
        }

        if (toker.ValidPosition() && !toker.ConsumeChar('\n')) {
            throw IOError("unexpected char in neighbors list");
        }
    }

    if (!exited_preemptively) {
        while (toker.TestChar('%')) {
            toker.SkipLine();
        }

        if (toker.ValidPosition()) {
            std::cerr << "Warning: ignoring extra lines at the end of the input file\n";
        }
    }

    return toker.Marked();
}

template <typename FormatCB, typename NodeCB, typename EdgeCB>
std::size_t Parse(MappedFileToker& toker, FormatCB&& format_cb, NodeCB&& node_cb, EdgeCB&& edge_cb) {
    static_assert(std::is_invocable_v<FormatCB, Format>);

    const Format format = ParseHeader(toker);
    format_cb(format);
    return ParseBody<NodeCB, EdgeCB>(
        toker, std::forward<NodeCB>(node_cb), std::forward<EdgeCB>(edge_cb), format.number_of_nodes,
        format.has_node_weights, format.has_edge_weights);
}
} // namespace

MetisReader::MetisReader(const std::string& filename) : toker_(filename) {}

std::pair<SInt, SInt> MetisReader::ReadSize() {
    toker_.Reset();
    const auto [n, m, has_node_weights, has_edge_weights] = ParseHeader(toker_);
    return {n, m};
}

Graph MetisReader::Read(
    const SInt from_vertex, const SInt to_vertex, const SInt to_edge, const GraphRepresentation representation) {
    SInt current_node = 0;
    SInt current_edge = 0;

    toker_.Reset();
    const auto [global_n, global_m, has_node_weights, has_edge_weights] = ParseHeader(toker_);

    if (cached_first_vertex_pos_ > 0 && cached_first_vertex_ <= from_vertex) {
        current_node = cached_first_vertex_;
        current_edge = cached_first_edge_;
        toker_.Seek(cached_first_vertex_pos_);
    }

    Graph       graph;
    std::size_t end_position = 0;
    if (representation == GraphRepresentation::EDGE_LIST) {
        end_position = ParseBody(
            toker_,
            [&, has_node_weights = has_node_weights](const SInt weight) {
                if (current_node >= to_vertex || current_edge >= to_edge) {
                    return false;
                }

                if (has_node_weights && current_node >= from_vertex) {
                    graph.vertex_weights.push_back(weight);
                }

                ++current_node;
                return true;
            },
            [&, has_edge_weights = has_edge_weights](const SInt weight, const SInt to) {
                ++current_edge;
                if (current_node >= from_vertex + 1) {
                    if (has_edge_weights) {
                        graph.edge_weights.push_back(weight);
                    }
                    graph.edges.emplace_back(current_node - 1, to);
                }
            },
            global_n, has_node_weights, has_edge_weights);
    } else if (representation == GraphRepresentation::CSR) {
        end_position = ParseBody(
            toker_,
            [&, has_node_weights = has_node_weights](const SInt weight) {
                if (current_node >= to_vertex || current_edge >= to_edge) {
                    return false;
                }

                if (current_node >= from_vertex) {
                    graph.xadj.push_back(graph.adjncy.size());
                    if (has_node_weights) {
                        graph.vertex_weights.push_back(weight);
                    }
                }

                ++current_node;
                return true;
            },
            [&, has_edge_weights = has_edge_weights](const SInt weight, const SInt to) {
                ++current_edge;
                if (current_node >= from_vertex + 1) {
                    if (has_edge_weights) {
                        graph.edge_weights.push_back(weight);
                    }
                    graph.adjncy.push_back(to);
                }
            },
            global_n, has_node_weights, has_edge_weights);
        graph.xadj.push_back(graph.adjncy.size());
    } else {
        __builtin_unreachable();
    }

    // Cache current position for chunked IO
    cached_first_vertex_     = current_node;
    cached_first_edge_       = current_edge;
    cached_first_vertex_pos_ = end_position;

    graph.vertex_range   = {from_vertex, current_node};
    graph.representation = representation;
    return graph;
}

SInt MetisReader::FindNodeByEdge(const SInt edge) {
    SInt current_node = 0;
    SInt current_edge = 0;

    toker_.Reset();
    const std::size_t end_position = Parse(
        toker_, [](auto) {},
        [&](SInt) {
            if (current_edge < edge) {
                ++current_node;
                return true;
            }
            return false;
        },
        [&](SInt, SInt) { ++current_edge; });

    cached_first_vertex_     = current_node;
    cached_first_edge_       = current_edge;
    cached_first_vertex_pos_ = end_position;

    return current_node;
}

std::unique_ptr<GraphReader> MetisFactory::CreateReader(const InputGraphConfig& config, PEID, PEID) const {
    return std::make_unique<MetisReader>(config.filename);
}

std::unique_ptr<GraphWriter> MetisFactory::CreateWriter(
    const OutputGraphConfig& config, Graph& graph, const GraphInfo info, const PEID rank, const PEID size) const {
    return std::make_unique<MetisWriter>(config, graph, info, rank, size);
}
} // namespace kagen
