#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

#include "../CLI11.h"
#include "kagen/definitions.h"

using namespace kagen;

namespace {
struct Edge {
    SInt  u;
    SInt  v;
    SSInt weight;

    bool operator<(const SInt other_v) const {
        return v < other_v;
    }
    bool operator<=(const SInt other_v) const {
        return v <= other_v;
    }
};

template <typename ForwardIterator, typename T>
ForwardIterator BinaryFind(ForwardIterator begin, ForwardIterator end, const T& val) {
    const auto i = std::lower_bound(begin, end, val);
    return (i != end && (*i <= val)) ? i : end;
}

bool IsCommentLine(const std::string& line) {
    const auto first_letter = line.find_first_not_of(' ');
    return first_letter != std::string::npos && line[first_letter] == '%';
}
} // namespace

int main(int argc, char* argv[]) {
    std::string filename;
    bool        quiet                         = false;
    bool        use_32bit_types               = false;
    bool        allow_self_loops              = false;
    bool        allow_directed                = false;
    bool        allow_multi_edges             = false;
    bool        allow_negative_edge_weights   = false;
    bool        allow_negative_vertex_weights = false;

    CLI::App app("chkmetis");
    app.add_option("input graph", filename, "Input graph")->check(CLI::ExistingFile)->required();
    app.add_flag("-q,--quiet", quiet, "Suppress any output to stdout.");
    app.add_flag("--32", use_32bit_types, "Assume 32 bit data types for IDs and weights.");
    app.add_flag("--self-loops", allow_self_loops, "Allow self loops");
    app.add_flag("--directed", allow_directed, "Allow directed graphs.");
    app.add_flag("--multi-edges", allow_multi_edges, "Allow multi edges.");
    app.add_flag("--negative-edge-weights", allow_negative_edge_weights, "Allow negative edge weights.");
    app.add_flag("--negative-vertex-weights", allow_negative_vertex_weights, "Allow negative vertex weights.");
    CLI11_PARSE(app, argc, argv);

    bool nonfatal_error = false;

    std::ifstream in(filename);
    std::string   line;
    while (std::getline(in, line) && IsCommentLine(line)) {
        // skip line ...
    }

    SInt n      = 0;
    SInt m      = 0;
    SInt format = 0;
    std::stringstream(line) >> n >> m >> format;

    if (n < 0) {
        if (!quiet) {
            std::cerr << "Number of vertices cannot be negative: " << n << " (stopping)\n";
        }
        std::exit(1);
    }
    if (m < 0) {
        if (!quiet) {
            std::cerr << "Number of edges cannot be negative: " << m << " (stopping)\n";
        }
        std::exit(1);
    }
    if (m > n * (n - 1) / 2) {
        if (!quiet) {
            std::cerr << "There are too many edges in the graph: with " << n << " vertices, there can be at most "
                      << (n * (n - 1) / 2) << " undirected edges in the graph, but there are " << m
                      << " undirected edges (stopping)\n";
        }
        std::exit(1);
    }
    if (use_32bit_types
        && (n > std::numeric_limits<std::int32_t>::max() || m * 2 > std::numeric_limits<std::int32_t>::max())) {
        if (!quiet) {
            std::cerr << "The graph has too many vertices or edges for 32 bit data types (stopping)\n";
        }
        std::exit(1);
    }

    if (format != 0 && format != 1 && format != 10 && format != 11) {
        if (!quiet) {
            std::cerr << "Graph format " << format << " is unsupported, should be 0, 1, 10 or 11 (stopping)\n";
        }
        std::exit(1);
    }

    const bool has_vertex_weights = format / 10;
    const bool has_edge_weights   = format % 10;
    m *= 2;

    std::vector<Edge> edges;
    edges.reserve(m);
    SSInt total_node_weight = 0;
    SSInt total_edge_weight = 0;

    std::cout << "Reading graph ..." << std::endl;

    SInt cur_u = 0;
    while (std::getline(in, line)) {
        if (IsCommentLine(line)) {
            continue;
        }
        std::stringstream node(line);

        if (has_vertex_weights) {
            SSInt weight;
            node >> weight;
            total_node_weight += weight;

            if (!allow_negative_vertex_weights && weight < 0) {
                if (!quiet) {
                    std::cerr << "Weight of vertex " << cur_u + 1 << " cannot be negative (stoppping)\n";
                }
                std::exit(1);
            } else if (use_32bit_types && weight > std::numeric_limits<std::int32_t>::max()) {
                if (!quiet) {
                    std::cerr << "Weight of vertex " << cur_u + 1 << " is too large for 32 bit data types (stopping)\n";
                }
                std::exit(1);
            }
        }

        SInt v;
        while (node >> v) {
            --v;

            if (v < 0) {
                if (!quiet) {
                    std::cerr << "Neighbor " << v + 1 << " of vertex " << cur_u + 1
                              << " cannot be negative (stopping)\n";
                }
                std::exit(1);
            } else if (v >= n) {
                if (!quiet) {
                    std::cerr << "Neighbor " << v + 1 << " of vertex " << cur_u + 1
                              << " is larger than the number of vertices specified (stopping)\n";
                }
                std::exit(1);
            } else if (!allow_self_loops && v == cur_u) {
                if (!quiet) {
                    std::cerr << "Vertex " << cur_u + 1 << " contains a self-loop (stopping)\n";
                }
                std::exit(1);
            }

            SSInt weight{1};
            if (has_edge_weights) {
                node >> weight;

                if (!allow_negative_edge_weights && weight <= 0) {
                    if (!quiet) {
                        std::cerr << "Edge incident to vertex " << cur_u + 1
                                  << " has negative edge weight (stopping)\n";
                    }
                    std::exit(1);
                } else if (use_32bit_types && weight > std::numeric_limits<std::int32_t>::max()) {
                    if (!quiet) {
                        std::cerr << "Edge incident to vertex " << cur_u + 1
                                  << " has edge weight larger than 32 bits (stopping)\n";
                    }
                    std::exit(1);
                }
            }

            edges.push_back({cur_u, v, weight});
        }
        ++cur_u;
    }

    if (cur_u != n) {
        if (!quiet) {
            std::cerr << "Number of lines does not match the number of vertices specified in the header line: expected "
                      << n << " vertices, scanned " << cur_u << " lines\n";
        }
        nonfatal_error = true;
    }
    if (static_cast<SInt>(edges.size()) != m) {
        if (!quiet) {
            std::cerr << "Number of edges does not match the number of edges specified in the header line: expected "
                      << m << " egdes, scanned " << edges.size() << " edges\n";
        }
        nonfatal_error = true;
    }

    if (use_32bit_types && total_node_weight > std::numeric_limits<std::int32_t>::max()) {
        if (!quiet) {
            std::cerr << "The total vertex weight exceeds 32 bits\n";
        }
        nonfatal_error = true;
    }
    if (use_32bit_types && total_edge_weight > std::numeric_limits<int32_t>::max()) {
        if (!quiet) {
            std::cerr << "The total edge weight exceeds 32 bits\n";
        }
        nonfatal_error = true;
    }

    if (!allow_multi_edges || !allow_directed) {
        std::cout << "Sorting edges ..." << std::endl;

        std::sort(edges.begin(), edges.end(), [](const auto& a, const auto& b) {
            return a.u < b.u || (a.u == b.u && a.v < b.v);
        });
    }

    if (!allow_multi_edges) {
        std::cout << "Checking multi edges ..." << std::endl;

        for (std::size_t i = 1; i < edges.size(); ++i) {
            const Edge& prev = edges[i - 1];
            const Edge& cur  = edges[i];
            if (prev.u == cur.u && prev.v == cur.v) {
                if (!quiet) {
                    std::cerr << "Duplicate edge " << cur.u + 1 << " --> " << cur.v + 1 << " with weights "
                              << cur.weight << " and " << prev.weight << " (stopping)\n";
                }
                std::exit(1);
            }
        }
    }

    if (!allow_directed) {
        std::cout << "Checking reverse edges ..." << std::endl;

        std::vector<SInt> nodes(n + 1);
        for (SInt i = 0, j = 0; i < n; ++i) {
            while (j < m && edges[j].u == i) {
                ++j;
            }
            nodes[i + 1] = j;
        }
        for (const auto& [u, v, weight]: edges) {
            if (u > v) {
                continue;
            }
            const auto end      = edges.begin() + nodes[v + 1];
            const auto rev_edge = BinaryFind(edges.begin() + nodes[v], end, u);

            if (rev_edge == end) {
                if (!quiet) {
                    std::cerr << "Missing reverse edge of edge " << u + 1 << " --> " << v + 1 << " (stopping)\n";
                }
                std::exit(1);
            }
            if (weight != rev_edge->weight) {
                if (!quiet) {
                    std::cerr << "Inconsistent edge weights: edge " << u + 1 << " --> " << v + 1 << " has edge weight "
                              << weight << ", but the reverse edge has edge weight " << rev_edge->weight
                              << " (stopping)\n";
                }
                std::exit(1);
            }
        }
    }

    if (!nonfatal_error) {
        std::cout << "Graph OK" << std::endl;
    }
    return nonfatal_error;
}

