#include "app/CLI11.h"

#include "kagen/context.h"
#include "kagen/io.h"
#include "kagen/kagen.h"
#include "kagen/tools/utils.h"

#include <mpi.h>

#include <algorithm>
#include <filesystem>

using namespace kagen;

void RemoveMultiEdges(Graph& graph) {
    std::sort(graph.edges.begin(), graph.edges.end());
    auto it = std::unique(graph.edges.begin(), graph.edges.end());
    graph.edges.erase(it, graph.edges.end());
}

void AddReverseEdges(Graph& graph) {
    if (!graph.edge_weights.empty()) {
        throw std::runtime_error("not implemented");
    }

    for (auto& [from, to]: graph.edges) {
        if (from > to) {
            std::swap(from, to);
        }
    }

    RemoveMultiEdges(graph);

    const std::size_t old_size = graph.edges.size();
    for (std::size_t i = 0; i < old_size; ++i) {
        const auto& [from, to] = graph.edges[i];
        if (from != to) {
            graph.edges.emplace_back(to, from);
        }
    }
}

PEID FindOwner(const SInt node, const std::vector<SInt>& vertex_distribution) {
    auto it = std::upper_bound(vertex_distribution.begin(), vertex_distribution.end(), node);
    return static_cast<PEID>(std::distance(vertex_distribution.begin(), it)) - 1;
}

std::string BufferFilename(const std::string& tmp_directory, const PEID fake_rank, const PEID fake_size) {
    return tmp_directory + "/pangraph_" + std::to_string(fake_rank) + "_" + std::to_string(fake_size) + ".buf";
}

struct Config {
    int         num_chunks    = 1;
    std::string tmp_directory = std::filesystem::temp_directory_path();

    bool quiet = false;

    bool remove_self_loops = false;
    bool add_reverse_edges = false;

    SInt num_vertices = 0;

    bool sort_edges          = false;
    bool drop_edge_weights   = false;
    bool drop_vertex_weights = false;

    bool  weight_vertices_by_degree = false;
    SSInt weight_deg0_vertices      = 0;
};

void RemoveSelfLoops(Graph& graph) {
    if (!graph.edge_weights.empty()) {
        throw IOError("edge weight support is not implemented");
    }
    graph.edges.erase(
        std::remove_if(
            graph.edges.begin(), graph.edges.end(), [](const auto& edge) { return edge.first == edge.second; }),
        graph.edges.end());
}

void DistributeToExternalBuffers(
    const Graph& graph, const std::vector<SInt>& vertex_distribution, const int from_chunk, const Config& config) {
    if (!graph.edge_weights.empty()) {
        throw IOError("edge weight support is not implemented");
    }

    std::vector<Edgelist> sendbufs(config.num_chunks);

    for (const auto& [from, to]: graph.edges) {
        if (config.remove_self_loops && from == to) {
            continue;
        }

        sendbufs[FindOwner(from, vertex_distribution)].emplace_back(from, to);
    }

    for (int to_chunk = 0; to_chunk < config.num_chunks; ++to_chunk) {
        std::ofstream out(BufferFilename(config.tmp_directory, from_chunk, to_chunk), std::ios::binary);
        const SInt    size = sendbufs[to_chunk].size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));
        out.write(
            reinterpret_cast<const char*>(sendbufs[to_chunk].data()), sizeof(typename Edgelist::value_type) * size);
    }
}

inline SInt ReadSInt(std::ifstream& in) {
    SInt size;
    in.read(reinterpret_cast<char*>(&size), sizeof(size));
    return size;
}

Graph RestoreFromExternalBuffers(
    const std::vector<SInt>& vertex_distribution, const int to_chunk, const Config& config) {
    SInt total_size = 0;

    for (int from_chunk = 0; from_chunk < config.num_chunks; ++from_chunk) {
        std::ifstream in(BufferFilename(config.tmp_directory, from_chunk, to_chunk), std::ios::binary);
        if (!in) {
            std::cerr << "Error: buffer file does not exist\n";
            std::exit(1);
        }
        total_size += ReadSInt(in);
    }

    Graph graph;
    graph.representation      = GraphRepresentation::EDGE_LIST;
    graph.vertex_range.first  = vertex_distribution[to_chunk];
    graph.vertex_range.second = vertex_distribution[to_chunk + 1];
    graph.edges.resize(total_size);

    SInt position = 0;
    for (int from_chunk = 0; from_chunk < config.num_chunks; ++from_chunk) {
        std::ifstream in(BufferFilename(config.tmp_directory, from_chunk, to_chunk), std::ios::binary);
        const SInt    size = ReadSInt(in);
        in.read(reinterpret_cast<char*>(graph.edges.data() + position), sizeof(typename Edgelist::value_type) * size);
        position += size;
    }

    return graph;
}

SInt FindNumberOfVertices(const Config& config, const Graph& graph) {
    if (config.num_vertices > 0) {
        return config.num_vertices;
    }
    if (graph.vertex_range.first == 0 && graph.vertex_range.second > 0) {
        return graph.NumberOfLocalVertices();
    }

    SInt global_n = 0;
    for (const auto& [from, to]: graph.edges) {
        global_n = std::max(global_n, std::max(from, to));
    }
    return global_n + 1;
}

SInt FindNumberOfVertices(const InputGraphConfig& in_config, Config config) {
    if (config.num_vertices > 0) {
        return config.num_vertices;
    }

    const auto reader0 = CreateGraphReader(in_config.format, in_config, 0, config.num_chunks);
    if ((reader0->Deficits() & ReaderDeficits::UNKNOWN_NUM_VERTICES) == 0) {
        return reader0->ReadSize().first;
    }

    if (!config.quiet) {
        std::cout << "The graph format does not specify the number of vertices in the graph.\n";
        std::cout << "However, this information is required for the conversion.\n";
        std::cout << "Counting the number of vertices (skip this step by providing --num-vertices=<n>) ... "
                  << std::endl;
    }

    SInt global_n = 0;

    for (int chunk = 0; chunk < config.num_chunks; ++chunk) {
        const auto reader        = CreateGraphReader(in_config.format, in_config, 0, 1);
        const auto reported_size = reader->ReadSize();

        if (!config.quiet) {
            std::cout << "Reading " << in_config.filename << " (chunk " << chunk + 1 << " of " << config.num_chunks
                      << "): " << std::flush;
        }

        const auto [from, to] = ComputeRange(reported_size.first, config.num_chunks, chunk);
        if (!config.quiet) {
            std::cout << "[" << from << ", " << to << ") ... " << std::flush;
        }

        Graph graph = reader->Read(from, to, std::numeric_limits<SInt>::max(), GraphRepresentation::EDGE_LIST);
        for (const auto& [from, to]: graph.edges) {
            global_n = std::max(global_n, std::max(from, to));
        }

        if (!config.quiet) {
            std::cout << "OK" << std::endl;
        }
    }

    return global_n + 1;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    PEID size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 1) {
        std::cerr << "Error: must be run sequentially\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    OutputGraphConfig out_config;
    InputGraphConfig  in_config;
    Config            config;

    CLI::App app("pangraph: external graph format converter");
    app.add_option("-C,--chunks", config.num_chunks)->capture_default_str();
    app.add_option("-T,--tmp-directory", config.tmp_directory, "Directory for external memory buffers.")
        ->capture_default_str();
    app.add_flag("-q,--quiet", config.quiet, "Suppress any output to stdout.");

    auto set_all_input_widths = [&in_config](const auto width) {
        in_config.width        = width;
        in_config.vtx_width    = width;
        in_config.adjncy_width = width;
        in_config.vwgt_width   = width;
        in_config.adjwgt_width = width;
    };

    app.add_option("--input-filename", in_config.filename, "Input graph")->check(CLI::ExistingFile)->required();
    app.add_option("--input-format", in_config.format, "Input graph format")
        ->transform(CLI::CheckedTransformer(GetInputFormatMap()))
        ->capture_default_str();
    app.add_option_function<SInt>("--input-width", set_all_input_widths, "Input width in bits.")->capture_default_str();
    app.add_option("--input-vtx-width", in_config.vtx_width, "")->capture_default_str();
    app.add_option("--input-adjncy-width", in_config.adjncy_width, "")->capture_default_str();
    app.add_option("--input-vwgt-width", in_config.vwgt_width, "")->capture_default_str();
    app.add_option("--input-adjwgt-width", in_config.adjwgt_width, "")->capture_default_str();

    app.add_option(
        "-n,--num-vertices", config.num_vertices,
        "Providing the number of vertices can speed up the conversion for some input formats.");

    auto set_all_output_widths = [&out_config](const auto width) {
        out_config.width        = width;
        out_config.vtx_width    = width;
        out_config.adjncy_width = width;
        out_config.vwgt_width   = width;
        out_config.adjwgt_width = width;
    };

    app.add_option("--output-filename", out_config.filename, "Output graph")->required();
    app.add_option("--output-format", out_config.formats, "Output graph format")
        ->transform(CLI::CheckedTransformer(GetOutputFormatMap()))
        ->required()
        ->capture_default_str();
    app.add_option_function<SInt>("--output-width", set_all_output_widths, "Output width in bits.")
        ->capture_default_str();
    app.add_option("--output-vtx-width", out_config.vtx_width, "")->capture_default_str();
    app.add_option("--output-adjncy-width", out_config.adjncy_width, "")->capture_default_str();
    app.add_option("--output-vwgt-width", out_config.vwgt_width, "")->capture_default_str();
    app.add_option("--output-adjwgt-width", out_config.adjwgt_width, "")->capture_default_str();

    app.add_flag("--remove-self-loops", config.remove_self_loops, "Remove self loops from the input graph.")
        ->capture_default_str();
    app.add_flag(
           "--add-reverse-edges", config.add_reverse_edges,
           "Add reverse edges to the input graph, such that the output graph is undirected.")
        ->capture_default_str();
    app.add_flag("--sort-edges", config.sort_edges, "Sort outgoing edges by target vertex ID.")->capture_default_str();
    app.add_flag("--drop-edge-weights", config.drop_edge_weights, "Drop edge weights.")->capture_default_str();
    app.add_flag("--drop-vertex-weights", config.drop_vertex_weights, "Drop vertex weights.")->capture_default_str();

    app.add_flag(
           "--weight-vertices-by-degree", config.weight_vertices_by_degree,
           "Introduce artificial vertex weights: use degree of vertices as vertex weight.")
        ->capture_default_str();
    app.add_option(
        "--deg0-weight", config.weight_deg0_vertices, "[--weight-vertices-by-degree] Weight for isolated vertices.");
    CLI11_PARSE(app, argc, argv);

    // Create output file to make sure that we can write there -- otherwise, we might waste a lot of time with no
    // results ...
    if (std::ofstream out(out_config.filename); !out) {
        std::cerr << "Error: cannot write to " << out_config.filename << "\n";
        std::exit(1);
    }

    std::vector<SInt> vertex_distribution(config.num_chunks + 1);
    GraphInfo         info;

    // We only need to count the number of vertices in advance if we work with external memory -- otherwise, we will
    // integrate this information after reading the graph
    if (config.num_chunks > 1) {
        info.global_n = FindNumberOfVertices(in_config, config);
        for (int chunk = 0; chunk < config.num_chunks; ++chunk) {
            vertex_distribution[chunk + 1] = ComputeRange(info.global_n, config.num_chunks, chunk).second;
        }
    }

    const auto reader        = CreateGraphReader(in_config.format, in_config, 0, 1);
    const auto reported_size = reader->ReadSize();
    Graph      in_memory_graph;

    // We transform the graph from its input format to the requested output formats in up to three steps:
    //
    // (1) Read in the original graph, chunk by chunk; each chunk "simulates" a PE / MPI rank. To handle all file
    // formats, distribute the graph to external memory buffers (source chunk x target chunk): each edge is owned by the
    // MPI rank that owns its tail vertex.
    //
    // If we have to make the graph undirected (symmetric), create edges in both directions; this might create
    // multi-edges.
    //
    // (2) Read in the external memory buffers, chunk by chunk, to count the overall number of edges after removing any
    // multi-edges.
    //
    // (3) Read in the external memory buffers and write the output graphs after remving any multi-edges.

    // @todo: if the graph is already symmetric / does not have to become symmetric, and does not have to be
    // redistributed for a vertex-centric output, we could skip the redistribution buffers.

    for (int chunk = 0; chunk < config.num_chunks; ++chunk) {
        if (!config.quiet) {
            std::cout << "Reading " << in_config.filename << " (chunk " << chunk + 1 << " of " << config.num_chunks
                      << "): " << std::flush;
        }

        const auto [from, to] = ComputeRange(reported_size.first, config.num_chunks, chunk);
        if (!config.quiet) {
            std::cout << "[" << from << ", " << to << ") ... " << std::flush;
        }

        Graph graph = reader->Read(from, to, std::numeric_limits<SInt>::max(), GraphRepresentation::EDGE_LIST);

        if (config.drop_edge_weights) {
            if (!config.quiet) {
                std::cout << "dropping edge weights ... " << std::flush;
            }
            graph.edge_weights.clear();
        }

        if (config.drop_vertex_weights) {
            if (!config.quiet) {
                std::cout << "dropping vertex weights ... " << std::flush;
            }

            graph.vertex_weights.clear();
        }

        if (config.add_reverse_edges) {
            if (!config.quiet) {
                std::cout << "adding reverse edges ... " << std::flush;
            }

            AddReverseEdges(graph);
        }

        info.global_m += graph.edges.size();
        info.has_vertex_weights |= !graph.vertex_weights.empty();
        info.has_edge_weights |= !graph.edge_weights.empty();

        if (config.num_chunks > 1) {
            if (!config.quiet) {
                std::cout << "distributing to external buffers ... " << std::flush;
            }

            DistributeToExternalBuffers(graph, vertex_distribution, chunk, config);
        } else {
            if (config.remove_self_loops) {
                if (!config.quiet) {
                    std::cout << "filtering self loops ... " << std::flush;
                }
                RemoveSelfLoops(graph);
                info.global_m = graph.edges.size();
            }

            in_memory_graph = std::move(graph);

            info.global_n                = FindNumberOfVertices(config, in_memory_graph);
            in_memory_graph.vertex_range = {0, info.global_n};
            vertex_distribution.back()   = info.global_n;
        }

        if (!config.quiet) {
            std::cout << "OK" << std::endl;
        }
    }

    if (config.add_reverse_edges) {
        info.global_m = 0;

        for (int chunk = 0; chunk < config.num_chunks; ++chunk) {
            if (!config.quiet) {
                std::cout << "Counting edges (chunk " << chunk + 1 << " of " << config.num_chunks << ") ... "
                          << std::flush;
                if (config.num_chunks > 1) {
                    std::cout << "reading ... " << std::flush;
                }
            }

            Graph graph = config.num_chunks > 1 ? RestoreFromExternalBuffers(vertex_distribution, chunk, config)
                                                : std::move(in_memory_graph);

            if (!config.quiet) {
                std::cout << "filtering duplicates ... " << std::flush;
            }

            RemoveMultiEdges(graph);
            info.global_m += graph.edges.size();

            if (config.num_chunks == 1) {
                in_memory_graph = std::move(graph);
            }

            if (!config.quiet) {
                std::cout << "OK" << std::endl;
            }
        }
    }

    const std::string base_filename = out_config.filename;
    for (const FileFormat& format: out_config.formats) {
        const auto& factory = GetGraphFormatFactory(format);

        // If there are multiple output formats, append the default extension of the each file format to avoid
        // conflicts
        if (out_config.formats.size() > 1 && !factory->DefaultExtensions().empty()) {
            out_config.filename = base_filename + "." + factory->DefaultExtensions().front();
        }

        bool continue_with_next_pass = true;
        for (int pass = 0; continue_with_next_pass; ++pass) {
            SInt offset_n = 0;
            SInt offset_m = 0;

            for (int chunk = 0; chunk < config.num_chunks; ++chunk) {
                if (!config.quiet) {
                    std::cout << "Writing " << out_config.filename << " (pass " << pass + 1 << ", chunk " << chunk + 1
                              << " of " << config.num_chunks << ") ... " << std::flush;
                    if (config.num_chunks > 1) {
                        std::cout << "reading ... " << std::flush;
                    }
                }

                Graph graph = config.num_chunks > 1 ? RestoreFromExternalBuffers(vertex_distribution, chunk, config)
                                                    : std::move(in_memory_graph);

                if (config.add_reverse_edges) {
                    if (!config.quiet) {
                        std::cout << "filtering duplicates ... " << std::flush;
                    }
                    RemoveMultiEdges(graph);
                }

                if (config.weight_vertices_by_degree || config.sort_edges) {
                    if (!config.quiet) {
                        std::cout << "sorting ... " << std::flush;
                    }
                    std::sort(graph.edges.begin(), graph.edges.end());
                }

                if (!config.quiet) {
                    std::cout << "writing ... " << std::flush;
                }

                if (config.weight_vertices_by_degree) {
                    if (!config.quiet) {
                        std::cout << "weighting ... " << std::flush;
                    }

                    graph.vertex_weights.resize(graph.NumberOfLocalVertices());

                    for (SInt e = 0, u = 0; u < graph.NumberOfLocalVertices(); ++u) {
                        SInt degree = 0;
                        while (e < graph.edges.size() && graph.edges[e].first == u) {
                            ++degree;
                            ++e;
                        }

                        const SSInt weight = (degree == 0) ? config.weight_deg0_vertices : static_cast<SSInt>(degree);
                        graph.vertex_weights[u] = weight;
                    }

                    info.has_vertex_weights = true;
                }

                GraphInfo pass_info = info;
                pass_info.local_n   = graph.NumberOfLocalVertices();
                pass_info.local_m   = graph.NumberOfLocalEdges();
                pass_info.offset_n  = offset_n;
                pass_info.offset_m  = offset_m;
                offset_n += pass_info.local_n;
                offset_m += pass_info.local_m;

                continue_with_next_pass = factory->CreateWriter(out_config, graph, pass_info, chunk, config.num_chunks)
                                              ->Write(pass, out_config.filename);

                if (!continue_with_next_pass && config.num_chunks > 1) {
                    if (!config.quiet) {
                        std::cout << "cleanup ... " << std::flush;
                    }
                    for (int from_chunk = 0; from_chunk < config.num_chunks; ++from_chunk) {
                        const std::string filename = BufferFilename(config.tmp_directory, from_chunk, chunk);
                        std::remove(filename.c_str());
                    }
                }

                if (config.num_chunks == 1) {
                    in_memory_graph = std::move(graph);
                }

                if (!config.quiet) {
                    std::cout << "OK" << std::endl;
                }
            }
        }
    }

    return MPI_Finalize();
}
