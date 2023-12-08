#include "app/CLI11.h"

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/io.h"
#include "kagen/kagen.h"
#include "kagen/tools/postprocessor.h"
#include "kagen/tools/utils.h"

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

    bool sort_edges = false;
};

void RemoveSelfLoops(Graph& graph) {
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

SInt FindNumberOfVertices(GraphReader& reader, Config config) {
    if (config.num_vertices > 0) {
        return config.num_vertices;
    }
    if ((reader.Deficits() & ReaderDeficits::UNKNOWN_NUM_VERTICES) == 0) {
        return reader.ReadSize().first;
    }

    std::cerr << "Error: please provide the number of vertices in the graph via the --num-vertices=<n> argument.\n";
    std::exit(1);
}

int main(int argc, char* argv[]) {
    OutputGraphConfig out_config;
    InputGraphConfig  in_config;
    Config            config;

    CLI::App app("pangraph: external graph format converter");
    app.add_option("-C,--chunks", config.num_chunks)->capture_default_str();
    app.add_option("-T,--tmp-directory", config.tmp_directory, "Directory for external memory buffers.")
        ->capture_default_str();

    app.add_option("--input-filename", in_config.filename, "Input graph")->check(CLI::ExistingFile)->required();
    app.add_option("--input-format", in_config.format, "Input graph format")
        ->transform(CLI::CheckedTransformer(GetInputFormatMap()))
        ->required()
        ->capture_default_str();
    app.add_option("--input-width", in_config.width, "Input width in bits.")->capture_default_str();
    app.add_option("--output-format", out_config.formats, "Output graph format")
        ->transform(CLI::CheckedTransformer(GetOutputFormatMap()))
        ->required()
        ->capture_default_str();
    app.add_option("--output-filename", out_config.filename, "Output graph")->required();
    app.add_option("--output-width", out_config.width, "Output width in bits.")->capture_default_str();
    app.add_flag("-q,--quiet", config.quiet, "Suppress any output to stdout.");

    app.add_flag("--remove-self-loops", config.remove_self_loops, "Remove self loops from the input graph.")
        ->capture_default_str();
    app.add_flag(
           "--add-reverse-edges", config.add_reverse_edges,
           "Add reverse edges to the input graph, such that the output graph is undirected.")
        ->capture_default_str();

    app.add_option(
        "-n,--num-vertices", config.num_vertices,
        "Providing the number of vertices can speed up the conversion for some input formats.");
    app.add_flag("--sort-edges", config.sort_edges, "Sort outgoing edges by target vertex ID.")->capture_default_str();
    CLI11_PARSE(app, argc, argv);

    const auto reader0 = CreateGraphReader(in_config.format, in_config, 0, config.num_chunks);
    GraphInfo  info;
    info.global_n = FindNumberOfVertices(*reader0, config);

    // Create output file to make sure that we can write there
    if (std::ofstream out(out_config.filename); !out) {
        std::cerr << "Error: cannot write to " << out_config.filename << "\n";
        std::exit(1);
    }

    std::vector<SInt> vertex_distribution(config.num_chunks + 1);
    for (int chunk = 0; chunk < config.num_chunks; ++chunk) {
        vertex_distribution[chunk + 1] = ComputeRange(info.global_n, config.num_chunks, chunk).second;
    }

    const auto reader        = CreateGraphReader(in_config.format, in_config, 0, 1);
    auto       reported_size = reader->ReadSize();

    Graph in_memory_graph;

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

        if (config.add_reverse_edges) {
            if (!config.quiet) {
                std::cout << "adding reverse edges ... " << std::flush;
            }
            AddReverseEdges(graph);
        }

        for (const auto& [from, to]: graph.edges) {
            if (from >= info.global_n || to >= info.global_n) {
                std::cout << "ERROR" << std::endl;
                std::cerr << "Error: edge (" << from << ", " << to << ") is out of bounds (expected at most "
                          << info.global_n << " vertices)\n";
                std::exit(1);
            }
        }

        info.global_m += graph.edges.size();
        info.has_vertex_weights |= !graph.vertex_weights.empty();
        info.has_edge_weights |= !graph.edge_weights.empty();

        if (config.num_chunks > 1) {
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
        }

        if (!config.quiet) {
            std::cout << "OK" << std::endl;
        }
    }

    if (config.add_reverse_edges) {
        info.global_m = 0;

        for (int chunk = 0; chunk < config.num_chunks; ++chunk) {
            if (!config.quiet) {
                std::cout << "Counting edges (chunk " << chunk + 1 << " of " << config.num_chunks
                          << ") ... reading ... " << std::flush;
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
                              << " of " << config.num_chunks << ") ... reading ... " << std::flush;
                }

                Graph graph = config.num_chunks > 1 ? RestoreFromExternalBuffers(vertex_distribution, chunk, config)
                                                    : std::move(in_memory_graph);

                if (config.add_reverse_edges) {
                    if (!config.quiet) {
                        std::cout << "filtering duplicates ... " << std::flush;
                    }
                    RemoveMultiEdges(graph);
                }
                if (config.sort_edges) {
                    if (!config.quiet) {
                        std::cout << "sorting ... " << std::flush;
                    }
                    std::sort(graph.edges.begin(), graph.edges.end());
                }

                if (!config.quiet) {
                    std::cout << "writing ... " << std::flush;
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

                if (!continue_with_next_pass) {
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

    return 0;
}
