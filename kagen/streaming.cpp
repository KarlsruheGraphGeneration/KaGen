#include "kagen/streaming.h"

#include "kagen/facade.h"
#include "kagen/io.h"

#include <fstream>
#include <iostream>

namespace kagen {
namespace {
void ForceSequentialExecution(MPI_Comm comm) {
    PEID size;
    MPI_Comm_size(comm, &size);
    if (size != 1) {
        std::cerr << "Error: streaming mode must be run sequentially\n";
        MPI_Barrier(comm);
        MPI_Abort(comm, 1);
    }
}

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

// @todo optimize
PEID FindOwner(const SInt node, const std::vector<SInt>& vertex_distribution) {
    auto it = std::upper_bound(vertex_distribution.begin(), vertex_distribution.end(), node);
    return static_cast<PEID>(std::distance(vertex_distribution.begin(), it)) - 1;
}

std::string BufferFilename(const std::string& tmp_directory, const PEID fake_rank, const PEID fake_size) {
    return tmp_directory + "/skagen_" + std::to_string(fake_rank) + "_" + std::to_string(fake_size) + ".buf";
}

void DistributeToExternalBuffers(
    const Graph& graph, const std::vector<SInt>& vertex_distribution, const int from_chunk,
    const PGeneratorConfig& config) {
    if (!graph.edge_weights.empty()) {
        throw IOError("edge weight support is not implemented");
    }

    std::vector<Edgelist> sendbufs(config.K);

    for (const auto& [from, to]: graph.edges) {
        if (config.streaming_remove_self_loops && from == to) {
            continue;
        }

        sendbufs[FindOwner(from, vertex_distribution)].emplace_back(from, to);
    }

    for (int to_chunk = 0; to_chunk < config.K; ++to_chunk) {
        std::ofstream out(BufferFilename(config.streaming_tmp_directory, from_chunk, to_chunk), std::ios::binary);
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
    const std::vector<SInt>& vertex_distribution, const int to_chunk, const PGeneratorConfig& config) {
    SInt total_size = 0;

    for (int from_chunk = 0; from_chunk < config.K; ++from_chunk) {
        std::ifstream in(BufferFilename(config.streaming_tmp_directory, from_chunk, to_chunk), std::ios::binary);
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
    for (int from_chunk = 0; from_chunk < config.K; ++from_chunk) {
        std::ifstream in(BufferFilename(config.streaming_tmp_directory, from_chunk, to_chunk), std::ios::binary);
        const SInt    size = ReadSInt(in);
        in.read(reinterpret_cast<char*>(graph.edges.data() + position), sizeof(typename Edgelist::value_type) * size);
        position += size;
    }

    return graph;
}
} // namespace

void GenerateStreamed(PGeneratorConfig config, MPI_Comm comm) {
    ForceSequentialExecution(comm);

    if (config.n == 0) {
        std::cerr << "Error: streaming mode requires the number of nodes to be given in advance\n";
        MPI_Abort(comm, 1);
    }

    if (!config.quiet && config.print_header) {
        PrintHeader(config);
    }

    GraphInfo info;
    info.global_n = config.n;

    std::vector<SInt> vertex_distribution(config.K + 1);
    for (PEID chunk = 0; chunk < config.K; ++chunk) {
        const SInt chunk_size = config.n / config.K;
        const SInt remainder  = config.n % config.K;
        const SInt from       = chunk * chunk_size + std::min<SInt>(chunk, remainder);

        vertex_distribution[chunk + 1] =
            std::min<SInt>(from + ((static_cast<SInt>(chunk) < remainder) ? chunk_size + 1 : chunk_size), config.n);
    }

    for (PEID chunk = 0; chunk < config.K; ++chunk) {
        auto factory = CreateGeneratorFactory(config.generator);
        try {
            config = factory->NormalizeParameters(config, 0, config.K, true);
        } catch (ConfigurationError& ex) {
            std::cerr << "Error: " << ex.what() << "\n";
            MPI_Barrier(comm);
            MPI_Abort(comm, 1);
        }

        auto generator = factory->Create(config, chunk, config.K);

        if (!config.quiet) {
            std::cout << "Generating chunk " << (chunk + 1) << " / " << config.K << " ... " << std::flush;
        }

        generator->Generate(GraphRepresentation::EDGE_LIST);
        Graph graph = generator->Take();

        if (config.streaming_add_reverse_edges) {
            if (!config.quiet) {
                std::cout << "adding reverse edges ... " << std::flush;
            }
            AddReverseEdges(graph);
        }

        info.global_m += graph.edges.size();
        info.has_vertex_weights |= !graph.vertex_weights.empty();
        info.has_edge_weights |= !graph.edge_weights.empty();

        DistributeToExternalBuffers(graph, vertex_distribution, chunk, config);

        if (!config.quiet) {
            std::cout << "OK" << std::endl;
        }
    }

    if (config.streaming_add_reverse_edges) {
        info.global_m = 0;

        for (int chunk = 0; chunk < config.K; ++chunk) {
            if (!config.quiet) {
                std::cout << "Counting edges (chunk " << chunk + 1 << " of " << config.K << ") ... reading ... "
                          << std::flush;
            }

            Graph graph = RestoreFromExternalBuffers(vertex_distribution, chunk, config);

            if (!config.quiet) {
                std::cout << "filtering duplicates ... " << std::flush;
            }

            RemoveMultiEdges(graph);
            info.global_m += graph.edges.size();

            if (!config.quiet) {
                std::cout << "OK" << std::endl;
            }
        }
    }

    OutputGraphConfig out_config    = config.output_graph;
    const std::string base_filename = out_config.filename;

    for (std::size_t i = 0; i < out_config.formats.size(); ++i) {
        const FileFormat& format  = out_config.formats[i];
        const auto&       factory = GetGraphFormatFactory(format);

        if ((out_config.formats.size() > 1 && !factory->DefaultExtensions().empty()) || out_config.extension) {
            out_config.filename = base_filename + "." + factory->DefaultExtensions().front();
        }

        if (std::ofstream out(out_config.filename); !out) {
            std::cerr << "Error: cannot write to " << out_config.filename << "\n";
            MPI_Abort(comm, 1);
        }

        bool continue_with_next_pass = true;
        for (int pass = 0; continue_with_next_pass; ++pass) {
            SInt offset_n = 0;
            SInt offset_m = 0;

            for (int chunk = 0; chunk < config.K; ++chunk) {
                if (!config.quiet) {
                    std::cout << "Writing " << out_config.filename << " (pass " << pass + 1 << ", chunk " << chunk + 1
                              << " of " << config.K << ") ... reading ... " << std::flush;
                }

                Graph graph = RestoreFromExternalBuffers(vertex_distribution, chunk, config);

                if (config.streaming_add_reverse_edges) {
                    if (!config.quiet) {
                        std::cout << "filtering duplicates ... " << std::flush;
                    }
                    RemoveMultiEdges(graph);
                }
                if (config.streaming_sort_edges) {
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

                continue_with_next_pass = factory->CreateWriter(out_config, graph, pass_info, chunk, config.K)
                                              ->Write(pass, out_config.filename);

                if (!continue_with_next_pass && i + 1 == out_config.formats.size()) {
                    if (!config.quiet) {
                        std::cout << "cleanup ... " << std::flush;
                    }
                    for (int from_chunk = 0; from_chunk < config.K; ++from_chunk) {
                        const std::string filename = BufferFilename(config.streaming_tmp_directory, from_chunk, chunk);
                        std::remove(filename.c_str());
                    }
                }

                if (!config.quiet) {
                    std::cout << "OK" << std::endl;
                }
            }
        }
    }
}
} // namespace kagen

