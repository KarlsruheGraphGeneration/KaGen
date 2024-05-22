#include "kagen/streaming_facade.h"

#include "kagen/definitions.h"
#include "kagen/factories.h"
#include "kagen/io.h"

#include <fstream>
#include <iostream>
#include <numeric>

namespace kagen {
namespace {
void FilterMultiEdges(Edgelist& edges) {
    std::sort(edges.begin(), edges.end());
    auto it = std::unique(edges.begin(), edges.end());
    edges.erase(it, edges.end());
}

void AddReverseEdges(Edgelist& edges) {
    for (auto& [from, to]: edges) {
        if (from > to) {
            std::swap(from, to);
        }
    }

    std::sort(edges.begin(), edges.end());
    auto it = std::unique(edges.begin(), edges.end());
    edges.erase(it, edges.end());

    edges.reserve(2 * edges.size());

    const std::size_t old_size = edges.size();
    for (std::size_t i = 0; i < old_size; ++i) {
        const auto& [from, to] = edges[i];
        if (from != to) {
            edges.emplace_back(to, from);
        }
    }
}

void AddNonlocalReverseEdges(Edgelist& edges, const std::pair<SInt, SInt>& local_vertices) {
    Edgelist additional_edges;

    for (auto& [from, to]: edges) {
        if (to < local_vertices.first || to >= local_vertices.second) {
            additional_edges.emplace_back(to, from);
        }
    }

    edges.insert(edges.end(), additional_edges.begin(), additional_edges.end());
}

std::string BufferFilename(const PEID chunk, const PGeneratorConfig& config) {
    return config.streaming.tmp_directory + "/sKaGen_" + std::to_string(chunk) + ".buf";
}

void SwapoutGraphChunk(
    const Edgelist& edges, const PEID chunk, const std::vector<SInt>& distribution, const PGeneratorConfig& config) {
    std::vector<SInt> index(config.streaming.num_chunks + 1);

    PEID cur = 0;
    for (const auto& [from, to]: edges) {
        while (cur + 1 < config.streaming.num_chunks && from >= distribution[cur + 1]) {
            ++cur;
        }
        ++index[cur + 1];
    }

    std::partial_sum(index.begin(), index.end(), index.begin());

    constexpr std::size_t edge_size = sizeof(typename Edgelist::value_type);

    std::cout << "Swapout edges: #edges=" << edges.size() << ", index=[";
    for (const SInt i: index) {
        std::cout << i << ", ";
    }
    std::cout << "\b\b]\n" << std::flush;

    const std::string filename = BufferFilename(chunk, config);

    // std::cout << "Edges for buffer " << filename << ":" << std::endl;
    // for (const auto& [u, v]: edges) {
    //     std::cout << u << " -> " << v << std::endl;
    // }

    std::ofstream out(filename, std::ios::binary);
    out.write(reinterpret_cast<const char*>(index.data()), sizeof(SInt) * index.size());
    out.write(reinterpret_cast<const char*>(edges.data()), edge_size * edges.size());
}

SInt CountEdges(const std::string& filename, const PEID chunk, const PGeneratorConfig& config) {
    std::vector<SInt> index(config.streaming.num_chunks + 1);

    SInt num_edges = 0;

    std::ifstream in(filename, std::ios::binary);
    in.read(reinterpret_cast<char*>(index.data()), sizeof(SInt) * index.size());

    const SInt first_edge         = index[chunk];
    const SInt first_invalid_edge = index[chunk + 1];
    num_edges += first_invalid_edge - first_edge;

    return num_edges;
}

void SwapinEdges(const std::string& filename, const PEID chunk, const PGeneratorConfig& config, Edgelist& append) {
    std::vector<SInt> index(config.streaming.num_chunks + 1);

    std::ifstream in(filename, std::ios::binary);
    in.read(reinterpret_cast<char*>(index.data()), sizeof(SInt) * index.size());

    const SInt first_edge         = index[chunk];
    const SInt first_invalid_edge = index[chunk + 1];
    const SInt num_edges          = first_invalid_edge - first_edge;

    if (num_edges > 0) {
        constexpr std::size_t edge_size = sizeof(typename Edgelist::value_type);

        const std::size_t old_size = append.size();
        append.resize(old_size + num_edges);

        in.seekg(first_edge * edge_size, std::ios_base::cur);
        in.read(reinterpret_cast<char*>(append.data() + old_size), edge_size * num_edges);

        // std::cout << "Edges from chunk " << filename << ":" << std::endl;
        // for (const auto& [u, v]: append) {
        //     std::cout << u << " -> " << v << std::endl;
        // }
    } else {
        std::cout << "No edges from chunk " << filename << std::endl;
    }
}

Graph SwapinGraphChunk(const PEID chunk, const std::vector<SInt>& distribution, const PGeneratorConfig& config) {
    // @todo count edges, then allocate just one buffer, sort when writing

    SInt num_edges = 0;
    for (int cur = 0; cur < config.streaming.num_chunks; ++cur) {
        const std::string filename = BufferFilename(cur, config);
        num_edges += CountEdges(filename, chunk, config);
    }

    Edgelist edges;
    edges.reserve(num_edges);

    for (int cur = 0; cur < config.streaming.num_chunks; ++cur) {
        const std::string filename = BufferFilename(cur, config);
        SwapinEdges(filename, chunk, config, edges);
    }

    FilterMultiEdges(edges);

    Graph graph;
    graph.vertex_range.first  = distribution[chunk];
    graph.vertex_range.second = distribution[chunk + 1];
    graph.edges               = std::move(edges);

    return graph;
}

template <typename Int>
std::pair<Int, Int> ComputeFirstLastElement(const Int rank, const Int size, const Int elements) {
    const Int div  = elements / size;
    const Int rem  = elements % size;
    const Int from = rank * div + std::min<Int>(rank, rem);
    const Int to   = from + ((rank < rem) ? div + 1 : div);
    return {from, to};
}

std::vector<SInt> CreateVertexDistribution(const SInt n, const PEID K) {
    std::vector<SInt> vertex_distribution(K + 1);
    for (PEID chunk = 0; chunk < K; ++chunk) {
        vertex_distribution[chunk + 1] = ComputeFirstLastElement<SInt>(chunk, K, n).second;
    }
    return vertex_distribution;
}
} // namespace

void GenerateStreamedToDisk(PGeneratorConfig config, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const bool output_error = (rank == ROOT);
    const bool output_info  = (rank == ROOT && !config.quiet);

    if (output_info && config.print_header) {
        PrintHeader(config);
    }

    if (config.n == 0) {
        if (output_error) {
            std::cerr << "Error: streaming mode requires the number of nodes to be given in advance\n";
        }
        MPI_Abort(comm, 1);
    }

    GraphInfo local_info;
    auto      vertex_distribution = CreateVertexDistribution(config.n, config.streaming.num_chunks);
    auto      generator_factory   = CreateGeneratorFactory(config.generator);

    std::cout << "Vertex distribution: [";
    for (const SInt i: vertex_distribution) {
        std::cout << i << ", ";
    }
    std::cout << "\b\b]\n" << std::flush;

    for (PEID chunk = rank; chunk < config.streaming.num_chunks; chunk += size) {
        try {
            config = generator_factory->NormalizeParameters(config, chunk, config.streaming.num_chunks, output_info);
            if (config.streaming.refuse_streaming_mode) {
                throw ConfigurationError("generator is not available in streaming mode");
            }
        } catch (ConfigurationError& ex) {
            if (output_error) {
                std::cerr << "Error: " << ex.what() << "\n";
            }
            MPI_Abort(comm, 1);
        }

        auto generator = generator_factory->Create(config, chunk, config.streaming.num_chunks);

        if (output_info) {
            std::cout << "Generating edges (" << (chunk + 1) << " / " << config.streaming.num_chunks << ") ... "
                      << std::flush;
        }

        Graph graph = generator->Generate(GraphRepresentation::EDGE_LIST)->Take();
        if (local_info.has_edge_weights) {
            if (output_error) {
                std::cerr << "Error: edge weights are not supported in streaming mode\n";
            }
            MPI_Abort(comm, 1);
        }

        Edgelist edges = std::move(graph.edges);

        if (config.streaming.fix_reverse_edges) {
            if (output_info) {
                std::cout << "fixing reverse edges ... " << std::flush;
            }
            AddReverseEdges(edges);
        }

        if (config.streaming.fix_nonlocal_reverse_edges) {
            if (output_info) {
                std::cout << "fixing nonlocal reverse edges ... " << std::flush;
            }
            AddNonlocalReverseEdges(edges, graph.vertex_range);
        }

        local_info.global_m += edges.size();
        local_info.has_vertex_weights |= !graph.vertex_weights.empty();

        if (output_info) {
            std::cout << "sorting edges ... " << std::flush;
        }

        std::sort(edges.begin(), edges.end());

        if (output_info) {
            std::cout << "writing to external buffer ... " << std::flush;
        }

        SwapoutGraphChunk(edges, chunk, vertex_distribution, config);

        if (output_info) {
            std::cout << "OK" << std::endl;
        }
    }

    if (output_info) {
        std::cout << "Waiting for other PEs ... " << std::flush;
    }
    MPI_Barrier(comm);
    if (output_info) {
        std::cout << "OK" << std::endl;
    }

    if (config.streaming.fix_reverse_edges) {
        local_info.global_m = 0;

        for (int chunk = rank; chunk < config.streaming.num_chunks; chunk += size) {
            if (output_info) {
                std::cout << "Counting edges (chunk " << chunk + 1 << " / " << config.streaming.num_chunks
                          << ") ... reading ... " << std::flush;
            }

            Graph graph = SwapinGraphChunk(chunk, vertex_distribution, config);

            if (output_info) {
                std::cout << "filtering duplicates ... " << std::flush;
            }

            FilterMultiEdges(graph.edges);
            local_info.global_m += graph.edges.size();

            if (!config.quiet) {
                std::cout << "OK" << std::endl;
            }
        }

        if (output_info) {
            std::cout << "Waiting for other PEs ... " << std::flush;
        }
        MPI_Barrier(comm);
        if (output_info) {
            std::cout << "OK" << std::endl;
        }
    }

    GraphInfo global_info(local_info, comm);
    global_info.global_n = config.n;

    OutputGraphConfig out_config    = config.output_graph;
    const std::string base_filename = out_config.filename;

    for (std::size_t i = 0; i < out_config.formats.size(); ++i) {
        const FileFormat& format  = out_config.formats[i];
        const auto&       factory = GetGraphFormatFactory(format);

        if ((out_config.formats.size() > 1 && !factory->DefaultExtensions().empty()) || out_config.extension) {
            out_config.filename = base_filename + "." + factory->DefaultExtensions().front();
        }

        if (rank == ROOT) {
            if (std::ofstream out(out_config.filename); !out) {
                std::cerr << "Error: cannot write to " << out_config.filename << "\n";
                MPI_Abort(comm, 1);
            }
        }

        bool continue_with_next_pass = true;
        for (int pass = 0; continue_with_next_pass; ++pass) {
            SInt offset_n = 0;
            SInt offset_m = 0;

            for (PEID chunk = rank; chunk < config.streaming.num_chunks; chunk += size) {
                if (output_info) {
                    std::cout << "Writing " << out_config.filename << " (pass " << pass + 1 << ", chunk " << chunk + 1
                              << " of " << config.streaming.num_chunks << ") ... reading ... " << std::flush;
                }

                Graph graph = SwapinGraphChunk(chunk, vertex_distribution, config);

                if (!config.quiet) {
                    std::cout << "writing ... " << std::flush;
                }

                GraphInfo pass_info = global_info;
                pass_info.local_n   = graph.NumberOfLocalVertices();
                pass_info.local_m   = graph.NumberOfLocalEdges();
                pass_info.offset_n  = offset_n;
                pass_info.offset_m  = offset_m;
                offset_n += pass_info.local_n;
                offset_m += pass_info.local_m;

                for (PEID pe = 0; pe < size; ++pe) {
                    MPI_Barrier(comm);
                    if (pe != rank) {
                        continue;
                    }

                    continue_with_next_pass =
                        factory->CreateWriter(out_config, graph, pass_info, chunk, config.streaming.num_chunks)
                            ->Write(pass, out_config.filename);
                }

                if (!continue_with_next_pass && i + 1 == out_config.formats.size()) {
                    if (output_info) {
                        std::cout << "cleanup ... " << std::flush;
                    }
                    for (int from_chunk = 0; from_chunk < config.streaming.num_chunks; ++from_chunk) {
                        const std::string filename = BufferFilename(chunk, config);
                        std::remove(filename.c_str());
                    }
                }

                if (output_info) {
                    std::cout << "OK" << std::endl;
                }
            }
        }

        if (output_info) {
            std::cout << "Waiting for other PEs ... " << std::flush;
        }
        MPI_Barrier(comm);
        if (output_info) {
            std::cout << "OK" << std::endl;
        }
    }
}
} // namespace kagen

