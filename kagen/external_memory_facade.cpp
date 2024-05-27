#include "kagen/external_memory_facade.h"

#include "kagen/definitions.h"
#include "kagen/factories.h"
#include "kagen/io.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>

namespace kagen {
namespace {
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

    for (const auto& [from, to]: edges) {
        if (to < local_vertices.first || to >= local_vertices.second) {
            additional_edges.emplace_back(to, from);
        }
    }

    edges.insert(edges.end(), additional_edges.begin(), additional_edges.end());
}

std::string BufferFilename(const PEID chunk, const PGeneratorConfig& config) {
    return config.external.tmp_directory + "/KaGen_" + std::to_string(chunk) + ".buf1";
}

std::string AggregatedBufferFilename(const PEID chunk, const PGeneratorConfig& config) {
    return config.external.tmp_directory + "/KaGen_" + std::to_string(chunk) + ".buf2";
}

void SwapoutGraphChunk(
    const Edgelist& edges, const PEID chunk, const std::vector<SInt>& distribution, const PGeneratorConfig& config) {
    std::vector<SInt> index(config.external.num_chunks + 1);

    PEID cur = 0;
    for (const auto& [from, to]: edges) {
        while (cur + 1 < config.external.num_chunks && from >= distribution[cur + 1]) {
            ++cur;
        }
        ++index[cur + 1];
    }

    std::partial_sum(index.begin(), index.end(), index.begin());

    const std::size_t edge_size = sizeof(typename Edgelist::value_type);
    const std::string filename  = BufferFilename(chunk, config);

    std::ofstream out(filename, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::ios_base::failure("cannot write to " + filename);
    }

    out.write(reinterpret_cast<const char*>(index.data()), sizeof(SInt) * index.size());
    out.write(reinterpret_cast<const char*>(edges.data()), edge_size * edges.size());
}

SInt CountEdges(const std::string& filename, const PEID chunk, const PGeneratorConfig& config) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        throw std::ios_base::failure("cannot read from " + filename);
    }

    std::vector<SInt> index(config.external.num_chunks + 1);
    in.read(reinterpret_cast<char*>(index.data()), sizeof(SInt) * index.size());

    return index[chunk + 1] - index[chunk];
}

void SwapinEdges(const std::string& filename, const PEID chunk, const PGeneratorConfig& config, Edgelist& append) {
    std::vector<SInt> index(config.external.num_chunks + 1);

    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        throw std::ios_base::failure("cannot read from " + filename);
    }

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
    }
}

void SwapoutAggregatedGraphChunk(const PEID chunk, const Graph& graph, const PGeneratorConfig& config) {
    const std::string filename = AggregatedBufferFilename(chunk, config);
    std::ofstream     out(filename, std::ios_base::binary | std::ios_base::trunc);
    if (!out) {
        throw std::ios_base::failure("cannot write to " + filename);
    }

    const SInt num_edges = graph.edges.size();
    out.write(reinterpret_cast<const char*>(&num_edges), sizeof(SInt));
    out.write(reinterpret_cast<const char*>(graph.edges.data()), sizeof(typename Edgelist::value_type) * num_edges);
}

Edgelist SwapinAggregatedEdgelist(std::ifstream& in) {
    SInt num_edges;
    in.read(reinterpret_cast<char*>(&num_edges), sizeof(SInt));

    Edgelist edges(num_edges);
    in.read(reinterpret_cast<char*>(edges.data()), sizeof(typename Edgelist::value_type) * num_edges);
    return edges;
}

Graph SwapinGraphChunk(
    const PEID chunk, const std::vector<SInt>& distribution, const PGeneratorConfig& config, bool cache_aggregate,
    const bool output_info) {
    Edgelist edges;

    const std::string aggregated_filename = AggregatedBufferFilename(chunk, config);
    if (std::ifstream in(aggregated_filename, std::ios_base::binary); in) {
        if (output_info) {
            std::cout << "reading aggregated buffer ... " << std::flush;
        }
        edges           = SwapinAggregatedEdgelist(in);
        cache_aggregate = false;
    } else {
        if (output_info) {
            std::cout << "counting unfiltered edges ... " << std::flush;
        }

        SInt num_edges = 0;
        for (int cur = 0; cur < config.external.num_chunks; ++cur) {
            const std::string filename = BufferFilename(cur, config);
            num_edges += CountEdges(filename, chunk, config);
        }

        if (output_info) {
            std::cout << "allocating(" << num_edges << ") ... " << std::flush;
        }

        edges.reserve(num_edges);

        if (output_info) {
            std::cout << "reading ... " << std::flush;
        }

        for (int cur = 0; cur < config.external.num_chunks; ++cur) {
            const std::string filename = BufferFilename(cur, config);
            SwapinEdges(filename, chunk, config, edges);
        }

        if (output_info) {
            std::cout << "filtering ... " << std::flush;
        }

        // Code calling this function expects the edges to be sorted ...
        std::sort(edges.begin(), edges.end());

        // ... and to be a "valid" graph chunk, i.e., without duplicate edges
        auto it = std::unique(edges.begin(), edges.end());
        edges.erase(it, edges.end());
    }

    Graph graph;
    graph.vertex_range.first  = distribution[chunk];
    graph.vertex_range.second = distribution[chunk + 1];
    graph.edges               = std::move(edges);

    if (cache_aggregate) {
        if (output_info) {
            std::cout << "writing aggregated buffer(" << graph.edges.size() << ") ... " << std::flush;
        }

        SwapoutAggregatedGraphChunk(chunk, graph, config);
    }

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

void RemoveBufferFiles(const PGeneratorConfig& config, const bool output_info) {
    if (output_info) {
        std::cout << "Removing buffer files ... " << std::flush;
    }
    for (PEID chunk = 0; chunk < config.external.num_chunks; ++chunk) {
        const std::string filename = BufferFilename(chunk, config);
        std::remove(filename.c_str());
    }
    if (output_info) {
        std::cout << "OK" << std::endl;
    }
}

void RemoveAggregatedBufferFiles(const PGeneratorConfig& config, const bool output_info) {
    if (output_info) {
        std::cout << "Removing aggregated buffer files ... " << std::flush;
    }
    for (PEID chunk = 0; chunk < config.external.num_chunks; ++chunk) {
        const std::string filename = AggregatedBufferFilename(chunk, config);
        std::remove(filename.c_str());
    }
    if (output_info) {
        std::cout << "OK" << std::endl;
    }
}
} // namespace

void GenerateExternalMemoryToDisk(PGeneratorConfig config, MPI_Comm comm) {
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
            std::cerr << "Error: external mode requires the number of nodes to be given in advance\n";
        }
        MPI_Abort(comm, 1);
    }

    GraphInfo local_info;
    auto      vertex_distribution = CreateVertexDistribution(config.n, config.external.num_chunks);
    auto      generator_factory   = CreateGeneratorFactory(config.generator);

    for (PEID chunk = rank; chunk < config.external.num_chunks; chunk += size) {
        try {
            config = generator_factory->NormalizeParameters(config, chunk, config.external.num_chunks, output_info);
            if (config.external.refuse_external_mode) {
                throw ConfigurationError("generator is not available in external mode");
            }
        } catch (ConfigurationError& ex) {
            if (output_error) {
                std::cerr << "Error: " << ex.what() << "\n";
            }
            MPI_Abort(comm, 1);
        }

        auto generator = generator_factory->Create(config, chunk, config.external.num_chunks);

        if (output_info) {
            std::cout << "Generating edges (" << (chunk + 1) << "... / " << config.external.num_chunks << ") ... "
                      << std::flush;
        }

        Graph graph = generator->Generate(GraphRepresentation::EDGE_LIST)->Take();
        if (local_info.has_edge_weights) {
            if (output_error) {
                std::cerr << "Error: edge weights are not supported in external mode\n";
            }
            MPI_Abort(comm, 1);
        }

        Edgelist edges = std::move(graph.edges);

        if (config.external.fix_reverse_edges) {
            if (output_info) {
                std::cout << "fixing reverse edges ... " << std::flush;
            }
            AddReverseEdges(edges);
        }

        if (config.external.fix_nonlocal_reverse_edges) {
            if (output_info) {
                std::cout << "fixing nonlocal reverse edges ... " << std::flush;
            }
            AddNonlocalReverseEdges(edges, graph.vertex_range);
        }

        // Edge count must only be correct if none of the "reverse edge fixing"-options are active; otherwise, we will
        // have to recount them during an IO run anyways
        // @todo count during final IO, then update the edge count in the file for supporting formats
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

    if (config.external.fix_reverse_edges || config.external.fix_nonlocal_reverse_edges) {
        local_info.global_m = 0;

        for (int chunk = rank; chunk < config.external.num_chunks; chunk += size) {
            if (output_info) {
                std::cout << "Counting edges (chunk " << chunk + 1 << "... / " << config.external.num_chunks << ") ... "
                          << std::flush;
            }

            Graph graph = SwapinGraphChunk(
                chunk, vertex_distribution, config, config.external.cache_aggregated_chunks, output_info);
            local_info.global_m += graph.edges.size();

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

        if (rank == ROOT && config.external.cache_aggregated_chunks) {
            RemoveBufferFiles(config, output_info);
        }
        MPI_Barrier(comm);
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
            SInt last_round_offset_n = 0;
            SInt last_round_offset_m = 0;

            const PEID rounded_chunk_count = std::ceil(1.0 * config.external.num_chunks / size) * size;

            for (PEID chunk = rank; chunk < rounded_chunk_count; chunk += size) {
                if (output_info) {
                    std::cout << "Writing " << out_config.filename << " (pass " << pass + 1 << ", chunk " << chunk + 1
                              << "... of " << config.external.num_chunks << ") ... " << std::flush;
                }

                // @todo to determine whether we want to cache the aggregated bufers, we would have to know whether we
                // need multiple IO passes or not -- extend IO interface to give this information?
                Graph graph = chunk < config.external.num_chunks
                                  ? SwapinGraphChunk(chunk, vertex_distribution, config, false, output_info)
                                  : Graph{};

                GraphInfo pass_info = global_info;
                pass_info.local_n   = graph.NumberOfLocalVertices();
                pass_info.local_m   = graph.NumberOfLocalEdges();
                MPI_Exscan(&pass_info.local_n, &pass_info.offset_n, 1, KAGEN_MPI_SINT, MPI_SUM, comm);
                MPI_Exscan(&pass_info.local_m, &pass_info.offset_m, 1, KAGEN_MPI_SINT, MPI_SUM, comm);
                pass_info.offset_n += last_round_offset_n;
                pass_info.offset_m += last_round_offset_m;

                SInt this_round_n = 0;
                SInt this_round_m = 0;
                MPI_Allreduce(&pass_info.local_n, &this_round_n, 1, KAGEN_MPI_SINT, MPI_SUM, comm);
                MPI_Allreduce(&pass_info.local_m, &this_round_m, 1, KAGEN_MPI_SINT, MPI_SUM, comm);
                last_round_offset_n += this_round_n;
                last_round_offset_m += this_round_m;

                if (output_info) {
                    std::cout << "writing(" << this_round_m << "...) ... " << std::flush;
                }

                for (PEID pe = 0; pe < size; ++pe) {
                    MPI_Barrier(comm);
                    if (pe != rank) {
                        continue;
                    }

                    continue_with_next_pass =
                        chunk < config.external.num_chunks
                        && factory->CreateWriter(out_config, graph, pass_info, chunk, config.external.num_chunks)
                               ->Write(pass, out_config.filename);
                }

                if (output_info) {
                    std::cout << "OK" << std::endl;
                }
            }

            // If not all PEs participate in the last round, we need to make sure that they all know whether we are done
            // or not
            MPI_Allreduce(MPI_IN_PLACE, &continue_with_next_pass, 1, MPI_C_BOOL, MPI_LOR, comm);
        }

        if (output_info) {
            std::cout << "Waiting for other PEs ... " << std::flush;
        }
        MPI_Barrier(comm);
        if (output_info) {
            std::cout << "OK" << std::endl;
        }
    }

    if (rank == ROOT) {
        RemoveBufferFiles(config, output_info);
        RemoveAggregatedBufferFiles(config, output_info);
    }
}
} // namespace kagen

