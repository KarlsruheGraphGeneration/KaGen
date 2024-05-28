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

std::string BufferFilename(const PEID pe, const PGeneratorConfig& config, const bool index) {
    return config.external.tmp_directory + "/KaGen.pe" + std::to_string(pe) + "." + (index ? "index" : "edges");
}

std::string AggregatedBufferFilename(const PEID pe, const PGeneratorConfig& config) {
    return config.external.tmp_directory + "/KaGen.chunk" + std::to_string(pe) + ".aggregated";
}

SInt SwapoutGraphChunk(
    std::ofstream& out_index, std::ofstream& out_edges, const SInt offset, const Edgelist& edges,
    const std::vector<SInt>& distribution, const PGeneratorConfig& config) {
    std::vector<SInt> index(config.external.num_chunks + 1);

    PEID cur = 0;
    for (const auto& [from, to]: edges) {
        while (cur + 1 < config.external.num_chunks && from >= distribution[cur + 1]) {
            ++cur;
        }
        ++index[cur + 1];
    }

    index.front() = offset;
    std::partial_sum(index.begin(), index.end(), index.begin());

    const std::size_t edge_size = sizeof(typename Edgelist::value_type);
    out_index.write(reinterpret_cast<const char*>(index.data()), sizeof(SInt) * index.size());
    out_edges.write(reinterpret_cast<const char*>(edges.data()), edge_size * edges.size());

    return offset + edges.size();
}

void SwapinEdges(std::ifstream& in, Edgelist& append, const std::pair<SInt, SInt>& range) {
    auto [first_edge, first_invalid_edge] = range;

    const SInt num_edges = first_invalid_edge - first_edge;

    if (num_edges > 0) {
        constexpr std::size_t edge_size = sizeof(typename Edgelist::value_type);

        const std::size_t old_size = append.size();
        append.resize(old_size + num_edges);

        in.seekg(first_edge * edge_size);
        in.read(reinterpret_cast<char*>(append.data() + old_size), edge_size * num_edges);
    }
}

void SwapoutAggregatedGraphChunk(std::ofstream& out, const Graph& graph) {
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

PEID ComputeNumChunksOnPE(const PEID pe, const PEID size, const PEID num_chunks) {
    return num_chunks / size + (pe < num_chunks % size);
}

using Indices = std::vector<std::vector<std::pair<SInt, SInt>>>;

std::vector<std::vector<std::pair<SInt, SInt>>>
ReadMyIndices(const PEID rank, const PEID size, const PGeneratorConfig& config) {
    const PEID num_chunks    = config.external.num_chunks;
    const PEID my_num_chunks = ComputeNumChunksOnPE(rank, size, num_chunks);

    // Total memory across all PEs: 16 bytes * num_chunks^2.
    // This should be fine for realistic number of chunks, e.g., for num_chunks = 16 384, we need 4 GB of memory.
    Indices my_indices(my_num_chunks, std::vector<std::pair<SInt, SInt>>(num_chunks));

    for (PEID pe = 0; pe < size; ++pe) {
        const PEID pe_num_chunks = ComputeNumChunksOnPE(pe, size, num_chunks);

        const std::vector<SInt> pe_indices = [&] {
            std::vector<SInt> pe_indices(1ul * (num_chunks + 1) * pe_num_chunks);

            const std::string filename = BufferFilename(pe, config, true);
            std::ifstream     in(filename, std::ios_base::binary);
            in.read(reinterpret_cast<char*>(pe_indices.data()), sizeof(SInt) * pe_indices.size());

            return pe_indices;
        }();

        for (PEID to_chunk = rank; to_chunk < num_chunks; to_chunk += size) {
            for (PEID from_chunk = pe; from_chunk < num_chunks; from_chunk += size) {
                const PEID nth_chunk_on_pe = from_chunk / size;

                const SInt i = nth_chunk_on_pe * (num_chunks + 1) + to_chunk;

                my_indices[to_chunk / size][from_chunk] = {pe_indices[i], pe_indices[i + 1]};
            }
        }
    }

    return my_indices;
}

Graph SwapinGraphChunk(
    const PEID chunk, const Indices& my_indices, const std::vector<SInt>& distribution, std::ifstream* aggregate_in,
    std::ofstream* aggregate_out, const PGeneratorConfig& config, MPI_Comm comm) {
    PEID size;
    PEID rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    const bool output_info = (rank == ROOT && !config.quiet);

    Edgelist edges;

    if (aggregate_in) {
        if (output_info) {
            std::cout << "reading aggregated buffer ... " << std::flush;
        }
        edges = SwapinAggregatedEdgelist(*aggregate_in);
    } else {
        const PEID my_nth_chunk = chunk / size;
        SInt       num_edges    = 0;
        for (PEID from_chunk = 0; from_chunk < config.external.num_chunks; ++from_chunk) {
            const auto& [first_edge, first_invalid_edge] = my_indices[my_nth_chunk][from_chunk];
            num_edges += first_invalid_edge - first_edge;
        }

        if (output_info) {
            std::cout << "allocating ... " << std::flush;
        }

        edges.reserve(num_edges);

        if (output_info) {
            std::cout << "reading ... " << std::flush;
        }

        {
            std::vector<std::ifstream> buffer_ins;
            buffer_ins.reserve(size);
            for (PEID pe = 0; pe < size; ++pe) {
                buffer_ins.emplace_back(BufferFilename(pe, config, false), std::ios_base::binary);
            }

            for (PEID from_chunk = 0; from_chunk < config.external.num_chunks; ++from_chunk) {
                const PEID pe = from_chunk % size;
                SwapinEdges(buffer_ins[pe], edges, my_indices[my_nth_chunk][from_chunk]);
            }
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
    graph.representation      = GraphRepresentation::EDGE_LIST;
    graph.vertex_range.first  = distribution[chunk];
    graph.vertex_range.second = distribution[chunk + 1];
    graph.edges               = std::move(edges);

    if (!aggregate_in && aggregate_out) {
        if (output_info) {
            std::cout << "writing aggregated buffer ... " << std::flush;
        }

        SwapoutAggregatedGraphChunk(*aggregate_out, graph);
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

void RemoveBufferFiles(const PEID rank, const PGeneratorConfig& config, const bool output_info) {
    if (output_info) {
        std::cout << "Removing buffer files ... " << std::flush;
    }
    const std::string filename_index = BufferFilename(rank, config, true);
    std::remove(filename_index.c_str());
    const std::string filename_edges = BufferFilename(rank, config, false);
    std::remove(filename_edges.c_str());
    if (output_info) {
        std::cout << "OK" << std::endl;
    }
}

void RemoveAggregatedBufferFiles(const PEID rank, const PGeneratorConfig& config, const bool output_info) {
    if (rank != ROOT) {
        return;
    }

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

std::pair<PGeneratorConfig, GraphInfo>
GenerateChunks(PGeneratorConfig config, const std::vector<SInt>& vertex_distribution, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const bool output_error = (rank == ROOT);
    const bool output_info  = (rank == ROOT && !config.quiet);

    GraphInfo local_info;
    auto      generator_factory = CreateGeneratorFactory(config.generator);

    std::ofstream out_index(BufferFilename(rank, config, true), std::ios_base::binary | std::ios_base::trunc);
    std::ofstream out_edges(BufferFilename(rank, config, false), std::ios_base::binary | std::ios_base::trunc);
    if (!out_index || !out_edges) {
        if (output_error) {
            std::cerr << "Error: cannot create files in " << config.external.tmp_directory << "\n";
        }
        MPI_Abort(comm, 1);
    }

    SInt cur_edge_offset = 0;

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
        local_info.local_n += graph.NumberOfLocalVertices();
        local_info.local_m += graph.NumberOfLocalEdges();
        local_info.has_vertex_weights |= !graph.vertex_weights.empty();

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

        if (output_info) {
            std::cout << "sorting edges ... " << std::flush;
        }
        std::sort(edges.begin(), edges.end());

        if (output_info) {
            std::cout << "writing to external buffer ... " << std::flush;
        }
        cur_edge_offset = SwapoutGraphChunk(out_index, out_edges, cur_edge_offset, edges, vertex_distribution, config);

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

    return {config, local_info};
}

GraphInfo DeduplicateEdges(
    GraphInfo info, PGeneratorConfig config, const Indices& my_indices, const std::vector<SInt>& vertex_distribution,
    MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const bool output_info = (rank == ROOT && !config.quiet);

    info.local_m = 0;
    info.local_n = 0;

    const std::string              agg_filename = AggregatedBufferFilename(rank, config);
    std::unique_ptr<std::ofstream> agg_out =
        config.external.cache_aggregated_chunks
            ? std::make_unique<std::ofstream>(agg_filename, std::ios_base::binary | std::ios_base::trunc)
            : nullptr;

    for (int chunk = rank; chunk < config.external.num_chunks; chunk += size) {
        if (output_info) {
            std::cout << "Counting edges (chunk " << chunk + 1 << "... / " << config.external.num_chunks << ") ... "
                      << std::flush;
        }

        Graph graph = SwapinGraphChunk(chunk, my_indices, vertex_distribution, nullptr, agg_out.get(), config, comm);
        info.local_n += graph.NumberOfLocalVertices();
        info.local_m += graph.NumberOfLocalEdges();

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

    return info;
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

    auto vertex_distribution = CreateVertexDistribution(config.n, config.external.num_chunks);

    GraphInfo local_info;
    std::tie(config, local_info) = GenerateChunks(config, vertex_distribution, comm);

    const auto my_indices = ReadMyIndices(rank, size, config);

    if (config.external.fix_reverse_edges || config.external.fix_nonlocal_reverse_edges) {
        local_info = DeduplicateEdges(local_info, config, my_indices, vertex_distribution, comm);
        RemoveBufferFiles(rank, config, output_info);
    }

    const GraphInfo global_info(local_info, comm);

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

            // Read aggregated buffers for each pass from beginning
            const std::string              agg_filename = AggregatedBufferFilename(rank, config);
            std::unique_ptr<std::ifstream> agg_in       = [&] {
                if (std::ifstream in(agg_filename, std::ios_base::binary); in) {
                    return std::make_unique<std::ifstream>(std::move(in));
                }
                return std::unique_ptr<std::ifstream>(nullptr);
            }();

            // In contrast to the other loops, we have to do some communication during the "write to disk" loop
            // Thus, we have to make sure that all PEs participate in this loop: we achieve this by rounding the number
            // of chunks to the next multiple of the number of PEs and make sure that the "dummy rounds" on some PEs are
            // a no-op in terms of IO
            const PEID rounded_chunk_count = std::ceil(1.0 * config.external.num_chunks / size) * size;

            for (PEID chunk = rank; chunk < rounded_chunk_count; chunk += size) {
                if (output_info) {
                    std::cout << "Writing " << out_config.filename << " (pass " << pass + 1 << ", chunk " << chunk + 1
                              << "... / " << config.external.num_chunks << ") ... " << std::flush;
                }

                // @todo to determine whether we want to cache the aggregated buffers, we would have to know whether we
                // need multiple IO passes or not -- extend IO interface to give this information?
                Graph graph = [&] {
                    if (chunk < config.external.num_chunks) {
                        return SwapinGraphChunk(
                            chunk, my_indices, vertex_distribution, agg_in.get(), nullptr, config, comm);
                    }

                    return Graph{};
                }();

                // global_info contains information about the graph distribution on a "per PE" level, but we need this
                // information on a "per chunk" level
                // We have to compute this information by our self:
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
                    std::cout << "writing(m=" << this_round_m << "...) ... " << std::flush;
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
    }

    RemoveBufferFiles(rank, config, output_info);
    RemoveAggregatedBufferFiles(rank, config, output_info);
}
} // namespace kagen

