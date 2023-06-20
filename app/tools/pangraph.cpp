#include <mpi.h>

#include "../CLI11.h"
#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/io.h"
#include "kagen/kagen.h"
#include "kagen/tools/postprocessor.h"
#include "kagen/tools/utils.h"

#include <algorithm>
#include <filesystem>

using namespace kagen;

std::string CreateGraphChunkFilename(const std::string& tmp_directory, const PEID fake_rank) {
    return tmp_directory + "/" + std::to_string(fake_rank) + ".buf";
}

void WriteGraphChunk(const std::string& filename, const Graph& graph) {
    std::ofstream out(filename, std::ios::binary | std::ios::trunc);

    const SInt from = graph.vertex_range.first;
    const SInt to   = graph.vertex_range.second;
    out.write(reinterpret_cast<const char*>(&from), sizeof(from));
    out.write(reinterpret_cast<const char*>(&to), sizeof(to));

    const std::size_t edges_size          = graph.edges.size();
    const std::size_t vertex_weights_size = graph.vertex_weights.size();
    const std::size_t edge_weights_size   = graph.edge_weights.size();
    out.write(reinterpret_cast<const char*>(&edges_size), sizeof(edges_size));
    out.write(reinterpret_cast<const char*>(&vertex_weights_size), sizeof(vertex_weights_size));
    out.write(reinterpret_cast<const char*>(&edge_weights_size), sizeof(edge_weights_size));

    out.write(reinterpret_cast<const char*>(graph.edges.data()), sizeof(typename EdgeList::value_type) * edges_size);
    out.write(
        reinterpret_cast<const char*>(graph.vertex_weights.data()),
        sizeof(typename VertexWeights::value_type) * vertex_weights_size);
    out.write(
        reinterpret_cast<const char*>(graph.edge_weights.data()),
        sizeof(typename EdgeWeights::value_type) * edge_weights_size);
}

Graph ReadGraphChunk(const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);

    auto read_sint = [&in] {
        SInt size;
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        return size;
    };

    Graph graph;
    graph.vertex_range.first  = read_sint();
    graph.vertex_range.second = read_sint();
    graph.edges.resize(read_sint());
    graph.vertex_weights.resize(read_sint());
    graph.edge_weights.resize(read_sint());

    in.read(reinterpret_cast<char*>(graph.edges.data()), sizeof(typename EdgeList::value_type) * graph.edges.size());
    in.read(
        reinterpret_cast<char*>(graph.vertex_weights.data()),
        sizeof(typename VertexWeights::value_type) * graph.vertex_weights.size());
    in.read(
        reinterpret_cast<char*>(graph.edge_weights.data()),
        sizeof(typename EdgeWeights::value_type) * graph.edge_weights.size());
    return graph;
}

void RemoveSelfLoops(Graph& graph) {
    if (!graph.edge_weights.empty()) {
        throw std::runtime_error("not implemented");
    }

    auto it = std::remove_if(graph.edges.begin(), graph.edges.end(), [](const auto& e) { return e.first == e.second; });
    graph.edges.erase(it, graph.edges.end());
}

void RemoveMultiEdges(Graph& graph) {
    if (!graph.edge_weights.empty()) {
        throw std::runtime_error("not implemented");
    }

    graph.SortEdgelist();
    auto it = std::unique(graph.edges.begin(), graph.edges.end());
    graph.edges.erase(it, graph.edges.end());
}

void ProcessEdges(Graph& graph, const bool remove_self_loops, const bool add_reverse_edges) {
    if (remove_self_loops) {
        RemoveSelfLoops(graph);
    }
    if (add_reverse_edges) {
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
}

PEID FindOwner(const SInt node, const std::vector<SInt>& vertex_distribution) {
    auto it = std::upper_bound(vertex_distribution.begin(), vertex_distribution.end(), node);
    return static_cast<PEID>(std::distance(vertex_distribution.begin(), it)) - 1;
}

std::string SendbufFilename(const std::string& tmp_directory, const PEID fake_rank, const PEID fake_size) {
    return tmp_directory + "/sendbuf_from" + std::to_string(fake_rank) + "_to" + std::to_string(fake_size) + ".edges";
}

void RedistributeToExternalMemory(
    const Graph& graph, const std::vector<SInt>& vertex_distribution, const PEID fake_rank, const PEID fake_size,
    const std::string& tmp_directory) {
    if (!graph.edge_weights.empty()) {
        throw std::runtime_error("not implemented");
    }

    std::vector<EdgeList> sendbufs(fake_size);
    for (const auto& [from, to]: graph.edges) {
        sendbufs[FindOwner(from, vertex_distribution)].emplace_back(from, to);
        sendbufs[FindOwner(to, vertex_distribution)].emplace_back(to, from);
    }

    for (PEID fake_pe = 0; fake_pe < fake_size; ++fake_pe) {
        std::ofstream out(SendbufFilename(tmp_directory, fake_rank, fake_pe), std::ios::binary);

        const SInt size = sendbufs[fake_pe].size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));
        out.write(
            reinterpret_cast<const char*>(sendbufs[fake_pe].data()), sizeof(typename EdgeList::value_type) * size);
    }
}

Graph RestoreFromExternalMemory(
    const std::vector<SInt>& vertex_distribution, const PEID fake_rank, const PEID fake_size,
    const std::string& tmp_directory) {
    SInt total_size = 0;

    auto read_sint = [](std::ifstream& in) {
        SInt size;
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        return size;
    };

    for (PEID fake_pe = 0; fake_pe < fake_size; ++fake_pe) {
        std::ifstream in(SendbufFilename(tmp_directory, fake_pe, fake_rank), std::ios::binary);
        total_size += read_sint(in);
    }

    Graph graph;
    graph.vertex_range.first  = vertex_distribution[fake_rank];
    graph.vertex_range.second = vertex_distribution[fake_rank + 1];
    graph.edges.resize(total_size);

    SInt position = 0;
    for (PEID fake_pe = 0; fake_pe < fake_size; ++fake_pe) {
        std::ifstream in(SendbufFilename(tmp_directory, fake_pe, fake_rank), std::ios::binary);
        const SInt    size = read_sint(in);
        in.read(reinterpret_cast<char*>(graph.edges.data() + position), sizeof(typename EdgeList::value_type) * size);
    }

    return graph;
}

void WriteExternallyRedistributedGraph(
    const FileFormatFactory& factory, const OutputGraphConfig& config, const GraphInfo info,
    const std::vector<SInt>& vertex_distribution, const int num_chunks, const bool add_reverse_edges,
    const std::string& tmp_directory) {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    const PEID fake_size = size * num_chunks;

    bool continue_with_next_pass = true;
    for (int pass = 0; continue_with_next_pass; ++pass) {
        for (PEID pe = 0; pe < size; ++pe) {
            if (pe == rank) {
                for (int chunk = 0; chunk < num_chunks; ++chunk) {
                    const PEID fake_rank = rank * num_chunks + chunk;

                    Graph graph = RestoreFromExternalMemory(vertex_distribution, fake_rank, fake_size, tmp_directory);
                    ProcessEdges(graph, false, add_reverse_edges);

                    auto writer             = factory.CreateWriter(config, graph, info, fake_rank, fake_size);
                    continue_with_next_pass = writer->Write(pass, config.filename);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}

void WriteExternalGraph(
    const FileFormatFactory& factory, const OutputGraphConfig& config, const GraphInfo info,
    const std::string& tmp_directory) {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    bool continue_with_next_pass = true;
    for (int pass = 0; continue_with_next_pass; ++pass) {
        for (PEID pe = 0; pe < size; ++pe) {
            if (pe == rank) {
                Graph graph = ReadGraphChunk(CreateGraphChunkFilename(tmp_directory, rank));

                auto writer             = factory.CreateWriter(config, graph, info, rank, size);
                continue_with_next_pass = writer->Write(pass, config.filename);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}

int main(int argc, char* argv[]) {
    OutputGraphConfig out_config;
    InputGraphConfig  in_config;
    in_config.width = 64;

    // General options
    bool quiet = false;

    // External memory options
    int         num_chunks    = 1;
    std::string tmp_directory = std::filesystem::temp_directory_path();

    // Transformation options
    bool remove_self_loops = false;
    bool add_reverse_edges = false;

    CLI::App app("pangraph: distributed and/or external graph format converter");
    app.add_option("-C,--chunks", num_chunks)->capture_default_str();
    app.add_option("-T,--tmp-directory", tmp_directory, "Directory for external memory buffers.")
        ->capture_default_str();

    app.add_option("--input-filename", in_config.filename, "Input graph")->check(CLI::ExistingFile)->required();
    app.add_option("--input-format", in_config.format, "Input graph format")
        ->transform(CLI::CheckedTransformer(GetInputFormatMap()))
        ->required()
        ->capture_default_str();
    app.add_option("--output-format", out_config.formats, "Output graph format")
        ->transform(CLI::CheckedTransformer(GetOutputFormatMap()))
        ->required()
        ->capture_default_str();
    app.add_option("--output-filename", out_config.filename, "Output graph")->required();
    app.add_flag("-q,--quiet", quiet, "Suppress any output to stdout.");

    app.add_flag("--remove-self-loops", remove_self_loops, "Remove self loops from the input graph.")
        ->capture_default_str();
    app.add_flag(
           "--add-reverse-edges", add_reverse_edges,
           "Add reverse edges to the input graph, such that the output graph is undirected.")
        ->capture_default_str();
    CLI11_PARSE(app, argc, argv);

    MPI_Init(&argc, &argv);
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    const PEID fake_size = num_chunks * size;

    const bool output                  = !quiet && rank == ROOT;
    bool       requires_redistribution = false;

    Graph     graph;
    GraphInfo info;

    for (int chunk = 0; chunk < num_chunks; ++chunk) {
        if (output) {
            std::cout << "Reading " << in_config.filename;
            if (num_chunks > 1) {
                std::cout << " (chunk " << chunk << " of " << num_chunks << ")";
            }
            std::cout << " ..." << std::endl;
        }

        const PEID fake_rank  = rank * num_chunks + chunk;
        const auto reader     = CreateGraphReader(in_config.format, in_config, fake_rank, fake_size);
        const auto [n, m]     = reader->ReadSize();
        const auto [from, to] = ComputeRange(n, fake_size, fake_rank);
        graph = reader->Read(from, to, std::numeric_limits<SInt>::max(), GraphRepresentation::EDGE_LIST);
        requires_redistribution |= (reader->Deficits() & ReaderDeficits::REQUIRES_REDISTRIBUTION);

        ProcessEdges(graph, remove_self_loops, add_reverse_edges);

        for (const auto& [from, to]: graph.edges) {
            info.global_n = std::max(info.global_n, std::max(from, to) + 1);
        }
        info.global_n += graph.NumberOfLocalVertices();
        info.global_m += graph.NumberOfLocalEdges();
        info.has_vertex_weights |= !graph.vertex_weights.empty();
        info.has_edge_weights |= !graph.edge_weights.empty();

        if (num_chunks > 1) {
            WriteGraphChunk(CreateGraphChunkFilename(tmp_directory, fake_rank), graph);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Allreduce(MPI_IN_PLACE, &info.global_n, 1, KAGEN_MPI_SINT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &info.global_m, 1, KAGEN_MPI_SINT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &info.has_vertex_weights, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &info.has_edge_weights, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

    std::vector<SInt> vertex_distribution(fake_size + 1);
    for (int fake_rank = 0; fake_rank < fake_size; ++fake_rank) {
        vertex_distribution[fake_rank + 1] = ComputeRange(info.global_n, fake_size, fake_rank).second;
    }

    // @todo optimize the second case
    if (requires_redistribution || (num_chunks > 1 && add_reverse_edges)) {
        if (num_chunks == 1) {
            if (output) {
                std::cout << "Redistributing in internal memory ... " << std::flush;
            }
            AddReverseEdgesAndRedistribute(graph.edges, graph.vertex_range, add_reverse_edges, MPI_COMM_WORLD);
            if (output) {
                std::cout << "OK" << std::endl;
            }
        } else {
            if (output) {
                std::cout << "Redistributing in external memory ... " << std::flush;
            }
            for (int chunk = 0; chunk < num_chunks; ++chunk) {
                const PEID fake_rank = rank * num_chunks + chunk;

                graph = ReadGraphChunk(CreateGraphChunkFilename(tmp_directory, fake_rank));
                RedistributeToExternalMemory(graph, vertex_distribution, fake_rank, fake_size, tmp_directory);
            }
            if (output) {
                std::cout << "OK" << std::endl;
            }
        }
    }

    const std::string base_filename = out_config.filename;
    for (const FileFormat& format: out_config.formats) {
        const auto& factory = GetGraphFormatFactory(format);

        // If there are multiple output formats, append the default extension of the each file format to avoid
        // conflicts
        if (out_config.formats.size() > 1) {
            out_config.filename = base_filename + "." + factory->DefaultExtension();
        }

        if (num_chunks == 1) {
            auto writer = factory->CreateWriter(out_config, graph, info, rank, size);
            WriteGraph(*writer, out_config, output, MPI_COMM_WORLD);
        } else if (requires_redistribution || (num_chunks > 1 && add_reverse_edges)) {
            WriteExternallyRedistributedGraph(
                *factory, out_config, info, vertex_distribution, num_chunks, add_reverse_edges, tmp_directory);
        } else {
            WriteExternalGraph(*factory, out_config, info, tmp_directory);
        }
    }

    return MPI_Finalize();
}
