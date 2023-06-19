#include <mpi.h>

#include "../CLI11.h"
#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/io.h"
#include "kagen/kagen.h"
#include "kagen/tools/utils.h"

#include <filesystem>

using namespace kagen;

std::string create_tmp_chunk_filename(const std::string& tmp_directory, const PEID chunk) {
    return tmp_directory + "/" + std::to_string(chunk) + ".edges";
}

void dump_edgelist(const std::string& filename, const EdgeList& edges) {
    std::ofstream     out(filename, std::ios::binary | std::ios::trunc);
    const std::size_t size = edges.size();
    out.write(reinterpret_cast<const char*>(&size), sizeof(size));
    out.write(reinterpret_cast<const char*>(edges.data()), sizeof(typename EdgeList::value_type) * size);
}

EdgeList read_edgelist(const std::string& filename) {
    std::ifstream     in(filename, std::ios::binary);
    const std::size_t size = [&in]() {
        std::size_t size;
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        return size;
    }();
    EdgeList edges(size);
    in.read(reinterpret_cast<char*>(edges.data()), sizeof(typename EdgeList::value_type) * size);
    return edges;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    OutputGraphConfig output;
    InputGraphConfig  input;
    input.width = 64;

    // General options
    bool quiet = false;

    // External memory options
    int         num_chunks    = 1;
    std::string tmp_directory = std::filesystem::temp_directory_path();

    // Transformation options
    bool remove_self_loops  = false;
    bool remove_multi_edges = false;
    bool add_reverse_edges  = false;

    CLI::App app("pangraph: distributed and/or external graph format converter");
    app.add_option("-C,--chunks", num_chunks)->capture_default_str();
    app.add_option("-T,--tmp-directory", tmp_directory, "Directory for external memory buffers.")
        ->capture_default_str();

    app.add_option("--input-filename", input.filename, "Input graph")->check(CLI::ExistingFile)->required();
    app.add_option("--input-format", input.format, "Input graph format")
        ->transform(CLI::CheckedTransformer(GetInputFormatMap()))
        ->required()
        ->capture_default_str();
    app.add_option("--output-format", output.formats, "Output graph format")
        ->transform(CLI::CheckedTransformer(GetOutputFormatMap()))
        ->required()
        ->capture_default_str();
    app.add_option("--output-filename", output.filename, "Output graph")->required();
    app.add_flag("-q,--quiet", quiet, "Suppress any output to stdout.");

    app.add_flag("--remove-self-loops", remove_self_loops, "Remove self loops from the input graph.")
        ->capture_default_str();
    app.add_flag("--remove-multi-edges", remove_multi_edges, "Remove multi edges from the input graph.");
    app.add_flag(
           "--add-reverse-edges", add_reverse_edges,
           "Add reverse edges to the input graph, such that the output graph is undirected.")
        ->capture_default_str();
    CLI11_PARSE(app, argc, argv);

    Graph graph;
    for (int chunk = 0; chunk < num_chunks; ++chunk) {
        const PEID fake_size = size * num_chunks;
        const PEID fake_rank = rank * num_chunks + chunk;

        const auto reader     = CreateGraphReader(input.format, input, fake_rank, fake_size);
        const auto [n, m]     = reader->ReadSize();
        const auto [from, to] = ComputeRange(n, fake_size, fake_rank);
        graph = reader->Read(from, to, std::numeric_limits<SInt>::max(), GraphRepresentation::EDGE_LIST);

        if (num_chunks > 1) {
            dump_edgelist(create_tmp_chunk_filename(tmp_directory, chunk), graph.edges);
        }
    }

    const std::string base_filename = output.filename;
    for (const FileFormat& format: output.formats) {
        const auto& factory = GetGraphFormatFactory(format);

        // If there are multiple output formats, append the default extension of the each file format to avoid conflicts
        if (output.formats.size() > 1) {
            output.filename = base_filename + "." + factory->DefaultExtension();
        }

        auto writer = factory->CreateWriter(output, graph, MPI_COMM_WORLD);
        if (writer != nullptr) {
            writer->Write(!quiet);
        } else if (!quiet && rank == 0) {
            std::cout << "Warning: invalid file format " << format << " for writing; skipping\n";
        }
    }

    return MPI_Finalize();
}
