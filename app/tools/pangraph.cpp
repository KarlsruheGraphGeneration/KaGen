#include <mpi.h>

#include "../CLI11.h"
#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/kagen.h"

using namespace kagen;

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    OutputGraphConfig output;
    InputGraphConfig  input;
    input.width = 64;

    bool quiet              = false;
    bool remove_self_loops  = false;
    bool remove_multi_edges = false;
    bool add_reverse_edges  = false;

    CLI::App app("pangraph: distributed and/or external graph format converter");
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
}
