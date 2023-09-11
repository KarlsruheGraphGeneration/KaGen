#include "app/CLI11.h"
#include "app/tools/strutils.h"

#include "kagen/context.h"
#include "kagen/generators/file/file_graph.h"
#include "kagen/io.h"
#include "kagen/kagen.h"
#include "kagen/tools/statistics.h"
#include "kagen/tools/utils.h"

#include <mpi.h>

#include <iostream>

using namespace kagen;

struct Statistics {
    std::string name;
    SInt        n;
    SInt        m;
    SInt        min_deg;
    LPFloat     avg_deg;
    SInt        max_deg;
};

Graph LoadGraph(const PGeneratorConfig& config) {
    const PEID rank = GetCommRank(MPI_COMM_WORLD);
    const PEID size = GetCommSize(MPI_COMM_WORLD);

    FileGraphFactory factory;
    const auto       normalized_config = factory.NormalizeParameters(config, rank, size, false);
    auto             loader            = factory.Create(normalized_config, rank, size);
    loader->Generate(GraphRepresentation::EDGE_LIST);
    loader->Finalize(MPI_COMM_WORLD);
    return loader->Take();
}

void PrintHeader() {
    std::cout << "Graph,";
    std::cout << "N,";
    std::cout << "M,";
    std::cout << "MinDeg,";
    std::cout << "AvgDeg,";
    std::cout << "MaxDeg";
    std::cout << std::endl;
}

void PrintRow(const Statistics& stats) {
    std::cout << stats.name << ",";
    std::cout << stats.n << ",";
    std::cout << stats.m << ",";
    std::cout << stats.min_deg << ",";
    std::cout << stats.avg_deg << ",";
    std::cout << stats.max_deg;
    std::cout << std::endl;
}

Statistics GenerateInternal(const PGeneratorConfig& config) {
    Graph graph = LoadGraph(config);

    Statistics stats;
    stats.n = FindNumberOfGlobalNodes(graph.vertex_range, MPI_COMM_WORLD);
    stats.m = FindNumberOfGlobalEdges(graph.edges, MPI_COMM_WORLD);

    const auto degree_stats = ReduceDegreeStatistics(graph.edges, stats.n, MPI_COMM_WORLD);
    stats.min_deg           = degree_stats.min;
    stats.avg_deg           = degree_stats.mean;
    stats.max_deg           = degree_stats.max;

    return stats;
}

Statistics GenerateExternal(const PGeneratorConfig& config, const int num_chunks) {
    if (GetCommSize(MPI_COMM_WORLD) > 1) {
        std::cerr << "Error: external statistics generation is only supported for a single MPI process\n";
        std::exit(1);
    }

    Statistics stats;

    const auto reader        = CreateGraphReader(config.input_graph.format, config.input_graph, 0, 1);
    auto       reported_size = reader->ReadSize();

    std::vector<SInt> degrees;

    for (int chunk = 0; chunk < num_chunks; ++chunk) {
        const auto [from, to] = ComputeRange(reported_size.first, num_chunks, chunk);
        Graph graph = reader->Read(from, to, std::numeric_limits<SInt>::max(), GraphRepresentation::EDGE_LIST);

        for (const auto& [from, to]: graph.edges) {
            while (degrees.size() <= from) {
                degrees.push_back(0);
            }
            ++degrees[from];
        }

        stats.m += graph.edges.size();
    }

    const auto [min_it, max_it] = std::minmax_element(degrees.begin(), degrees.end());

    stats.n       = degrees.size();
    stats.min_deg = *min_it;
    stats.avg_deg = 1.0 * stats.m / stats.n;
    stats.max_deg = *max_it;

    return stats;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    std::vector<std::string> input_filenames;
    bool                     do_strip_extension = false;
    bool                     do_no_header       = false;
    bool                     do_header_only     = false;
    int                      num_chunks         = 1;
    PGeneratorConfig         config;

    CLI::App app("graphstats: compute some basic statistics on a graph");

    CLI::Option_group* group = app.add_option_group("Options");
    group->require_option(1);
    group->add_option("input filenames", input_filenames)->check(CLI::ExistingFile);
    group->add_flag("--header-only", do_header_only);

    app.add_option("-f,--format", config.input_graph.format, "File format of the input file(s).")
        ->transform(CLI::CheckedTransformer(GetInputFormatMap()));
    app.add_flag(
        "--strip-extension", do_strip_extension,
        "If set, print the filename in the Graph column without file extension.");
    app.add_flag("-H,--no-header", do_no_header, "If set, do not print the CSV header line.");
    app.add_option(
        "-C,--num-chunks", num_chunks,
        "If set, compute the statistics externally by splitting the graph into this many chunks; some statistics might "
        "not be available in this mode. Still requires O(n) memory.");
    CLI11_PARSE(app, argc, argv);

    // Catch special case: only print CSV header line
    if ((do_header_only || !do_no_header) && GetCommRank(MPI_COMM_WORLD) == ROOT) {
        PrintHeader();
    }
    if (do_header_only) {
        return MPI_Finalize();
    }

    for (const auto& filename: input_filenames) {
        config.input_graph.filename = filename;

        Statistics stats;
        if (num_chunks == 1) {
            stats = GenerateInternal(config);
        } else {
            stats = GenerateExternal(config, num_chunks);
        }

        stats.name = ExtractFilename(config.input_graph.filename);
        if (do_strip_extension) {
            stats.name = StripExtension(stats.name);
        }

        if (GetCommRank(MPI_COMM_WORLD) == ROOT) {
            PrintRow(stats);
        }
    }

    return MPI_Finalize();
}
