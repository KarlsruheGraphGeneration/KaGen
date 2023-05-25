#include <mpi.h>

#include "kagen/context.h"
#include "kagen/generators/file/file_graph.h"
#include "kagen/kagen.h"
#include "kagen/tools/statistics.h"

#include "../CLI11.h"

#include <iostream>

using namespace kagen;

PEID get_rank() {
    PEID rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

PEID get_size() {
    PEID size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

Graph load_graph(const PGeneratorConfig& config) {
    const PEID rank = get_rank();
    const PEID size = get_size();

    FileGraphFactory factory;
    const auto       normalized_config = factory.NormalizeParameters(config, rank, size, false);
    auto             loader            = factory.Create(normalized_config, rank, size);

    loader->Generate(GraphRepresentation::EDGE_LIST);
    loader->Finalize(MPI_COMM_WORLD);
    return loader->Take();
}

struct Statistics {
    std::string name;
    SInt        n;
    SInt        m;
    SInt        min_deg;
    LPFloat     avg_deg;
    SInt        max_deg;
};

void print_csv_header() {
    std::cout << "Graph,";
    std::cout << "N,";
    std::cout << "M,";
    std::cout << "MinDeg,";
    std::cout << "AvgDeg,";
    std::cout << "MaxDeg";
    std::cout << std::endl;
}

void print_csv_row(const Statistics& stats) {
    std::cout << stats.name << ",";
    std::cout << stats.n << ",";
    std::cout << stats.m << ",";
    std::cout << stats.min_deg << ",";
    std::cout << stats.avg_deg << ",";
    std::cout << stats.max_deg;
    std::cout << std::endl;
}

std::string extract_filename(const std::string& filename) {
    const auto pos = filename.find_last_of('/');
    if (pos == std::string::npos) {
        return filename;
    }
    return filename.substr(pos + 1);
}

std::string strip_extension(const std::string& filename) {
    const auto pos = filename.find_last_of('.');
    if (pos == std::string::npos) {
        return filename;
    }
    return filename.substr(0, pos);
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    std::vector<std::string> input_filenames;
    bool                     do_strip_extension = false;
    bool                     do_no_header       = false;
    bool                     do_header_only     = false;
    PGeneratorConfig         config;

    CLI::App app("graphstats: compute some basic statistics on a graph");

    CLI::Option_group* group = app.add_option_group("Options");
    group->require_option(1);
    group->add_option("input filenames", input_filenames)->check(CLI::ExistingFile);
    group->add_flag("--header-only", do_header_only);

    app.add_option("-f,--format", config.input_graph.format)->transform(CLI::CheckedTransformer(GetInputFormatMap()));
    app.add_flag("--strip-extension", do_strip_extension);
    app.add_flag("-H,--no-header", do_no_header);
    CLI11_PARSE(app, argc, argv);

    // Catch special case: only print CSV header line
    if ((do_header_only || !do_no_header) && get_rank() == ROOT) {
        print_csv_header();
    }
    if (do_header_only) {
        return MPI_Finalize();
    }

    for (const auto& filename: input_filenames) {
        config.input_graph.filename = filename;
        Graph graph                 = load_graph(config);

        Statistics stats;
        stats.name = extract_filename(config.input_graph.filename);
        if (do_strip_extension) {
            stats.name = strip_extension(stats.name);
        }
        stats.n = FindNumberOfGlobalNodes(graph.vertex_range, MPI_COMM_WORLD);
        stats.m = FindNumberOfGlobalEdges(graph.edges, MPI_COMM_WORLD);

        const auto degree_stats = ReduceDegreeStatistics(graph.edges, stats.n, MPI_COMM_WORLD);
        stats.min_deg           = degree_stats.min;
        stats.avg_deg           = degree_stats.mean;
        stats.max_deg           = degree_stats.max;

        if (get_rank() == ROOT) {
            print_csv_row(stats);
        }
    }

    return MPI_Finalize();
}
