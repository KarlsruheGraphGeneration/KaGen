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

struct Configuration {
    std::vector<std::string> input_filenames;
    FileFormat               input_format = FileFormat::EXTENSION;

    int  num_chunks      = 1;
    bool header_only     = false;
    bool omit_header     = false;
    bool strip_extension = false;

    std::vector<SInt> count_num_deg_nodes;
};

struct Statistics {
    std::string name;

    SInt    n       = 0;
    SInt    m       = 0;
    SInt    min_deg = 0;
    LPFloat avg_deg = 0.0;
    SInt    max_deg = 0;

    std::vector<SInt> num_deg_nodes;
};

void PrintHeader(const Configuration& config) {
    ((void)config);

    std::cout << "Graph,";
    std::cout << "N,";
    std::cout << "M,";
    for (const SInt deg: config.count_num_deg_nodes) {
        std::cout << "NumDeg" << deg << "Nodes,";
    }
    std::cout << "MinDeg,";
    std::cout << "AvgDeg,";
    std::cout << "MaxDeg";
    std::cout << std::endl;
}

void PrintRow(const Statistics& stats) {
    std::cout << stats.name << ",";
    std::cout << stats.n << ",";
    std::cout << stats.m << ",";
    for (const SInt num_deg_nodes: stats.num_deg_nodes) {
        std::cout << num_deg_nodes << ",";
    }
    std::cout << stats.min_deg << ",";
    std::cout << stats.avg_deg << ",";
    std::cout << stats.max_deg;
    std::cout << std::endl;
}

struct StatisticsComputator {
    StatisticsComputator(const Configuration& config) : config_(config) {}

    void operator()(const GraphFragment& fragment) {
        if (degrees_.size() < fragment.graph.vertex_range.second) {
            degrees_.resize(fragment.graph.vertex_range.second, 0);
        }

        const auto& graph = fragment.graph;

        stats_.m += graph.edges.size();

        for (const auto& [from, to]: graph.edges) {
            while (degrees_.size() <= from) {
                degrees_.push_back(0);
            }
            ++degrees_[from];
        }
    }

    Statistics Finalize(const Graph& graph) {
        FinalizeStreamingStatistics();
        return std::move(stats_);
    }

    Statistics Finalize() {
        FinalizeStreamingStatistics();
        return std::move(stats_);
    }

private:
    void FinalizeStreamingStatistics() {
        stats_.n = degrees_.size();
        FinalizeDegreeStatistics();
    }

    void FinalizeDegreeStatistics() {
        stats_.min_deg = std::numeric_limits<SInt>::max();
        stats_.max_deg = std::numeric_limits<SInt>::min();

        stats_.num_deg_nodes.resize(config_.count_num_deg_nodes.size());
        std::fill(stats_.num_deg_nodes.begin(), stats_.num_deg_nodes.end(), 0);

        for (SInt node = 0; node < stats_.n; ++node) {
            const SInt deg = degrees_[node];
            stats_.min_deg = std::min(stats_.min_deg, deg);
            stats_.max_deg = std::max(stats_.max_deg, deg);

            for (std::size_t i = 0; i < config_.count_num_deg_nodes.size(); ++i) {
                stats_.num_deg_nodes[i] += (deg == config_.count_num_deg_nodes[i]);
            }
        }

        stats_.avg_deg = 1.0 * stats_.m / stats_.n;
    }

    const Configuration& config_;

    std::vector<SInt> degrees_;

    Statistics stats_;
};

Statistics ComputeStatistics(const Configuration& stats_config, const PGeneratorConfig& kagen_config) {
    StatisticsComputator computator(stats_config);

    auto reader = CreateGraphReader(kagen_config.input_graph.format, kagen_config.input_graph, 0, 1);

    GraphFragment first_fragment = ReadGraphFragment(
        *reader, GraphRepresentation::EDGE_LIST, kagen_config.input_graph, 0, stats_config.num_chunks);
    computator(first_fragment);

    for (int chunk = 1; chunk < stats_config.num_chunks; ++chunk) {
        const GraphFragment fragment = ReadGraphFragment(
            *reader, GraphRepresentation::EDGE_LIST, kagen_config.input_graph, chunk, stats_config.num_chunks);
        computator(fragment);
    }

    if (stats_config.num_chunks == 1) {
        const Graph graph = FinalizeGraphFragment(std::move(first_fragment), false, MPI_COMM_WORLD);
        return computator.Finalize(graph);
    } else {
        return computator.Finalize();
    }
}

Configuration parse_cli_arguments(int argc, char* argv[]) {
    Configuration config;

    CLI::App app("graphstats: compute basic graph statistics");

    CLI::Option_group* group = app.add_option_group("Options");
    group->require_option(1);
    group->add_option("input filenames", config.input_filenames)->check(CLI::ExistingFile);
    group->add_flag("--header-only", config.header_only);

    app.add_option(
        "-C,--num-chunks", config.num_chunks,
        "If set, compute the statistics externally by splitting the graph into this many chunks; some statistics might "
        "not be available in this mode. Still requires O(n) memory.");
    app.add_option("-f,--format", config.input_format, "File format of the input file(s).")
        ->transform(CLI::CheckedTransformer(GetInputFormatMap()));
    app.add_flag(
        "--strip-extension", config.strip_extension,
        "If set, print the filename in the Graph column without file extension.");
    app.add_flag("-H,--omit-header", config.omit_header, "If set, do not print the CSV header line.");

    app.add_option("--count-degree", config.count_num_deg_nodes, "Count the number of nodes with this degree.");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        (app).exit(e);
        std::exit(1);
    }

    return config;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    if (GetCommSize(MPI_COMM_WORLD) != 1) {
        std::cerr << "must be run with just one MPI process\n";
        return MPI_Finalize();
    }

    Configuration config = parse_cli_arguments(argc, argv);

    if (config.header_only || !config.omit_header) {
        PrintHeader(config);
    }
    if (config.header_only) {
        return MPI_Finalize();
    }

    for (const auto& filename: config.input_filenames) {
        PGeneratorConfig kagen_config;
        kagen_config.input_graph.filename = filename;
        kagen_config.input_graph.format   = config.input_format;

        Statistics stats = ComputeStatistics(config, kagen_config);
        stats.name       = ExtractFilename(filename);
        if (config.strip_extension) {
            stats.name = StripExtension(stats.name);
        }

        PrintRow(stats);
    }

    return MPI_Finalize();
}
