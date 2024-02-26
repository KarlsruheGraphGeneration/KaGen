/*******************************************************************************
 * app/generate_kagen.cpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Daniel Funke <funke@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/facade.h"
#include "kagen/io.h"

#include <mpi.h>

#include "CLI11.h"
#include <iostream>

using namespace kagen;

void PrintVersion() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == ROOT) {
        std::cout << BuildDescription() << std::endl;
    }
}

void SetupCommandLineArguments(CLI::App& app, PGeneratorConfig& config) {
    auto log_cb = [&](SInt& result) {
        return [&](const SInt value) {
            result = static_cast<SInt>(1) << value;
        };
    };

    auto add_option_n = [&](CLI::App* cmd) {
        auto* opt_log_n = cmd->add_option_function<SInt>("-N,--log-nodes", log_cb(config.n), "Logarithm value");
        auto* opt_n     = cmd->add_option("-n,--num-nodes", config.n, "Exact value");
        auto* group     = cmd->add_option_group("Number of nodes");
        group->add_options(opt_log_n, opt_n);
        group->require_option(1);
        group->silent();
        return group;
    };
    auto add_option_m = [&](CLI::App* cmd) {
        auto* opt_log_m = cmd->add_option_function<SInt>("-M,--log-edges", log_cb(config.m), "Logarithm value");
        auto* opt_m     = cmd->add_option("-m,--num-edges", config.m, "Exact value");
        auto* group     = cmd->add_option_group("Number of edges");
        group->add_options(opt_log_m, opt_m);
        group->require_option(1);
        group->silent();
        return group;
    };
    auto add_option_p = [&](CLI::App* cmd) {
        return cmd->add_option("-p,--prob", config.p, "Edge probability");
    };
    auto add_option_self_loops = [&](CLI::App* cmd) {
        return cmd->add_flag("--self-loops", config.self_loops, "Allow self loops");
    };
    auto add_option_directed = [&](CLI::App* cmd) {
        return cmd->add_flag("--directed", config.directed, "Generate a directed graph");
    };
    auto add_option_r = [&](CLI::App* cmd) {
        return cmd->add_option("-r,--radius", config.r, "Edge radius");
    };
    auto add_option_x = [&](CLI::App* cmd) {
        auto* opt_log_x = cmd->add_option_function<SInt>("-X,--grid-log-x", log_cb(config.grid_x), "Logarithm value");
        auto* opt_x     = cmd->add_option("-x,--grid-x", config.grid_x, "Exact value");
        auto* group     = cmd->add_option_group("Grid x dimension");
        group->add_options(opt_log_x, opt_x);
        group->require_option(1);
        group->silent();
        return group;
    };
    auto add_option_y = [&](CLI::App* cmd) {
        auto* opt_log_y = cmd->add_option_function<SInt>("-Y,--grid-log-y", log_cb(config.grid_y), "Logarithm value");
        auto* opt_y     = cmd->add_option("-y,--grid-y", config.grid_y, "Exact value");
        auto* group     = cmd->add_option_group("Grid y dimension");
        group->add_options(opt_log_y, opt_y);
        group->require_option(1);
        group->silent();
        return group;
    };
    auto add_option_z = [&](CLI::App* cmd) {
        auto* opt_log_z = cmd->add_option_function<SInt>("-Z,--grid-log-z", log_cb(config.grid_z), "Logarithm value");
        auto* opt_z     = cmd->add_option("-z,--grid-z", config.grid_z, "Exact value");
        auto* group     = cmd->add_option_group("Grid z dimension");
        group->add_options(opt_log_z, opt_z);
        group->require_option(1);
        group->silent();
        return group;
    };
    auto add_option_min_deg = [&](CLI::App* cmd) {
        return cmd->add_option("-d,--min-deg", config.min_degree, "Minimal vertex degree");
    };
    auto add_option_gamma = [&](CLI::App* cmd) {
        return cmd->add_option("-g,--gamma", config.plexp, "Power-law exponent");
    };
    auto add_option_avg_deg = [&](CLI::App* cmd) {
        return cmd->add_option("-d,--avg-deg", config.avg_degree, "Average vertex degree");
    };

    // Use 40 characters width for help
    auto formatter = std::make_shared<CLI::Formatter>();
    formatter->column_width(40);
    app.formatter(formatter);

    // Enable help-all
    app.set_help_all_flag("--help-all", "Show all help options");

    // Require exactly one subcommand == graph generator
    app.require_subcommand(1, 1);
    app.ignore_case();
    app.ignore_underscore();

    // Allow toplevel options to occur after the generator name
    app.fallthrough();

    // General parameters
    app.add_flag("-q,--quiet", config.quiet, "Quiet mode");
    app.add_flag(
           "-v,--version", [&](auto) { PrintVersion(); }, "Print KaGen version")
        ->trigger_on_parse();
    app.add_flag("-V,--validate", config.validate_simple_graph)
        ->description(
            R"(Validate that the generated graph is undirected and does not have duplicated edges or inconsistent edge weights.
This is mostly useful for experimental graph generators or when using KaGen to load graphs from disk.)");
    app.add_flag(
        "--skip-postprocessing", config.skip_postprocessing,
        "Skip postprocessing (repair inconsistency due to floating point inaccuracies etc.)");
    app.add_option("-s,--seed", config.seed, "Seed for PRNG (must be the same on all PEs)");
    auto* stats_group = app.add_option_group("Statistics output");
    stats_group->add_option("--stats", config.statistics_level)
        ->transform(CLI::CheckedTransformer(GetStatisticsLevelMap()).description(""))
        ->description(
            R"(Controls the level of statistics that are computed for the generated graph. Possible levels are:
  - none:     do not output any statistics
  - basic:    only output very basic statistics
  - advanced: also output information about the degree distribution)");
    stats_group->add_flag(
        "-S",
        [&config](const auto count) {
            if (count == 1) {
                config.statistics_level = StatisticsLevel::BASIC;
            } else if (count > 1) {
                config.statistics_level = StatisticsLevel::ADVANCED;
            }
        },
        "Controls how much statistics on the generated graph gets calculated; pass flag multiple times to increase "
        "extend");
    stats_group->silent();

    // Generator parameters
    app.add_option("-k,--num-chunks", config.k, "Number of chunks used for graph generation");
    app.add_option(
        "--automatic-num-chunks-imbalance-threshold", config.max_vertex_imbalance,
        "Controls the trade-off between vertex imbalance and number of chunks when deducing the number of chunks "
        "automatically");
    app.add_flag("-C,--coordinates", config.coordinates, "Generate coordinates (geometric generators only)");

    { // Options string
        auto* cmd = app.add_subcommand(
            "options",
            "Generate graph as specified by an options string; see library documentation for further details");
        cmd->add_option_function<std::string>(
               "options", [&config](const std::string& options) { config = CreateConfigFromString(options, config); })
            ->required();
    }

    { // GNM_DIRECTED
        auto* cmd =
            app.add_subcommand("gnm-directed", "Directed Erdos-Renyi Graph")->alias("gnm_directed")->callback([&] {
                config.generator = GeneratorType::GNM_DIRECTED;
            });
        add_option_n(cmd)->required();
        add_option_m(cmd)->required();
        add_option_self_loops(cmd);
    }

    { // GNM_UNDIRECTED
        auto* cmd = app.add_subcommand("gnm-undirected", "Undirected Erdos-Renyi Graph");
        cmd->alias("gnm_undirected");
        cmd->callback([&] { config.generator = GeneratorType::GNM_UNDIRECTED; });
        add_option_n(cmd)->required();
        add_option_m(cmd)->required();
        add_option_self_loops(cmd);
    }

    { // GNP_DIRECTED
        auto* cmd = app.add_subcommand("gnp-directed", "Directed Erdos-Renyi Graph");
        cmd->alias("gnp_directed");
        cmd->callback([&] { config.generator = GeneratorType::GNP_DIRECTED; });
        add_option_n(cmd)->required();
        add_option_p(cmd)->required();
        add_option_self_loops(cmd);
    }

    { // GNP_UNDIRECTED
        auto* cmd = app.add_subcommand("gnp-undirected", "Undirected Erdos-Renyi Graph");
        cmd->alias("gnp_undirected");
        cmd->callback([&] { config.generator = GeneratorType::GNP_UNDIRECTED; });
        add_option_n(cmd)->required();
        add_option_p(cmd)->required();
        add_option_self_loops(cmd);
    }

    { // RGG2D
        auto* cmd = app.add_subcommand("rgg2d", "2D Random Geometric Graph");
        cmd->alias("rgg_2d")->alias("rgg-2d");
        cmd->callback([&] { config.generator = GeneratorType::RGG_2D; });

        auto* params = cmd->add_option_group("Parameters");
        add_option_n(params);
        add_option_r(params);
        add_option_m(params);
        params->require_option(2);
        params->silent();
    }

    { // RGG3D
        auto* cmd = app.add_subcommand("rgg3d", "3D Random Geometric Graph");
        cmd->alias("rgg_3d")->alias("rgg-3d");
        cmd->callback([&] { config.generator = GeneratorType::RGG_3D; });

        auto* params = cmd->add_option_group("Parameters");
        add_option_n(params);
        add_option_r(params);
        add_option_m(params);
        params->require_option(2);
        params->silent();
    }

#ifdef KAGEN_CGAL_FOUND
    { // RDG2D
        auto* cmd = app.add_subcommand("rdg2d", "2D Random Delaunay Graph");
        cmd->alias("rdg_2d")->alias("rdg-2d");
        cmd->callback([&] { config.generator = GeneratorType::RDG_2D; });
        cmd->add_flag(
            "--periodic", config.periodic,
            "Enables the periodic boundary condition. Can yield unexpected results when using less than 9 PEs.");

        auto* params = cmd->add_option_group("Parameters");
        add_option_n(params);
        add_option_m(params);
        params->require_option(1);
        params->silent();
    }

    { // RDG3D
        auto* cmd = app.add_subcommand("rdg3d", "3D Random Delaunay Graph");
        cmd->alias("rdg_3d")->alias("rdg-3d");
        cmd->callback([&] { config.generator = GeneratorType::RDG_3D; });

        auto* params = cmd->add_option_group("Parameters");
        add_option_n(params);
        add_option_m(params);
        params->require_option(1);
        params->silent();
    }
#endif // KAGEN_CGAL_FOUND

    { // GRID 2D
        auto* cmd = app.add_subcommand("grid2d", "2D Grid Graph");
        cmd->alias("grid_2d")->alias("grid-2d");
        cmd->callback([&] { config.generator = GeneratorType::GRID_2D; });
        cmd->add_flag(
            "--periodic", config.periodic,
            "Enables the periodic boundary condition. Can yield unexpected results when using less than 9 PEs.");

        add_option_x(cmd)->required();
        add_option_y(cmd)->required();

        auto* params = cmd->add_option_group("Parameters");
        add_option_p(params);
        add_option_m(params);
        params->require_option(1);
        params->silent();
    }

    { // GRID 3D
        auto* cmd = app.add_subcommand("grid3d", "3D Grid Graph");
        cmd->alias("grid_3d")->alias("grid-3d");
        cmd->callback([&] { config.generator = GeneratorType::GRID_3D; });
        cmd->add_flag(
            "--periodic", config.periodic,
            "Enables the periodic boundary condition. Can yield unexpected results when using less than 9 PEs.");

        add_option_x(cmd)->required();
        add_option_y(cmd)->required();
        add_option_z(cmd)->required();

        auto* params = cmd->add_option_group("Parameters");
        add_option_p(params);
        add_option_m(params);
        params->require_option(1);
        params->silent();
    }

    { // PATH_DIRECTED
        auto* cmd = app.add_subcommand("path-directed", "Directed Path Graph")->alias("path_directed");
        cmd->callback([&] { config.generator = GeneratorType::PATH_DIRECTED; });
        cmd->add_flag(
            "--periodic", config.periodic,
            "Enables the periodic boundary condition. If enabled the generated path is a cycle.");
        cmd->add_flag(
            "--permute", config.permute,
            "Enables the permuation of vertices. If enabled the path is permuted following a pseudo-random "
            "permutation.");

        auto* params = cmd->add_option_group("Parameters");
        add_option_n(params);
        params->require_option(1);
        params->silent();
    }

    { // BA
        auto* cmd = app.add_subcommand("ba", "Barabassi Graph");
        cmd->callback([&] { config.generator = GeneratorType::BA; });
        add_option_directed(cmd);
        add_option_self_loops(cmd);

        auto* params = cmd->add_option_group("Parameters");
        add_option_n(params);
        add_option_m(params);
        add_option_min_deg(params);
        params->require_option(2);
        params->silent();
    }

    { // KRONECKER
        auto* cmd = app.add_subcommand("kronecker", "Kronecker Graph");
        cmd->callback([&] { config.generator = GeneratorType::KRONECKER; });
        add_option_self_loops(cmd);
        add_option_directed(cmd);
        add_option_n(cmd)->required();
        add_option_m(cmd)->required();
    }

    { // RHG
        auto* cmd = app.add_subcommand("rhg", "Random Hyperbolic Graph");
        cmd->callback([&] { config.generator = GeneratorType::RHG; });
        cmd->add_flag("--query-both", config.query_both, "Generate reverse cut edges communication-free (slow!)");
        cmd->add_flag("--hp-floats,!--no-hp-floats", config.hp_floats, "Use 80 bit floating point numbers");
        add_option_gamma(cmd)->required();

        auto* params = cmd->add_option_group("Parameters");
        add_option_n(params);
        add_option_avg_deg(params);
        add_option_m(params);
        params->require_option(2);
        params->silent();
    }

    { // RMAT
        auto* cmd = app.add_subcommand("rmat", "R-MAT Graph")->alias("r-mat");
        cmd->callback([&] { config.generator = GeneratorType::RMAT; });
        add_option_self_loops(cmd);
        add_option_directed(cmd);
        add_option_n(cmd);
        add_option_m(cmd);
        cmd->add_option("-a", config.rmat_a, "Probability for block a");
        cmd->add_option("-b", config.rmat_b, "Probability for block b");
        cmd->add_option("-c", config.rmat_c, "Probability for block c");
    }

    { // ImageMesh
        auto* cmd = app.add_subcommand("image", "Mesh graphs based on images")->alias("imagemesh")->alias("image-mesh");
        cmd->callback([&] { config.generator = GeneratorType::IMAGE_MESH; });
        cmd->add_option("--filename", config.image_mesh.filename, "Input image filename")
            ->required()
            ->check(CLI::ExistingFile);
        cmd->add_option("--weight-model", config.image_mesh.weight_model)
            ->transform(CLI::CheckedTransformer(GetImageMeshWeightModelMap()).description(""))
            ->description(R"(The following weight models are available:
  - l2:        sqrt(dR^2 + dG^2 + dB^2)
  - inv-l2:    \sqrt{3} * 255 + 1 - <l2>
  - ratio:     (Rmax / Rmin * Gmax / Gmin * Bmax / Bmin)
  - inv-ratio: 1 / <ratio>))");
        cmd->add_option("--weight-multiplier", config.image_mesh.weight_multiplier, "Multiplier for edge weights");
        cmd->add_option("--weight-offset", config.image_mesh.weight_offset, "Static offset for edge weights");
        cmd->add_option(
            "--min-weight-threshold", config.image_mesh.weight_min_threshold,
            "Only keep edges with weight more than this.");
        cmd->add_option(
            "--max-weight-threshold", config.image_mesh.weight_max_threshold,
            "Only keep edges with weight less than this.");
        cmd->add_option("--neighborhood", config.image_mesh.neighborhood, "Neighborhood size (4, 8 or 24)")
            ->check(CLI::IsMember({4, 8, 24}));
        cmd->add_option("--max-grid-x", config.image_mesh.max_grid_x, "Number of grid columns");
        cmd->add_option("--max-grid-y", config.image_mesh.max_grid_y, "Number of grid rows");
        cmd->add_option("--grid-x", config.image_mesh.grid_x, "Number of grid columns that are assigned to PEs");
        cmd->add_option("--grid-y", config.image_mesh.grid_y, "Number of grid rows that are assigned to PEs");
        cmd->add_option("--cols-per-pe", config.image_mesh.cols_per_pe, "Number of columns assigned to the same PE");
        cmd->add_option("--rows-per-pe", config.image_mesh.rows_per_pe, "Number of rows assigned to the same PE");
    }

    { // Graph from file
        auto* cmd =
            app.add_subcommand("file", "Loads a static graph from disk")->alias("staticgraph")->alias("static-graph");
        cmd->alias("file"); // @deprecated
        cmd->callback([&] { config.generator = GeneratorType::FILE; });
        cmd->add_option("--filename", config.input_graph.filename, "Input graph filename")
            ->required()
            ->check(CLI::ExistingFile);
        cmd->add_option("--distribution", config.input_graph.distribution)
            ->transform(CLI::CheckedTransformer(GetGraphDistributionMap()).description(""))
            ->description(R"(The following options for how to distribute the static graph across PEs are available:
  - balance-vertices: assign roughly the same number of nodes to each PE
  - balance-edges:    assign roughly the same number of edges to each PE by assigning consecutive vertices to a PE until the number of incident edges is >= m/<nproc>)");
        cmd->add_option("--input-format", config.input_graph.format)
            ->transform(CLI::CheckedTransformer(GetInputFormatMap()).description(""))
            ->description(R"(The following file formats are supported:
  - metis:          text format used by METIS
  - parhip:         binary format used by ParHIP
  - plain-edgelist: text file containing one edge per line, separated by spaces or tabs, starting at 0)");
    }

    // IO options
    app.add_option("-o,--output", config.output_graph.filename, "Output filename");
    app.add_option("-f,--output-format", config.output_graph.formats)
        ->transform(CLI::CheckedTransformer(GetOutputFormatMap()).description(""))
        ->description(R"(File formats for the generated graph, available formats are:
  - noop:            do not save the generated graph
  - edgelist:        text file containing the generated edges
  - binary-edgelist: binary file containing the generated edges
  - metis:           format used by METIS
  - hmetis:          format used by hMETIS
  - parhip:          binary format used by ParHIP
  - dot:             GraphViz format
  - xtrapulp:        format used by XtraPuLP
  - coordinates:     text file containing x y z coordinates)");
    app.add_option("--output-header", config.output_graph.header)
        ->transform(CLI::CheckedTransformer(GetOutputHeaderMap()).description(""))
        ->description(
            R"(When using distributed output: controls which PEs add a file header to their output file, possible values are:
  - never:  no PE outputs a file header
  - root:   only the root PE outputs a file header
  - always: every PE outputs a file header)");
    app.add_flag(
        "--distributed-output", [&config](auto) { config.output_graph.distributed = true; },
        "Output one file for each PE");
    app.add_flag(
        "--64",
        [&config](auto) {
            config.output_graph.width = 64;
            config.input_graph.width  = 64;
        },
        "Use 64 bit data types for the {binary-edge-list, xtrapulp} formats.");
    app.add_flag(
        "--32",
        [&config](auto) {
            config.output_graph.width = 32;
            config.input_graph.width  = 32;
        },
        "Use 32 bit data types for the {binary-edge-list, xtrapulp} formats.");
    app.add_flag(
        "--extension", config.output_graph.extension, "Always append a default extension to the output filename.");
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    // Parse parameters
    PGeneratorConfig config;
    CLI::App         app("KaGen: Karlsruhe Graph Generator");
    SetupCommandLineArguments(app, config);
    CLI11_PARSE(app, argc, argv);

    // Coordinates output format implies --coordinates
    if (std::find(config.output_graph.formats.begin(), config.output_graph.formats.end(), FileFormat::COORDINATES)
        != config.output_graph.formats.end()) {
        config.coordinates = true;
    }

    // If use more than one output format, always make output filenames distinct by appending the default extension
    if (config.output_graph.formats.size() > 1) {
        config.output_graph.extension = true;
    }

    // Run KaGen
    auto graph = Generate(config, GraphRepresentation::EDGE_LIST, MPI_COMM_WORLD);

    // Write resulting graph to disk
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const std::string base_filename = config.output_graph.filename;
    for (const FileFormat& format: config.output_graph.formats) {
        const auto& factory = GetGraphFormatFactory(format);

        // Append default filename
        const std::string filename   = (config.output_graph.extension && !factory->DefaultExtensions().empty())
                                           ? base_filename + "." + factory->DefaultExtensions().front()
                                           : base_filename;
        config.output_graph.filename = filename;

        GraphInfo info(graph, MPI_COMM_WORLD);
        auto      writer = factory->CreateWriter(config.output_graph, graph, info, rank, size);
        if (writer != nullptr) {
            WriteGraph(*writer.get(), config.output_graph, rank == ROOT && !config.quiet, MPI_COMM_WORLD);
        } else if (!config.quiet && rank == ROOT) {
            std::cout << "Warning: invalid file format " << format << " for writing; skipping\n";
        }
    }

    return MPI_Finalize();
}
