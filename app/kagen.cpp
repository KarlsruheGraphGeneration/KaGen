/*******************************************************************************
 * app/generate_kagen.cpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Daniel Funke <funke@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#include <iostream>

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/facade.h"
#include "kagen/io/io.h"
#include "kagen/tools/statistics.h"

#include "CLI11.h"

using namespace kagen;

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
    app.add_flag("-V,--validate-simple", config.validate_simple_graph, "Validate that the generated graph is simple");
    app.add_flag(
        "--skip-postprocessing", config.skip_postprocessing, "Skip postprocessing (repair graph, fix vertex ranges)");
    app.add_flag("-s,--seed", config.seed, "Seed for PRNG");
    auto* stats_group = app.add_option_group("Statistics output");
    stats_group
        ->add_option(
            "--stats", config.statistics_level, "Controls how much statistics on the generated graph gets calculated")
        ->transform(CLI::CheckedTransformer(GetStatisticsLevelMap()));
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
    app.add_flag("-C,--coordinates", config.coordinates, "Generate coordinates (geometric generators only)");

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
            "-p,--periodic", config.periodic,
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
        add_option_x(cmd)->required();
        add_option_y(cmd)->required();
    }

    { // GRID 3D
        auto* cmd = app.add_subcommand("grid3d", "3D Grid Graph");
        cmd->alias("grid_3d")->alias("grid-3d");
        cmd->callback([&] { config.generator = GeneratorType::GRID_3D; });
        add_option_x(cmd)->required();
        add_option_y(cmd)->required();
        add_option_z(cmd)->required();
    }

    { // BA
        auto* cmd = app.add_subcommand("ba", "Barabassi Graph");
        cmd->callback([&] { config.generator = GeneratorType::BA; });
        add_option_n(cmd)->required();
        add_option_min_deg(cmd)->required();
    }

    { // KRONECKER
        auto* cmd = app.add_subcommand("kronecker", "Kronecker Graph");
        cmd->callback([&] { config.generator = GeneratorType::KRONECKER; });
        add_option_n(cmd)->required();
        add_option_m(cmd)->required();
    }

    { // RHG
        auto* cmd = app.add_subcommand("rhg", "Random Hyperbolic Graph");
        cmd->callback([&] { config.generator = GeneratorType::RHG; });
        cmd->add_flag("--query-both", config.query_both, "(Experimental)");
        add_option_gamma(cmd)->required();

        auto* params = cmd->add_option_group("Parameters");
        add_option_n(params);
        add_option_avg_deg(params);
        add_option_m(params);
        params->require_option(2);
        params->silent();
    }

    // IO options
    app.add_option("-o,--output", config.output_file, "Output filename");
    app.add_option("-f,--output-format", config.output_format, "Output format")
        ->transform(CLI::CheckedTransformer(GetOutputFormatMap()));
    app.add_option("--output-header", config.output_header, "Output file header")
        ->transform(CLI::CheckedTransformer(GetOutputHeaderMap()));
    app.add_flag(
        "--distributed-output", [&config](auto) { config.output_single_file = false; }, "Output one file for each PE");

    // coordinates output format implies --coordinates
    if (config.output_format == OutputFormat::COORDINATES) {
        config.coordinates = true;
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    PEID rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Parse parameters
    PGeneratorConfig config;
    CLI::App         app("KaGen: Karlsruhe Graph Generator");
    SetupCommandLineArguments(app, config);
    CLI11_PARSE(app, argc, argv);

    // Only output on root if not in quiet mode
    const bool output = !config.quiet && rank == ROOT;
    if (output) {
        std::cout << "###############################################################################\n";
        std::cout << "#                         _  __      ____                                     #\n";
        std::cout << "#                        | |/ /__ _ / ___| ___ _ __                           #\n";
        std::cout << "#                        | ' // _` | |  _ / _ \\ '_ \\                          #\n";
        std::cout << "#                        | . \\ (_| | |_| |  __/ | | |                         #\n";
        std::cout << "#                        |_|\\_\\__,_|\\____|\\___|_| |_|                         #\n";
        std::cout << "#                         Karlsruhe Graph Generation                          #\n";
        std::cout << "#                                                                             #\n";
        std::cout << "###############################################################################\n";
        std::cout << config;
    }

    // Generate graph
    auto [edges, vertex_range, coordinates] = Generate(config, MPI_COMM_WORLD);

    // Output graph
    WriteGraph(config, edges, vertex_range, coordinates, MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}
