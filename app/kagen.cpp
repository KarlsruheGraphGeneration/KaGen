/*******************************************************************************
 * app/generate_kagen.cpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Daniel Funke <funke@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
//#define DEL_STATS 1

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/facade.h"
#include "kagen/io/io.h"

#include "CLI11.h"

using namespace kagen;

void SetupCommandLineArguments(CLI::App& app, PGeneratorConfig& config) {
    auto add_option_n = [&](CLI::App* cmd) {
        return cmd->add_option("-n,--num-nodes", config.n, "Number of nodes");
    };
    auto add_option_m = [&](CLI::App* cmd) {
        return cmd->add_option("-m,--num-edges", config.m, "Number of edges");
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
        return cmd->add_option("-x,--grid-x", config.grid_x, "Grid x dimension");
    };
    auto add_option_y = [&](CLI::App* cmd) {
        return cmd->add_option("-y,--grid-y", config.grid_y, "Grid y dimension");
    };
    auto add_option_z = [&](CLI::App* cmd) {
        return cmd->add_option("-z,--grid-z", config.grid_z, "Gird z dimension");
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
    app.add_flag("-S,--seed", config.seed, "Seed for PRNG");

    // Generator parameters
    app.add_flag("-k,--num-chunks", config.k, "Number of chunks used for graph generation");

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
        add_option_n(cmd)->required();
        add_option_r(cmd)->required();
    }
    { // RGG3D
        auto* cmd = app.add_subcommand("rgg3d", "3D Random Geometric Graph");
        cmd->alias("rgg_3d")->alias("rgg-3d");
        cmd->callback([&] { config.generator = GeneratorType::RGG_3D; });
        add_option_n(cmd)->required();
        add_option_r(cmd)->required();
    }
#ifdef KAGEN_CGAL_FOUND
    { // RDG2D
        auto* cmd = app.add_subcommand("rdg2d", "2D Random Delaunay Graph");
        cmd->alias("rdg_2d")->alias("rdg-2d");
        cmd->callback([&] { config.generator = GeneratorType::RDG_2D; });
        add_option_n(cmd)->required();
    }
    { // RDG3D
        auto* cmd = app.add_subcommand("rdg3d", "3D Random Delaunay Graph");
        cmd->alias("rdg_3d")->alias("rdg-3d");
        cmd->callback([&] { config.generator = GeneratorType::RDG_2D; });
        add_option_n(cmd)->required();
    }
#endif // KAGEN_CGAL_FOUND
    {  // GRID 2D
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
        // @todo
    }
    { // RHG
        auto* cmd = app.add_subcommand("rhg", "Random Hyperbolic Graph");
        cmd->callback([&] { config.generator = GeneratorType::RHG; });
        add_option_n(cmd)->required();
        add_option_gamma(cmd)->required();
        add_option_avg_deg(cmd)->required();
    }

    // IO options
    app.add_option("-o,--output", config.output_file, "Output filename");
    app.add_option("-f,--output-format", config.output_format, "Output format")
        ->check(CLI::IsMember(GetOutputFormatMap()))
        ->transform(CLI::CheckedTransformer(GetOutputFormatMap()).description(""));
    app.add_option("--output-header", config.output_header, "Output file header")
        ->check(CLI::IsMember(GetOutputHeaderMap()))
        ->transform(CLI::CheckedTransformer(GetOutputHeaderMap()).description(""));
    app.add_flag("--single-file", config.output_single_file, "Collect graph in a single file");
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
        std::cout << config;
    }

    // Generate graph
    if (output) {
        std::cout << "Generating graph ..." << std::endl;
    }
    auto [edges, vertex_range] = Generate(config);

    // Print statistics

    // Output graph
    if (output) {
        std::cout << "Writing graph ..." << std::endl;
    }
    WriteGraph(config, edges, vertex_range);

    MPI_Finalize();
    return 0;
}
