/*******************************************************************************
 * app/generate_kagen.cpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Daniel Funke <funke@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#include "kagen.h"

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/external_memory_facade.h"
#include "kagen/in_memory_facade.h"

#include <mpi.h>

#include "CLI11.h"
#include <cmath>
#include <iostream>
#include <limits>

using namespace kagen;

void PrintVersion() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == ROOT) {
        std::cout << BuildDescription() << '\n';
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
    auto set_hypergraph_generator = [&](const GeneratorType generator) {
        return [&config, generator] {
            config.generator     = generator;
            config.is_hypergraph = true;
        };
    };

    auto add_hypergraph_nm = [&](CLI::App* cmd) {
        add_option_n(cmd)->required();
        add_option_m(cmd)->required();
    };

    auto set_random_radius = [&config](const std::string&) {
        config.random_radius = true;
    };

    auto add_option_pin_budget = [&](CLI::App* cmd) {
        auto* exact_option =
            cmd->add_option("--pin-budget", config.size_dist_pin_budget, "Expected total number of pins");
        auto* log_option = cmd->add_option_function<SInt>(
            "--PB,--log-pin-budget",
            [&](const SInt exponent) {
                config.size_dist_pin_budget = static_cast<SInt>(1) << exponent;
                config.random_radius        = true;
            },
            "Logarithm of the expected total number of pins");

        exact_option->each(set_random_radius);
        exact_option->excludes(log_option);
        log_option->excludes(exact_option);

        return log_option;
    };

    auto add_hyperedge_radius_options = [&](CLI::App* cmd) {
        auto* radius_mode = cmd->add_option_group("Hyperedge radius mode");

        auto* fixed_radius =
            radius_mode->add_option("--fixed-radius,-r,--radius", config.r, "Use a constant hyperedge radius");

        auto* radius_exponent = radius_mode->add_option(
            "--radius-exponent", config.hyperedge_radius_exponent,
            "Use variable hyperedge radii with the given distribution exponent");

        radius_exponent->each(set_random_radius);

        auto* pin_budget = add_option_pin_budget(radius_mode);

        radius_mode->require_option(1);
        radius_mode->silent();

        auto* variable_radius_params = cmd->add_option_group("Variable hyperedge radius parameters");

        auto* min_radius = variable_radius_params->add_option(
            "--minr,--radius-lower-bound", config.min_hyperedge_radius,
            "Optional lower bound for variable hyperedge radius");

        auto* max_radius = variable_radius_params->add_option(
            "--maxr,--radius-upper-bound", config.max_hyperedge_radius,
            "Optional upper bound for variable hyperedge radius");

        auto* min_size = variable_radius_params->add_option(
            "--lower-size-bound", config.size_dist_lower_bound, "Optional lower bound on expected hyperedge size");

        auto* max_size = variable_radius_params->add_option(
            "--upper-size-bound", config.size_dist_upper_bound, "Optional upper bound on expected hyperedge size");

        auto* radius_quantile = variable_radius_params->add_option(
            "--rq,--radius-quantile", config.quantile,
            "Quantile used to tune the spatial grid for variable hyperedge radii");

        radius_quantile->check(CLI::Range(0.0, 1.0));

        min_radius->excludes(fixed_radius);
        max_radius->excludes(fixed_radius);
        min_size->excludes(fixed_radius);
        max_size->excludes(fixed_radius);
        radius_quantile->excludes(fixed_radius);

        min_radius->excludes(min_size);
        min_size->excludes(min_radius);

        max_radius->excludes(max_size);
        max_size->excludes(max_radius);

        return radius_mode;
    };

    auto add_partial_cell_mode = [&](CLI::App* cmd) {
        cmd->add_flag_callback(
            "--exact,-e", [&config]() { config.partial_cell_mode = PartialCellMode::GenerateAndCheck; },
            "Generate and check all vertices in boundary cells");
        cmd->add_flag_callback(
            "--approx-range", [&config]() { config.partial_cell_mode = PartialCellMode::EstimateByCoverageRange; },
            "Approximate boundary cells by sampling one random pin range per cell");
        cmd->add_flag_callback(
            "--approx-floyd,-a", [&config]() { config.partial_cell_mode = PartialCellMode::EstimateByCoverageFloyd; },
            "Approximate boundary cells by Floyd-sampling vertices uniformly in cell");
    };

    auto add_power_of_two_budget = [&](CLI::App* parent, const std::string& group_name, const std::string& exact_flags,
                                       const std::string& log_flags, double& destination,
                                       const std::string& quantity_name,
                                       const bool         required = false) -> std::pair<CLI::Option*, CLI::Option*> {
        auto* group        = parent->add_option_group(group_name);
        auto* exact_option = group->add_option(exact_flags, destination, "Expected total number of " + quantity_name)
                                 ->ignore_case(false);

        auto* log_option =
            group
                ->add_option_function<SInt>(
                    log_flags,
                    [&](const SInt exponent) {
                        if (exponent < 0
                            || exponent > static_cast<SInt>(std::numeric_limits<double>::max_exponent - 1)) {
                            throw CLI::ValidationError(log_flags, "exponent is outside the supported range");
                        }

                        destination = std::ldexp(1.0, static_cast<int>(exponent));
                    },
                    "Base-2 logarithm of the expected total number of " + quantity_name)
                ->ignore_case(false);

        exact_option->check(CLI::PositiveNumber);
        log_option->check(CLI::NonNegativeNumber);
        exact_option->excludes(log_option);
        log_option->excludes(exact_option);

        if (required) {
            group->require_option(1);
        }
        group->silent();

        return {exact_option, log_option};
    };

    auto add_edge_budget = [&](CLI::App* cmd, const bool required = false) {
        return add_power_of_two_budget(
            cmd, "Expected hyperedges", "--edge-budget,--eb", "--log-edge-budget,--EB", config.edge_budget,
            "hyperedges", required);
    };

    auto add_hyper_size_bounds = [&](CLI::App* cmd) {
        auto* group = cmd->add_option_group("Hyperedge size bounds");
        group->add_option("-l,--lower-size-bound", config.size_dist_lower_bound, "Lower bound for hyperedge sizes");
        group->add_option("-u,--upper-size-bound", config.size_dist_upper_bound, "Upper bound for hyperedge sizes");
        group->silent();
        return group;
    };

    auto add_explicit_hyperedge_sizes = [&](CLI::App* cmd, const std::string& explicit_description,
                                            const bool required) {
        const auto min_hyperedge_size = static_cast<SInt>(2);
        const auto max_hyperedge_size = std::numeric_limits<SInt>::max();

        auto* group = cmd->add_option_group("Hyperedge sizes");
        auto* lower =
            group->add_option("-l,--lower-size-bound", config.size_dist_lower_bound, "Lower bound for hyperedge sizes")
                ->check(CLI::Range(min_hyperedge_size, max_hyperedge_size));
        auto* upper =
            group->add_option("-u,--upper-size-bound", config.size_dist_upper_bound, "Upper bound for hyperedge sizes")
                ->check(CLI::Range(min_hyperedge_size, max_hyperedge_size));
        auto* sizes = group->add_option("--sizes", config.cigam_sizes, explicit_description)
                          ->expected(1, -1)
                          ->check(CLI::Range(min_hyperedge_size, max_hyperedge_size));

        sizes->excludes(lower);
        sizes->excludes(upper);

        if (required) {
            group->require_option(1, 2);
        }
        group->silent();
        return group;
    };

    auto add_approximation_options = [&](CLI::App* cmd) {
        auto* group = cmd->add_option_group("Approximation");
        group->add_flag("--approx", config.approx, "Use approximate generation");
        group->add_flag("--exact", [&config](auto) { config.approx = false; }, "Use exact generation");
        group->silent();
    };

    auto add_hgnm_size_distribution = [&](CLI::App* cmd) {
        auto* group = cmd->add_option_group("Hyperedge size distribution");
        auto* alpha = group->add_option("--alpha,--decay", config.size_dist_alpha, "Geometric distribution alpha")
                          ->check(CLI::Range(0.0, 1.0));
        auto* deterministic = group->add_flag(
            "--deterministic-size-distribution", config.size_dist_deterministic,
            "Use deterministic size counts instead of sampled counts");
        auto* size_decay =
            group
                ->add_option(
                    "--size-decay", config.size_decay, "Decay parameter for deterministic hyperedge size counts")
                ->check(CLI::Range(0.0, 1.0));
        auto* pin_budget = add_option_pin_budget(group);

        pin_budget->excludes(alpha);
        pin_budget->excludes(deterministic);
        pin_budget->excludes(size_decay);
        size_decay->needs(deterministic);

        group->silent();
        return group;
    };

    auto add_hgnp_probabilities = [&](CLI::App* cmd) {
        auto* probabilities = cmd->add_option_group("Probabilities");
        add_option_p(probabilities);
        probabilities
            ->add_option(
                "--size-probs", config.size_probabilities,
                "Comma-separated probabilities by hyperedge size, starting at lower-size-bound")
            ->delimiter(',');
        probabilities
            ->add_option(
                "--size-probs-file", config.size_probabilities_file,
                "File with one probability per line, starting at lower-size-bound")
            ->check(CLI::ExistingFile);
        probabilities
            ->add_option(
                "--size-expected", config.size_expected_counts,
                "Expected hyperedge counts by size, starting at lower-size-bound")
            ->delimiter(',');
        probabilities
            ->add_option(
                "--size-expected-file", config.size_expected_counts_file,
                "File with one expected hyperedge count per line, starting at lower-size-bound")
            ->check(CLI::ExistingFile);

        auto* budget_params = probabilities->add_option_group("Budget Parameters");
        add_edge_budget(budget_params);
        add_option_pin_budget(budget_params);
        budget_params
            ->add_option("--sd,--size-decay", config.size_decay, "Geometric decay parameter for expected counts")
            ->check(CLI::Range(0.0, 1.0));

        probabilities->require_option(1);
        probabilities->silent();
        return probabilities;
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
    app.add_option(
           "--external-k,--external-num-chunks", config.external.num_chunks,
           "Number of chunks for generating the graph in external memory mode. Set to '1' to disable external "
           "memory "
           "mode.")
        ->capture_default_str();
    app.add_option(
           "--external-T,--external-tmp-directory", config.external.tmp_directory,
           "Directory for temporary buffer files (requires free space between ~1x to ~2x the final graph "
           "size).")
        ->capture_default_str();
    app.add_flag(
        "!--external-avoid-extra-writes", config.external.cache_aggregated_chunks,
        "Do not cache aggregated chunks in extra buffer files (less IO writes, more extra work).");

    app.add_flag("-q,--quiet", config.quiet, "Quiet mode");
    app.add_flag("-v,--version", [&](auto) { PrintVersion(); }, "Print KaGen version")->trigger_on_parse();
    app.add_flag("-V,--validate", config.validate_simple_graph)
        ->description(
            R"(Validate that the generated graph is undirected and does not have duplicated edges or inconsistent edge weights.
This is mostly useful for experimental graph generators or when using KaGen to load graphs from disk.)");
    app.add_flag(
        "--skip-postprocessing", config.skip_postprocessing,
        "Skip postprocessing (repair inconsistency due to floating point inaccuracies etc.)");
    app.add_option("-s,--seed", config.seed, "Seed for PRNG (must be the same on all PEs)");
    app.add_flag("--debug", config.debug, "Enable debug output");
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
        "Controls how much statistics on the generated graph gets calculated; pass flag multiple times to "
        "increase "
        "extend");
    stats_group->silent();

    // Generator parameters
    app.add_option("-k,--num-chunks", config.k, "Number of chunks used for graph generation");
    app.add_option(
        "--automatic-num-chunks-imbalance-threshold", config.max_vertex_imbalance,
        "Controls the trade-off between vertex imbalance and number of chunks when deducing the number of "
        "chunks "
        "automatically");
    app.add_flag("-C,--coordinates", config.coordinates, "Generate coordinates (geometric generators only)");

    // Edge weight generator parameters
    app.add_option(
           "--edgeweights-generator", config.edge_weights.generator_type,
           "Generator type to use for generating edge weights.")
        ->transform(CLI::CheckedTransformer(GetEdgeWeightGeneratorTypeMap(), CLI::ignore_case));
    app.add_option(
        "--edgeweights-range-begin", config.edge_weights.weight_range_begin,
        "(Included) begin of weight range used for edge weights, i.e., minimum edge weight.");
    app.add_option(
        "--edgeweights-range-end", config.edge_weights.weight_range_end,
        "(Excluded) end of weight range to be used for edge weights.");

    // Vertex weight generator parameters
    app.add_option(
           "--vertexweights-generator", config.vertex_weights.generator_type,
           "Generator type to use for generating vertex weights.")
        ->transform(CLI::CheckedTransformer(GetVertexWeightGeneratorTypeMap(), CLI::ignore_case));
    app.add_option(
        "--vertexweights-range-begin", config.vertex_weights.weight_range_begin,
        "(Included) begin of weight range used for vertex weights, i.e., minimum vertex weight.");
    app.add_option(
        "--vertexweights-range-end", config.vertex_weights.weight_range_end,
        "(Excluded) end of weight range to be used for vertex weights.");

    // Permute Options
    app.add_flag(
        "--permute", config.permute,
        "Enables the permuation of vertices. If enabled the graph is permuted following a pseudo-random "
        "permutation.");

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

    { // Hyper RGG2D
        auto* cmd = app.add_subcommand("hrgg2d", "2D Random Geometric Graph adapted for Hypergraphs");
        cmd->alias("hrgg_2d")->alias("hrgg-2d");
        cmd->callback(set_hypergraph_generator(GeneratorType::HRGG_2D));

        auto* params = cmd->add_option_group("Parameters");
        add_hypergraph_nm(params);
        params->silent();

        auto* hyperedge_dist_options = cmd->add_option_group("Hyperedge Distribution options");

        add_hyperedge_radius_options(hyperedge_dist_options);

        add_partial_cell_mode(cmd);
    }

    { // Hyper RHG
        auto* cmd = app.add_subcommand("hrhg", "Random Hyperbolic Hypergraph");
        cmd->callback(set_hypergraph_generator(GeneratorType::H_RHG));

        cmd->add_flag("--hp-floats,!--no-hp-floats", config.hp_floats, "Use 80 bit floating point numbers");
        add_option_gamma(cmd)->required();

        auto* params = cmd->add_option_group("Parameters");
        add_hypergraph_nm(params);
        params->silent();

        auto* hyperedge_dist_options = cmd->add_option_group("Hyperedge Distribution options");

        add_hyperedge_radius_options(hyperedge_dist_options);

        add_partial_cell_mode(cmd);
    }

    { // Hyper GNM
        auto* cmd = app.add_subcommand("hgnm", "Erdos-Renyi Hypergraph G(n,m)");
        cmd->alias("h-gnm")->alias("h_gnm");
        cmd->callback(set_hypergraph_generator(GeneratorType::H_GNM));

        auto* params = cmd->add_option_group("Parameters");
        add_hypergraph_nm(params);
        add_hyper_size_bounds(params);
        add_hgnm_size_distribution(params);
        add_approximation_options(cmd);

        cmd->add_flag("--fast", config.allow_duplicates, "Do not check for duplicate edges in order to speed up");
        params->silent();
    }

    { // Hyper GNP
        auto* cmd = app.add_subcommand("hgnp", "Erdos-Renyi Hypergraph G(n,p)");
        cmd->alias("h-gnp")->alias("h_gnp")->alias("hyper_gnp")->alias("hyper-gnp");
        cmd->callback(set_hypergraph_generator(GeneratorType::H_GNP));

        auto* params = cmd->add_option_group("Parameters");
        add_option_n(params)->required();
        add_explicit_hyperedge_sizes(params, "Explicit HGNP hyperedge sizes", false);
        add_approximation_options(cmd);

        cmd->add_flag("--fast", config.allow_duplicates, "Do not check for duplicate edges in order to speed up");
        params->silent();

        add_hgnp_probabilities(cmd);
    }
    { // CIGAM
        auto* cmd = app.add_subcommand("cigam", "CIGAM Hypergraph Generator");

        cmd->callback(set_hypergraph_generator(GeneratorType::CIGAM));

        add_option_n(cmd)->required();

        auto* size_group = add_explicit_hyperedge_sizes(cmd, "Explicit CIGAM hyperedge sizes", true);
        size_group->get_option("-l")->description("Minimum hyperedge size");
        size_group->get_option("-u")->description("Maximum hyperedge size");

        auto* model_group = cmd->add_option_group("CIGAM model");

        model_group->add_option("--lambda", config.cigam_lambda, "Truncated-exponential prestige parameter")
            ->required()
            ->check(CLI::PositiveNumber);

        model_group->add_option("--c", config.cigam_c, "Density parameters, one per layer")
            ->required()
            ->expected(1, -1)
            ->check(
                CLI::Validator(
                    [](const std::string& value) {
                        const long double c = std::stold(value);

                        if (!std::isfinite(c) || c <= 1.0L) {
                            return std::string(
                                "CIGAM density parameters "
                                "must be finite and greater than 1");
                        }

                        return std::string{};
                    },
                    "finite c > 1"));

        model_group
            ->add_option(
                "--breakpoints", config.cigam_breakpoints,
                "Layer breakpoints in (0,1], "
                "strictly increasing, last equal to 1")
            ->required()
            ->expected(1, -1)
            ->check(
                CLI::Validator(
                    [](const std::string& value) {
                        const long double breakpoint = std::stold(value);

                        if (!std::isfinite(breakpoint) || breakpoint <= 0.0L || breakpoint > 1.0L) {
                            return std::string(
                                "CIGAM breakpoints must be "
                                "finite and lie in (0,1]");
                        }

                        return std::string{};
                    },
                    "finite 0 < H <= 1"));

        model_group->silent();

        auto* mode_group = cmd->add_option_group("Generation mode");

        mode_group->add_option("--mode", config.cigam_mode, "CIGAM generation mode")
            ->transform(CLI::CheckedTransformer(GetCIGAMModeMap(), CLI::ignore_case))
            ->description(
                R"(Selects the CIGAM generation strategy:
  - paper:  sampled IID prestige ranks and paper-reference generation
  - exact:  deterministic quantile ranks with binomial block counts
  - approx: deterministic quantile ranks with Poissonized range counts)")
            ->capture_default_str();

        mode_group->silent();

        /*
         * Target output volume.
         *
         * These two groups are individually optional. Validation after parsing
         * determines whether the chosen combination is supported.
         */
        const auto [edge_budget, log_edge_budget] = add_edge_budget(cmd);

        /*
         * Neither target is mandatory at CLI parsing time because paper mode
         * may be controlled directly through c.
         */
        edge_budget->group("");
        log_edge_budget->group("");

        /*
         * Optional size weighting when multiple hyperedge sizes are used.
         */
        cmd->add_option("--size-decay", config.size_decay, "Geometric weighting across configured hyperedge sizes")
            ->check(CLI::Range(0.0, 1.0));

        cmd->add_flag("--fast", config.allow_duplicates, "Skip duplicate checking where supported");
    }
#ifdef KAGEN_CGAL_FOUND
    { // RDG2D
        auto* cmd = app.add_subcommand("rdg2d", "2D Random Delaunay Graph");
        cmd->alias("rdg_2d")->alias("rdg-2d");
        cmd->callback([&] { config.generator = GeneratorType::RDG_2D; });
        cmd->add_flag(
            "--periodic", config.periodic,
            "Enables the periodic boundary condition. Can yield unexpected results when using less than 9 "
            "PEs.");

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
            "Enables the periodic boundary condition. Can yield unexpected results when using less than 9 "
            "PEs.");

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
            "Enables the periodic boundary condition. Can yield unexpected results when using less than 9 "
            "PEs.");

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
        cmd->add_option("--redistribution", config.redistribution)
            ->transform(CLI::CheckedTransformer(GetGraphRedistributionMap()).description(""))
            ->description(R"(How to distribute the generated graph across PEs:
  - balance-vertices: assign roughly the same number of vertices to each PE
  - balance-edges:    assign roughly the same number of edges to each PE)");
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
            ->description(
                R"(The following options for how to distribute the static graph across PEs are available:
  - balance-vertices: assign roughly the same number of nodes to each PE
  - balance-edges:    assign roughly the same number of edges to each PE by assigning consecutive vertices to a PE until the number of incident edges is >= m/<nproc>
  - explicit:         explicitly specify the number of vertices on each PE through a text file specified via the --explicit-distribution=<filename> option)");
        cmd->add_option(
            "--explicit-distribution", config.input_graph.explicit_distribution_filename,
            "A text file containing the number of vertices on each PE, one line per PE. Only used when "
            "--distribution=explicit.");
        cmd->add_flag(
            "--explicit-distribution-is-prefix-sum", config.input_graph.explicit_distribution_is_prefix_sum,
            "Interpret the explicit distribution as prefix sum instead of vertices on each PE.");
        cmd->add_option("--input-format", config.input_graph.format)
            ->transform(CLI::CheckedTransformer(GetInputFormatMap()).description(""))
            ->description(R"(The following file formats are supported:
  - metis:          text format used by METIS
  - parhip:         binary format used by ParHIP
  - plain-edgelist: text file containing one edge per line, separated by spaces or tabs, starting at 0)");
        cmd->add_flag("--drop-vertex-weights", config.input_graph.drop_vertex_weights, "Drop vertex weights");
        cmd->add_flag("--drop-edge-weights", config.input_graph.drop_edge_weights, "Drop edge weights");

        cmd->add_option("--input-vtx-width", config.input_graph.vtx_width, "")->capture_default_str();
        cmd->add_option("--input-adjncy-width", config.input_graph.adjncy_width, "")->capture_default_str();
        cmd->add_option("--input-vwgt-width", config.input_graph.vwgt_width, "")->capture_default_str();
        cmd->add_option("--input-adjwgt-width", config.input_graph.adjwgt_width, "")->capture_default_str();
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

    app.add_option("--vtx-width", config.output_graph.vtx_width, "")->capture_default_str();
    app.add_option("--adjncy-width", config.output_graph.adjncy_width, "")->capture_default_str();
    app.add_option("--vwgt-width", config.output_graph.vwgt_width, "")->capture_default_str();
    app.add_option("--adjwgt-width", config.output_graph.adjwgt_width, "")->capture_default_str();

    auto set_all_widths = [&](const int width) {
        return [&config, width](auto) {
            config.output_graph.width = width;
            config.input_graph.width  = width;

            config.input_graph.vtx_width    = width;
            config.input_graph.adjncy_width = width;
            config.input_graph.vwgt_width   = width;
            config.input_graph.adjwgt_width = width;

            config.output_graph.vtx_width    = width;
            config.output_graph.adjncy_width = width;
            config.output_graph.vwgt_width   = width;
            config.output_graph.adjwgt_width = width;
        };
    };

    app.add_flag("--64", set_all_widths(64), "Use 64 bit data types for binary formats (input and output).");
    app.add_flag("--32", set_all_widths(32), "Use 32 bit data types for binary formats (input and output).");

    app.add_flag(
        "--extension", config.output_graph.extension, "Always append a default extension to the output filename.");
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    try {
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

        // If use more than one output format, always make output filenames distinct by appending the default
        // extension
        if (config.output_graph.formats.size() > 1) {
            config.output_graph.extension = true;
        }

        if (config.external.num_chunks > 1) {
            GenerateExternalMemoryToDisk(config, MPI_COMM_WORLD);
        } else {
            GenerateInMemoryToDisk(config, MPI_COMM_WORLD);
        }

        MPI_Finalize();
        return EXIT_SUCCESS;
    } catch (const CLI::ParseError& e) {
        std::cerr << "Fatal error: " << e.what() << '\n';
    } catch (const std::exception& e) {
        std::cerr << "Fatal std::exception: " << e.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    } catch (...) {
        std::cerr << "Fatal non-std exception\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
    return EXIT_FAILURE;
}