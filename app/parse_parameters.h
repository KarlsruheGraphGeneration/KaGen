/*******************************************************************************
 * app/parse_parameters.cpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _PARSE_PARAMETERS_H_
#define _PARSE_PARAMETERS_H_

#include "definitions.h"
#include "generator_config.h"
#include "tools/arg_parser.h"
#include <iostream>

#include <mpi.h>
#include <string.h>

namespace kagen {

void ParseParameters(int argn, char** argv, PEID, PEID size, PGeneratorConfig& generator_config) {
    ArgParser args(argn, argv);

    // Generator
    generator_config.generator = args.Get<std::string>("gen", "");

    if (args.IsSet("help") || argn < 2) {
        if (generator_config.generator == "gnm_undirected" || generator_config.generator == "gnm_directed") {
            std::cout << "================================================" << std::endl;
            std::cout << "========== Erdos-Renyi Graphs G(n,m) ===========" << std::endl;
            std::cout << "================================================" << std::endl;
            std::cout << "Parameters:" << std::endl;
            std::cout << "-n\t\t<number of vertices as a power of two>" << std::endl;
            std::cout << "-m\t\t<number of edges as a power of two>" << std::endl;
            std::cout << "-k\t\t<number of chunks>" << std::endl;
            std::cout << "-seed\t\t<seed for PRNGs>" << std::endl;
            std::cout << "-self_loops" << std::endl;
            std::cout << "\nExample:" << std::endl;
            std::cout << "mpirun -n 16 ./build/app/kagen -gen gnm_directed -n 20 -m 22 -self_loops -output tmp"
                      << std::endl;
        } else if (generator_config.generator == "gnp_undirected" || generator_config.generator == "gnp_directed") {
            std::cout << "================================================" << std::endl;
            std::cout << "========== Erdos-Renyi Graphs G(n,p) ===========" << std::endl;
            std::cout << "================================================" << std::endl;
            std::cout << "Parameters:" << std::endl;
            std::cout << "-n\t\t<number of vertices as a power of two>" << std::endl;
            std::cout << "-p\t\t<edge probability>" << std::endl;
            std::cout << "-k\t\t<number of chunks>" << std::endl;
            std::cout << "-seed\t\t<seed for PRNGs>" << std::endl;
            std::cout << "-self_loops" << std::endl;
            std::cout << "\nExample:" << std::endl;
            std::cout << "mpirun -n 16 ./build/app/kagen -gen gnp_directed -n 20 -p 0.001 -self_loops -output tmp"
                      << std::endl;
        } else if (generator_config.generator == "rgg_2d" || generator_config.generator == "rgg_3d") {
            std::cout << "================================================" << std::endl;
            std::cout << "======= Random Geometric Graphs RGG(n,d) ========" << std::endl;
            std::cout << "================================================" << std::endl;
            std::cout << "Parameters:" << std::endl;
            std::cout << "-n\t\t<number of vertices as a power of two>" << std::endl;
            std::cout << "-r\t\t<radius for vertices to be connected> (r <= 1.0)>" << std::endl;
            std::cout << "-k\t\t<number of chunks>" << std::endl;
            std::cout << "-seed\t\t<seed for PRNGs>" << std::endl;
            std::cout << "\nExample:" << std::endl;
            std::cout << "mpirun -n 16 ./build/app/kagen -gen rgg_3d -n 20 -r 0.00275 -output tmp" << std::endl;
        } else if (generator_config.generator == "rdg_2d" || generator_config.generator == "rdg_3d") {
            std::cout << "================================================" << std::endl;
            std::cout << "======== Random Delaunay Graphs RDG(n) =========" << std::endl;
            std::cout << "================================================" << std::endl;
            std::cout << "Parameters:" << std::endl;
            std::cout << "-n\t\t<number of vertices as a power of two>" << std::endl;
            std::cout << "-k\t\t<number of chunks>" << std::endl;
            std::cout << "-seed\t\t<seed for PRNGs>" << std::endl;
            std::cout << "\nExample:" << std::endl;
            std::cout << "mpirun -n 16 ./build/app/kagen -gen rdg_3d -n 20 -output tmp" << std::endl;
        } else if (generator_config.generator == "ba") {
            std::cout << "================================================" << std::endl;
            std::cout << "======= Barabassi-Albert Graphs BA(n,d) ========" << std::endl;
            std::cout << "================================================" << std::endl;
            std::cout << "Parameters:" << std::endl;
            std::cout << "-n\t\t<number of vertices as a power of two>" << std::endl;
            std::cout << "-md\t\t<min degree for each vertex>" << std::endl;
            std::cout << "-k\t\t<number of chunks>" << std::endl;
            std::cout << "-seed\t\t<seed for PRNGs>" << std::endl;
            std::cout << "\nExample:" << std::endl;
            std::cout << "mpirun -n 16 ./build/app/kagen -gen ba -n 20 -md 4 -output tmp" << std::endl;
        } else if (generator_config.generator == "rhg") {
            std::cout << "Parameters for Random Hyperbolic Graphs RHG(n,gamma,d)" << std::endl;
            std::cout << "================================================" << std::endl;
            std::cout << "=== Random Hyperbolic Graphs RHG(n,gamma,d) ====" << std::endl;
            std::cout << "================================================" << std::endl;
            std::cout << "Parameters:" << std::endl;
            std::cout << "-n\t\t<number of vertices as a power of two>" << std::endl;
            std::cout << "-gamma\t\t<power-law exponent>" << std::endl;
            std::cout << "-d\t\t<average degree>" << std::endl;
            std::cout << "-k\t\t<number of chunks>" << std::endl;
            std::cout << "-seed\t\t<seed for PRNGs>" << std::endl;
            std::cout << "\nExample:" << std::endl;
            std::cout << "mpirun -n 16 ./build/app/kagen -gen rhg -n 20 -d 8 -gamma 2.2 -output tmp" << std::endl;
        } else if (generator_config.generator == "rmat") {
            std::cout << "Parameters for Kronecker Graphs RMAT(n,m)" << std::endl;
            std::cout << "================================================" << std::endl;
            std::cout << "=========== Kronecker Graphs RMAT(n,m) =========" << std::endl;
            std::cout << "================================================" << std::endl;
            std::cout << "Parameters:" << std::endl;
            std::cout << "-n\t\t<number of vertices as a power of two>" << std::endl;
            std::cout << "-m\t\t<number of edges as a power of two>" << std::endl;
            std::cout << "-k\t\t<number of chunks>" << std::endl;
            std::cout << "-seed\t\t<seed for PRNGs>" << std::endl;
            std::cout << "\nExample:" << std::endl;
            std::cout << "mpirun -n 16 ./build/app/kagen -gen rmat -n 20 -m 22 -output tmp" << std::endl;
        } else if (generator_config.generator == "grid") {
            std::cout << "Parameters for 2D/3D Grid Graphs G(x,y(,z),periodic)" << std::endl;
            std::cout << "================================================" << std::endl;
            std::cout << "=========== Grid Graphs G(x,y(,z)) ================" << std::endl;
            std::cout << "================================================" << std::endl;
            std::cout << "Parameters:" << std::endl;
            std::cout << "-x\t\t<size of first dimension>" << std::endl;
            std::cout << "-y\t\t<size of second dimension>" << std::endl;
            std::cout << "-z\t\t<size of third dimension>" << std::endl;
            std::cout << "-p\t\t<probability of edge insertion>" << std::endl;
            std::cout << "-periodic\t\t<use periodic boundary condition>" << std::endl;
            std::cout << "-k\t\t<number of chunks>" << std::endl;
            std::cout << "-seed\t\t<seed for PRNGs>" << std::endl;
            std::cout << "\nExample:" << std::endl;
            std::cout << "mpirun -n 16 ./build/app/kagen -gen grid -x 16 -y 16 -output tmp" << std::endl;
        } else {
            std::cout << "================================================" << std::endl;
            std::cout << "==================== KaGen =====================" << std::endl;
            std::cout << "================================================" << std::endl;
            std::cout << "Usage:\t\t\tmpirun -n <num_proc> ./kagen -gen <generator> [additional parameters]"
                      << std::endl;
            std::cout << "Generators:\t\tgnm_directed|gnm_undirected|gnp_directed|gnp_undirected|rgg_2d|rgg_3d|rdg_2d|"
                         "rdg_3d|ba|rhg"
                      << std::endl;
            std::cout << "Additional help:\t./kagen -gen <generator> -help" << std::endl;
        }

        std::cout << std::endl;
        std::cout << "Output options:" << std::endl;
        std::cout << "-output\t\t<output file>" << std::endl;
        std::cout << "-omit_header\t\tIf set, output the graph without header line" << std::endl;
        std::cout << "-single_file\t\tIf set, if set, write graph to a single file" << std::endl;
        exit(0);
    }

    // Nodes
    bool exact_n = args.IsSet("exact_n");
    if (exact_n)
        generator_config.n = args.Get<ULONG>("n", 100);
    else
        generator_config.n = (ULONG)1 << args.Get<ULONG>("n", 3);

    // Blocks
    generator_config.k = args.Get<ULONG>("k", size);

    // RNG
    generator_config.seed        = args.Get<ULONG>("seed", 1);
    generator_config.hash_sample = args.Get<bool>("hash_sample", false);
    generator_config.use_binom   = args.IsSet("binom");

    // I/O
    generator_config.output_file  = args.Get<std::string>("output", "out");
    generator_config.debug_output = args.Get<std::string>("debug", "dbg");
    generator_config.dist_size    = args.Get<ULONG>("dist", 10);

    auto output_format = args.Get<std::string>("output_format", "edge_list");
    if (output_format == "edge_list") {
        generator_config.output_format = OutputFormat::EDGE_LIST;
    } else if (output_format == "binary_edge_list") {
        generator_config.output_format = OutputFormat::BINARY_EDGE_LIST;
    } else {
        std::cout << "Invalid output format: " << output_format << std::endl;
        std::cout << "Supported output formats: edge_list, binary_edge_list" << std::endl;
        std::exit(0);
    }
    generator_config.output_header      = !args.IsSet("omit_header");
    generator_config.output_single_file = args.IsSet("single_file");

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout << "IO options:" << std::endl;
        std::cout << "- Output filename: " << generator_config.output_file << std::endl;
        std::cout << "- Output format: ";
        switch (generator_config.output_format) {
            case OutputFormat::EDGE_LIST:
                std::cout << "edge list";
                break;

            case OutputFormat::BINARY_EDGE_LIST:
                std::cout << "binary edge list";
                break;

            default:
                std::cout << "undefined";
        }
        std::cout << std::endl;
        std::cout << "- Output header: " << generator_config.output_header << std::endl;
        std::cout << "- Output to a single file: " << generator_config.output_single_file << std::endl;
    }

    // Edges
    bool exact_m = args.IsSet("exact_m");
    if (exact_m)
        generator_config.m = args.Get<ULONG>("m", 0);
    else
        generator_config.m = (ULONG)1 << args.Get<ULONG>("m", 0);
    generator_config.p          = args.Get<double>("p", 0.0);
    generator_config.self_loops = args.IsSet("self_loops");

    // Radius/Edges
    generator_config.r = args.Get<double>("r", 0.125);

    // Average degree
    generator_config.avg_degree = args.Get<double>("d", 5.0);
    generator_config.plexp      = args.Get<double>("gamma", 2.6);

    // RHG
    generator_config.thres      = args.Get<ULONG>("t", 0);
    generator_config.query_both = args.Get<bool>("qb", false);

    // BA
    generator_config.min_degree = args.Get<ULONG>("md", 4);

    // GRID
    generator_config.grid_x   = args.Get<ULONG>("x", 1);
    generator_config.grid_y   = args.Get<ULONG>("y", 1);
    generator_config.grid_z   = args.Get<ULONG>("z", 1);
    generator_config.periodic = args.IsSet("periodic");

    // Floating-point precision
    generator_config.precision = args.Get<ULONG>("prec", 32);

    // Sampling algorithm
    generator_config.base_size = (ULONG)1 << args.Get<ULONG>("sk", 8);
    generator_config.hyp_base  = (ULONG)1 << args.Get<ULONG>("hk", 8);

    // Benchmarks
    generator_config.iterations = args.Get<ULONG>("i", 1);
}

} // namespace kagen
#endif
