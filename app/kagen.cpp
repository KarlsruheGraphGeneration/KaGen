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

#include "kagen/generator_config.h"
#include "kagen/io/generator_io.h"
#include "kagen/postprocessing.h"
#include "kagen/tools/benchmark.h"
#include "kagen/tools/timer.h"
#include "parse_parameters.h"

#if KAGEN_CGAL_FOUND
    #include "kagen/generators/geometric/delaunay/delaunay_2d.h"
    #include "kagen/generators/geometric/delaunay/delaunay_3d.h"
#endif // KAGEN_CGAL_FOUND

#include "kagen/generators/barabassi/barabassi.h"
#include "kagen/generators/geometric/rgg/rgg_2d.h"
#include "kagen/generators/geometric/rgg/rgg_3d.h"
#include "kagen/generators/gnm/gnm_directed.h"
#include "kagen/generators/gnm/gnm_undirected.h"
#include "kagen/generators/gnp/gnp_directed.h"
#include "kagen/generators/gnp/gnp_undirected.h"
#include "kagen/generators/grid/grid_2d.h"
#include "kagen/generators/grid/grid_3d.h"
#include "kagen/generators/hyperbolic/hyperbolic.h"
#include "kagen/generators/kronecker/kronecker.h"

using namespace kagen;

void OutputParameters(PGeneratorConfig& config, const PEID /* rank */, const PEID size) {
    std::cout << "________________________________________________________________________________\n";
    std::cout << "                       _    _             __                                    \n";
    std::cout << "                       /  ,'            /    )                                  \n";
    std::cout << "----------------------/_.'-------__----/----------__-----__---------------------\n";
    std::cout << "                     /  \\      /   )  /  --,    /___)  /   )                    \n";
    std::cout << "____________________/____\\____(___(__(____/____(___ __/___/_____________________\n";
    std::cout << "\n";
    std::cout << "Input Parameters:\n";
    std::cout << "+-- Generator:\n";
    std::cout << "|   +-- Type: ......................... " << GeneratorToString(config.generator) << "\n";
    std::cout << "|   +-- n: ............................ " << config.n << "\n";

    switch (config.generator) {
        case Generator::GNM_DIRECTED:
        case Generator::GNM_UNDIRECTED:
        case Generator::GNP_DIRECTED:
        case Generator::GNP_UNDIRECTED:
            std::cout << "|   +-- m: ............................ " << config.m << "\n";
            std::cout << "|   +-- p: ............................ " << config.p << "\n";
            std::cout << "|   +-- k: ............................ " << config.k << "\n";
            break;

        case Generator::RGG_2D:
        case Generator::RGG_3D:
            std::cout << "|   +-- r: ............................ " << config.r << "\n";
            std::cout << "|   +-- k: ............................ " << config.k << "\n";
            break;

        case Generator::RDG_2D:
        case Generator::RDG_3D:
            std::cout << "|   +-- k: ............................ " << config.k << "\n";
            break;

        case Generator::RHG:
            std::cout << "|   +-- d: ............................ " << config.avg_degree << "\n";
            std::cout << "|   +-- gamma: ........................ " << config.plexp << "\n";
            std::cout << "|   +-- k: ............................ " << config.k << "\n";
            break;

        case Generator::BA:
            std::cout << "|   +-- d: ............................ " << config.avg_degree << "\n";
            std::cout << "|   +-- k: ............................ " << config.k << "\n";
            break;

        case Generator::GRID_2D:
            std::cout << "|   +-- row: .......................... " << config.grid_x << "\n";
            std::cout << "|   +-- col: .......................... " << config.grid_y << "\n";
            std::cout << "|   +-- p: ............................ " << config.p << "\n";
            std::cout << "|   +-- k: ............................ " << config.k << "\n";
            break;

        case Generator::GRID_3D:
            std::cout << "|   +-- x: ............................ " << config.grid_x << "\n";
            std::cout << "|   +-- y: ............................ " << config.grid_y << "\n";
            std::cout << "|   +-- z: ............................ " << config.grid_z << "\n";
            std::cout << "|   +-- p: ............................ " << config.p << "\n";
            std::cout << "|   +-- k: ............................ " << config.k << "\n";
            break;

        case Generator::KRONECKER:
            // @todo
            break;

        default:
            std::cerr << "Error: unknown generator type\n";
            std::exit(1);
    }

    std::cout << "|   +-- s: ............................ " << config.seed << "\n";
    std::cout << "|   `-- P: ............................ " << size << "\n";
    std::cout << "+-- IO\n";
    std::cout << "|   +-- Filename: ..................... " << config.output_file << "\n";
    std::cout << "|   +-- Format: ....................... " << OutputFormatToString(config.output_format) << "\n";
    std::cout << "|   +-- Header: ....................... " << (config.output_header ? "Yes" : "No") << "\n";
    std::cout << "|   `-- Single file: .................. " << (config.output_single_file ? "Yes" : "No") << "\n";
    std::cout << "`-- Miscellaneous\n";
    std::cout << "    `-- Postprocessing: ............... " << PostprocessingToString(config.postprocessing) << "\n";
    std::cout << std::endl;
}

void PrintStatistics(const EdgeList& edge_list, const VertexRange& range) {
    const SInt num_local_edges = edge_list.size();

    SInt sum_edges;
    SInt min_edges;
    SInt max_edges;

    MPI_Reduce(&num_local_edges, &sum_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, MPI_COMM_WORLD);
    MPI_Reduce(&num_local_edges, &min_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, ROOT, MPI_COMM_WORLD);
    MPI_Reduce(&num_local_edges, &max_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, ROOT, MPI_COMM_WORLD);

    SInt local_max_node = range.second;
    SInt max_node;
    MPI_Reduce(&local_max_node, &max_node, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, ROOT, MPI_COMM_WORLD);

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == ROOT) {
        std::cout << "+-- Number of nodes: .................. " << max_node << "\n";
        std::cout << "`-- Number of edges\n";
        std::cout << "    +-- Total: ........................ " << sum_edges << "\n";
        std::cout << "    +-- Minimum: ...................... " << min_edges << "\n";
        std::cout << "    +-- Maximum: ...................... " << max_edges << "\n";
        std::cout << "    `-- Average: ...................... " << 1.0 * sum_edges / size << std::endl;
    }
}

template <typename Generator>
void RunGenerator(
    PGeneratorConfig& config, const PEID rank, const PEID size, Statistics& stats, Statistics& edge_stats,
    Statistics& edges) {
    // Start timers
    Timer  t;
    double local_time = 0.0;
    double total_time = 0.0;
    t.Restart();

    if (rank == ROOT) {
        std::cout << "Generating graph ..." << std::endl;
    }

    // Chunk distribution
    Generator gen(config, rank, size);
    gen.Generate();

    // Output
    local_time = t.Elapsed();
    MPI_Reduce(&local_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
    if (rank == ROOT) {
        stats.Push(total_time);
        edge_stats.Push(total_time / gen.IO().NumEdges());
        edges.Push(gen.IO().NumEdges());
    }
    PrintStatistics(gen.IO().GetEdges(), gen.GetVertexRange());
    if (rank == ROOT) {
        std::cout << std::endl;
    }

    if (config.postprocessing != Postprocessing::SKIP || config.validate_undirected_graph) {
        if (rank == ROOT) {
            std::cout << "Postprocessing ..." << std::endl;
        }
        Postprocess(config.postprocessing, gen);
        PrintStatistics(gen.IO().GetEdges(), gen.GetVertexRange());
        if (config.validate_undirected_graph) {
            Postprocess(Postprocessing::VALIDATE_UNDIRECTED, gen);
        }
        if (rank == ROOT) {
            std::cout << std::endl;
        }
    }

    if (rank == ROOT) {
        std::cout << "Writing edges ..." << std::endl;
    }
    gen.IO().OutputEdges();
    if (rank == ROOT) {
        std::cout << std::endl;
    }
}

int main(int argn, char** argv) {
    // Init MPI
    MPI_Init(&argn, &argv);
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Read command-line args
    PGeneratorConfig generator_config;
    ParseParameters(argn, argv, rank, size, generator_config);

    if (rank == ROOT)
        OutputParameters(generator_config, rank, size);

    // Statistics
    Statistics stats;
    Statistics edge_stats;
    Statistics edges;

    ULONG user_seed = generator_config.seed;
    for (ULONG i = 0; i < generator_config.iterations; ++i) {
        MPI_Barrier(MPI_COMM_WORLD);
        generator_config.seed = user_seed + i;

        switch (generator_config.generator) {
            case Generator::GNM_DIRECTED:
                RunGenerator<GNMDirected>(generator_config, rank, size, stats, edge_stats, edges);
                break;

            case Generator::GNM_UNDIRECTED:
                RunGenerator<GNMUndirected>(generator_config, rank, size, stats, edge_stats, edges);
                break;

            case Generator::GNP_DIRECTED:
                RunGenerator<GNPDirected>(generator_config, rank, size, stats, edge_stats, edges);
                break;

            case Generator::GNP_UNDIRECTED:
                RunGenerator<GNPUndirected>(generator_config, rank, size, stats, edge_stats, edges);
                break;

            case Generator::RGG_2D:
                RunGenerator<RGG2D>(generator_config, rank, size, stats, edge_stats, edges);
                break;

            case Generator::RGG_3D:
                RunGenerator<RGG3D>(generator_config, rank, size, stats, edge_stats, edges);
                break;

            case Generator::RDG_2D:
#ifdef KAGEN_CGAL_FOUND
                RunGenerator<Delaunay2D>(generator_config, rank, size, stats, edge_stats, edges);
                break;
#else  // KAGEN_CGAL_FOUND
                std::cerr << "Error: RDG 2D generator has been disabled\n";
                std::exit(1);
#endif // KAGEN_CGAL_FOUND

            case Generator::RDG_3D:
#ifdef KAGEN_CGAL_FOUND
                RunGenerator<Delaunay3D>(generator_config, rank, size, stats, edge_stats, edges);
                break;
#else  // KAGEN_CGAL_FOUND
                std::cerr << "Error: RDG 3D genertor has been disabled\n";
                std::exit(1);
#endif // KAGEN_CGAL_FOUND

            case Generator::RHG:
                RunGenerator<Hyperbolic>(generator_config, rank, size, stats, edge_stats, edges);
                break;

            case Generator::BA:
                RunGenerator<Barabassi>(generator_config, rank, size, stats, edge_stats, edges);
                break;

            case Generator::KRONECKER:
                RunGenerator<Kronecker>(generator_config, rank, size, stats, edge_stats, edges);
                break;

            case Generator::GRID_2D:
                RunGenerator<Grid2D>(generator_config, rank, size, stats, edge_stats, edges);
                break;

            case Generator::GRID_3D:
                RunGenerator<Grid3D>(generator_config, rank, size, stats, edge_stats, edges);
                break;

            case Generator::UNDEFINED:
                __builtin_unreachable();
        }
    }

    if (rank == ROOT) {
        std::cout << "RESULT runner=" << GeneratorToString(generator_config.generator) << " time=" << stats.Avg()
                  << " stddev=" << stats.Stddev() << " iterations=" << generator_config.iterations
                  << " edges=" << edges.Avg() << " time_per_edge=" << edge_stats.Avg() << std::endl;
    }

    MPI_Finalize();
    return 0;
}
