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

#include "benchmark.h"
#include "generator_config.h"
#include "io/generator_io.h"
#include "parse_parameters.h"
#include "postprocessing.h"
#include "timer.h"

#if KAGEN_CGAL_FOUND
    #include "geometric/delaunay/delaunay_2d.h"
    #include "geometric/delaunay/delaunay_3d.h"
#endif // KAGEN_CGAL_FOUND

#include "barabassi/barabassi.h"
#include "geometric/rgg/rgg_2d.h"
#include "geometric/rgg/rgg_3d.h"
#include "gnm/gnm_directed.h"
#include "gnm/gnm_undirected.h"
#include "gnp/gnp_directed.h"
#include "gnp/gnp_undirected.h"
#include "grid/grid_2d.h"
#include "grid/grid_3d.h"
#include "hyperbolic/hyperbolic.h"
#include "kronecker/kronecker.h"

using namespace kagen;

void OutputParameters(PGeneratorConfig& config, const PEID /* rank */, const PEID size) {
    switch (config.generator) {
        case Generator::GNM_DIRECTED:
        case Generator::GNM_UNDIRECTED:
        case Generator::GNP_DIRECTED:
        case Generator::GNP_UNDIRECTED:
            std::cout << "generate graph (n=" << config.n << ", m=" << config.m << " (p=" << config.p
                      << "), k=" << config.k << ", s=" << config.seed << ", P=" << size << ")" << std::endl;
            break;

        case Generator::RGG_2D:
        case Generator::RGG_3D:
            std::cout << "generate graph (n=" << config.n << ", r=" << config.r << ", k=" << config.k
                      << ", s=" << config.seed << ", P=" << size << ")" << std::endl;
            break;

        case Generator::RDG_2D:
        case Generator::RDG_3D:
            std::cout << "generate graph (n=" << config.n << ", k=" << config.k << ", s=" << config.seed
                      << ", P=" << size << ")" << std::endl;
            break;

        case Generator::RHG:
            std::cout << "generate graph (n=" << config.n << ", d=" << config.avg_degree << ", gamma=" << config.plexp
                      << ", k=" << config.k << ", s=" << config.seed << ", P=" << size << ")" << std::endl;
            break;

        case Generator::BA:
            std::cout << "generate graph (n=" << config.n << ", d=" << config.min_degree << ", k=" << config.k
                      << ", s=" << config.seed << ", P=" << size << ")" << std::endl;
            break;

        case Generator::GRID_2D:
            std::cout << "generate graph (row=" << config.grid_x << ", col=" << config.grid_y << ", p=" << config.p
                      << ", k=" << config.k << ", s=" << config.seed << ", P=" << size << ")" << std::endl;
            break;

        case Generator::GRID_3D:
            std::cout << "generate graph (x=" << config.grid_x << ", y=" << config.grid_y << ", z=" << config.grid_z
                      << ", p=" << config.p << ", k=" << config.k << ", s=" << config.seed << ", P=" << size << ")"
                      << std::endl;
            break;

        case Generator::KRONECKER:
            std::cout << "generate graph (?)" << std::endl; // @todo
            break;

        default:
            std::cerr << "Error: unknown generator type\n";
            std::exit(1);
    }

    std::cout << "IO options:" << std::endl;
    std::cout << "- Output filename: " << config.output_file << std::endl;
    std::cout << "- Output format: ";
    switch (config.output_format) {
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
    std::cout << "- Output header: " << config.output_header << std::endl;
    std::cout << "- Output to a single file: " << config.output_single_file << std::endl;
    std::cout << std::endl;
    std::cout << "Other options:" << std::endl;
    std::cout << "- Postprocessing: " << PostprocessingToString(config.postprocessing) << std::endl;
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

    if (config.postprocessing != Postprocessing::SKIP) {
        if (rank == ROOT) {
            std::cout << "postprocessing..." << std::endl;
        }
        Postprocess(config.postprocessing, gen);
    }

    if (rank == ROOT) {
        std::cout << "write output..." << std::endl;
    }
    gen.IO().OutputEdges();
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
