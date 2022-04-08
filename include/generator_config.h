/*******************************************************************************
 * include/generator_config.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <iostream>
#include <string>

#include "definitions.h"

namespace kagen {
enum class OutputFormat {
    EDGE_LIST,
    BINARY_EDGE_LIST,
};

enum class Generator {
    GNM_DIRECTED,
    GNM_UNDIRECTED,
    GNP_DIRECTED,
    GNP_UNDIRECTED,
    RGG_2D,
    RGG_3D,
    RDG_2D,
    RDG_3D,
    GRID_2D,
    GRID_3D,
    BA,
    KRONECKER,
    RHG,
    UNDEFINED
};

inline const char* GeneratorToString(const Generator gen) {
    switch (gen) {
        case Generator::GNM_DIRECTED:
            return "gnm_directed";
        case Generator::GNM_UNDIRECTED:
            return "gnm_undirected";
        case Generator::GNP_DIRECTED:
            return "gnp_directed";
        case Generator::GNP_UNDIRECTED:
            return "gnp_undirected";
        case Generator::RGG_2D:
            return "rgg_2d";
        case Generator::RGG_3D:
            return "rgg_3d";
        case Generator::RDG_2D:
            return "rdg_2d";
        case Generator::RDG_3D:
            return "rdg_3d";
        case Generator::GRID_2D:
            return "grid_2d";
        case Generator::GRID_3D:
            return "grid_3d";
        case Generator::BA:
            return "ba";
        case Generator::KRONECKER:
            return "kronecker";
        case Generator::RHG:
            return "rhg";
        default:
            return "undefined";
    }
}

inline Generator StringToGenerator(const std::string& name) {
    for (const Generator gen:
         {Generator::GNM_DIRECTED, Generator::GNM_UNDIRECTED, Generator::GNP_DIRECTED, Generator::GNP_UNDIRECTED,
          Generator::RGG_2D, Generator::RGG_3D, Generator::RDG_2D, Generator::RDG_3D, Generator::GRID_2D,
          Generator::GRID_3D, Generator::BA, Generator::KRONECKER, Generator::RHG}) {
        if (name == GeneratorToString(gen)) {
            return gen;
        }
    }

    return Generator::UNDEFINED;
}

// Configuration for the generator.
struct PGeneratorConfig {
    PGeneratorConfig() {}

    // Type of generator
    Generator generator;
    // Seed for the PRNG
    int seed;
    // Number of nodes/edges
    ULONG n, m;
    // Chunk size
    ULONG k;
    // Edge probability
    double p;
    // Edge radius
    double r;
    // Output filename
    std::string output_file;
    // Output format (text, binary)
    OutputFormat output_format;
    // Output file header
    bool output_header;
    // If set, write graph to a single file
    bool output_single_file;
    // Debug output
    std::string debug_output;
    // Use hash tryagain sampling
    bool hash_sample;
    // Allow self loops
    bool self_loops;
    // Power-law exponent
    double plexp;
    // Avg. degree
    double avg_degree;
    // RHG clique threshold
    double thres;
    // RHG query strategy
    double query_both;
    // BA minimum degree
    double min_degree;
    // Size of histogramm
    ULONG dist_size;
    // Use binomial approximation to hypergeometric
    bool use_binom;
    // Grid dimensions
    ULONG grid_x, grid_y, grid_z;
    // Use periodic boundary condition for grid generators
    bool periodic;
    // Floating-point precision
    ULONG precision;
    // Sampler base size
    ULONG base_size;
    ULONG hyp_base;
    // Benchmarks
    ULONG iterations;
};
} // namespace kagen
