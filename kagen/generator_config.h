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

inline const char* OutputFormatToString(const OutputFormat format) {
    switch (format) {
        case OutputFormat::EDGE_LIST:
            return "edge_list";

        case OutputFormat::BINARY_EDGE_LIST:
            return "binary_edge_list";

        default:
            return "undefined";
    }
}

inline OutputFormat StringToOutputFormat(const std::string& name) {
    for (OutputFormat format: {OutputFormat::EDGE_LIST, OutputFormat::BINARY_EDGE_LIST}) {
        if (name == OutputFormatToString(format)) {
            return format;
        }
    }

    std::cerr << "Error: invalid output format " << name << "\n";
    std::exit(1);
}

enum class GeneratorType {
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

inline const char* GeneratorTypeToString(const GeneratorType gen) {
    switch (gen) {
        case GeneratorType::GNM_DIRECTED:
            return "gnm_directed";
        case GeneratorType::GNM_UNDIRECTED:
            return "gnm_undirected";
        case GeneratorType::GNP_DIRECTED:
            return "gnp_directed";
        case GeneratorType::GNP_UNDIRECTED:
            return "gnp_undirected";
        case GeneratorType::RGG_2D:
            return "rgg_2d";
        case GeneratorType::RGG_3D:
            return "rgg_3d";
        case GeneratorType::RDG_2D:
            return "rdg_2d";
        case GeneratorType::RDG_3D:
            return "rdg_3d";
        case GeneratorType::GRID_2D:
            return "grid_2d";
        case GeneratorType::GRID_3D:
            return "grid_3d";
        case GeneratorType::BA:
            return "ba";
        case GeneratorType::KRONECKER:
            return "kronecker";
        case GeneratorType::RHG:
            return "rhg";
        default:
            return "undefined";
    }
}

inline GeneratorType StringToGeneratorType(const std::string& name) {
    for (const GeneratorType gen:
         {GeneratorType::GNM_DIRECTED, GeneratorType::GNM_UNDIRECTED, GeneratorType::GNP_DIRECTED,
          GeneratorType::GNP_UNDIRECTED, GeneratorType::RGG_2D, GeneratorType::RGG_3D, GeneratorType::RDG_2D,
          GeneratorType::RDG_3D, GeneratorType::GRID_2D, GeneratorType::GRID_3D, GeneratorType::BA,
          GeneratorType::KRONECKER, GeneratorType::RHG}) {
        if (name == GeneratorTypeToString(gen)) {
            return gen;
        }
    }

    if (name == "rmat") {
        return GeneratorType::KRONECKER;
    }

    return GeneratorType::UNDEFINED;
}

enum class Postprocessing {
    VALIDATE_RANGES,
    VALIDATE_RANGES_CONSECUTIVE,
    VALIDATE_UNDIRECTED,
    FIX_UNDIRECTED_EDGE_LIST,
    REDISTRIBUTE_GRAPH,
    SKIP
};

inline const char* PostprocessingToString(const Postprocessing postprocessing) {
    switch (postprocessing) {
        case Postprocessing::VALIDATE_RANGES:
            return "validate_ranges";

        case Postprocessing::VALIDATE_RANGES_CONSECUTIVE:
            return "validate_ranges_consecutive";

        case Postprocessing::VALIDATE_UNDIRECTED:
            return "validate_undirected";

        case Postprocessing::FIX_UNDIRECTED_EDGE_LIST:
            return "fix_undirected_edge_list";

        case Postprocessing::REDISTRIBUTE_GRAPH:
            return "redistribute";

        case Postprocessing::SKIP:
            return "skip";
    }

    __builtin_unreachable();
}

inline Postprocessing StringToPostprocessing(const std::string& name) {
    for (const Postprocessing postprocessing:
         {Postprocessing::VALIDATE_RANGES, Postprocessing::VALIDATE_RANGES_CONSECUTIVE,
          Postprocessing::VALIDATE_UNDIRECTED, Postprocessing::FIX_UNDIRECTED_EDGE_LIST,
          Postprocessing::REDISTRIBUTE_GRAPH, Postprocessing::SKIP}) {
        if (name == PostprocessingToString(postprocessing)) {
            return postprocessing;
        }
    }

    std::cerr << "Error: invalid postprocessing option " << name << "\n";
    std::exit(1);
}

// Configuration for the generator.
struct PGeneratorConfig {
    PGeneratorConfig() {}

    // Type of generator
    GeneratorType generator;
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
    // Optional postprocessing (validation, make it undirected, fix edge list etc)
    Postprocessing postprocessing;
    bool           validate_undirected_graph;
};
} // namespace kagen
