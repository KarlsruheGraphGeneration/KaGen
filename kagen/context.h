/*******************************************************************************
 * include/generator_config.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include <iostream>
#include <ostream>
#include <string>
#include <unordered_map>

#include "kagen/definitions.h"

namespace kagen {
enum class OutputFormat {
    NONE,
    EDGE_LIST,
    BINARY_EDGE_LIST,
    METIS,
    HMETIS,
};

std::unordered_map<std::string, OutputFormat> GetOutputFormatMap();

std::ostream& operator<<(std::ostream& out, OutputFormat output_format);

enum class OutputHeader {
    ALWAYS,
    ROOT,
    NEVER,
};

std::unordered_map<std::string, OutputHeader> GetOutputHeaderMap();

std::ostream& operator<<(std::ostream& out, OutputHeader output_header);

enum class GeneratorType {
    GNM_DIRECTED,
    GNM_UNDIRECTED,
    GNP_DIRECTED,
    GNP_UNDIRECTED,
    RGG_2D,
    RGG_3D,
#ifdef KAGEN_CGAL_FOUND
    RDG_2D,
    RDG_3D,
#endif // KAGEN_CGAL_FOUND
    GRID_2D,
    GRID_3D,
    BA,
    KRONECKER,
    RHG,
};

std::unordered_map<std::string, GeneratorType> GetGeneratorTypeMap();

std::ostream& operator<<(std::ostream& out, GeneratorType generator_type);

enum class StatisticsLevel : std::uint8_t {
    NONE     = 0,
    BASIC    = 1,
    ADVANCED = 2,
};

bool operator<=(StatisticsLevel a, StatisticsLevel b);

std::unordered_map<std::string, StatisticsLevel> GetStatisticsLevelMap();

std::ostream& operator<<(std::ostream& out, StatisticsLevel statistics_level);

// Configuration for the generator.
struct PGeneratorConfig {
    // General settings
    bool            quiet                 = false; // Disable all console output
    bool            validate_simple_graph = false; // Validate that the result is a simple graph
    StatisticsLevel statistics_level      = StatisticsLevel::BASIC;

    // Generator settings
    GeneratorType generator;          // Generator type
    SInt          n          = 0;     // Number of nodes
    SInt          m          = 0;     // Number of edges
    SInt          k          = 0;     // Number of chunks
    double        p          = 0.5;   // Edge probability
    double        r          = 0.125; // Edge radius
    bool          self_loops = false; // Allow self loops
    double        plexp      = 2.6;   // Power law exponent
    double        avg_degree = 5.0;   // Average degree
    double        thres      = 0.0;   // Clique threshold (RHG)
    bool          query_both = false; // Query strategy (RHG) -- should be set to false
    double        min_degree = 4.0;   // Minimum degree (BA)
    SInt          grid_x     = 0;     // Grid x dimension (Grid2D, Grid3D)
    SInt          grid_y     = 0;     // Grid y dimension (Grid2D, Grid3D)
    SInt          grid_z     = 0;     // Grid z dimension (Grid3D)
    bool          periodic   = false; // Use periodic boundary (Grid2D, Grid3D)

    bool coordinates = false; // Store vertex coordinates

    // Hashing / sampling settings
    int   seed        = 1;      // Seed for PRNG
    bool  hash_sample = false;  // Use hash tryagain sampling
    bool  use_binom   = false;  // Use binomial approximation to hypergeomtry
    ULONG precision   = 32;     // Floating-point precision
    ULONG base_size   = 1 << 8; // Sampler base size
    ULONG hyp_base    = 1 << 8;

    // IO settings
    OutputFormat output_format      = OutputFormat::EDGE_LIST; // Output format
    OutputHeader output_header      = OutputHeader::ROOT;      // PEs that print file headers
    std::string  output_file        = "out";                   // Output filename
    bool         output_single_file = true;                    // Collect all graphs in a single output file
};

std::ostream& operator<<(std::ostream& out, const PGeneratorConfig& config);
} // namespace kagen
