/*******************************************************************************
 * include/generator_config.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include "kagen/kagen.h"

#include <iostream>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>
#include <unordered_map>

namespace kagen {
enum class OutputHeader {
    ALWAYS,
    ROOT,
    NEVER,
};

std::unordered_map<std::string, OutputHeader> GetOutputHeaderMap();

std::ostream& operator<<(std::ostream& out, OutputHeader output_header);

struct ImageMeshConfig {
    std::string          filename             = "";
    ImageMeshWeightModel weight_model         = ImageMeshWeightModel::RATIO;
    double               weight_multiplier    = 255.0;
    double               weight_offset        = 0.0;
    double               weight_min_threshold = 1.0;
    double               weight_max_threshold = std::numeric_limits<double>::max();
    SInt                 neighborhood         = 4;
    SInt                 max_grid_x           = 0;
    SInt                 max_grid_y           = 0;
    SInt                 grid_x               = 0;
    SInt                 grid_y               = 0;
    SInt                 cols_per_pe          = 0;
    SInt                 rows_per_pe          = 0;
};

struct InputGraphConfig {
    std::string       filename     = "";
    FileFormat        format       = FileFormat::EXTENSION;
    GraphDistribution distribution = GraphDistribution::BALANCE_VERTICES;
    int               width        = 64;
};

struct OutputGraphConfig {
    std::string             filename    = "out";
    bool                    extension   = false;
    std::vector<FileFormat> formats     = {FileFormat::EDGE_LIST};
    OutputHeader            header      = OutputHeader::ROOT;
    bool                    distributed = false;
    int                     width       = 64;
};

// Configuration for the generator.
struct PGeneratorConfig {
    // General settings
    bool            quiet                 = false; // Disable all console output
    bool            validate_simple_graph = false; // Validate that the result is a simple graph
    StatisticsLevel statistics_level      = StatisticsLevel::BASIC;
    bool            skip_postprocessing   = false;
    bool            print_header          = true;

    // Generator settings
    GeneratorType generator;          // Generator type
    SInt          n          = 0;     // Number of nodes
    SInt          m          = 0;     // Number of edges
    SInt          k          = 0;     // Number of chunks
    double        p          = 0.0;   // Edge probability
    double        r          = 0.0;   // Edge radius
    bool          self_loops = false; // Allow self loops
    double        plexp      = 2.6;   // Power law exponent
    double        avg_degree = 0.0;   // Average degree
    double        thres      = 0.0;   // Clique threshold (RHG)
    bool          query_both = false; // Query strategy (RHG) -- should be set to false
    SInt          min_degree = 0.0;   // Minimum degree (BA)
    SInt          grid_x     = 0;     // Grid x dimension (Grid2D, Grid3D)
    SInt          grid_y     = 0;     // Grid y dimension (Grid2D, Grid3D)
    SInt          grid_z     = 0;     // Grid z dimension (Grid3D)
    bool          periodic   = false; // Use periodic boundary (Grid2D, Grid3D)
    int           hp_floats  = 0;     // Use 80 bit floating point numbers for RHG generator, 0 for auto
    double        rmat_a     = 0.0;
    double        rmat_b     = 0.0;
    double        rmat_c     = 0.0;
    bool          directed   = false;
    bool          permute    = false; // Permute node vertices

    double max_vertex_imbalance = 0.1; // RGG, RDG, RHG

    bool coordinates = false; // Store vertex coordinates

    // Image mesh generator settings
    ImageMeshConfig image_mesh{};

    // Settings for the static graph pseudo-generator
    InputGraphConfig input_graph{};

    // Hashing / sampling settings
    int  seed        = 1;      // Seed for PRNG
    bool hash_sample = false;  // Use hash tryagain sampling
    bool use_binom   = false;  // Use binomial approximation to hypergeomtry
    SInt precision   = 32;     // Floating-point precision
    SInt base_size   = 1 << 8; // Sampler base size
    SInt hyp_base    = 1 << 8;

    OutputGraphConfig output_graph{};
};

std::ostream& operator<<(std::ostream& out, const PGeneratorConfig& config);

PGeneratorConfig CreateConfigFromString(const std::string& options_str, PGeneratorConfig config = {});

template <typename Enum>
std::string StringifyEnum(Enum value) {
    std::stringstream ss;
    ss << value;
    return ss.str();
}
} // namespace kagen
