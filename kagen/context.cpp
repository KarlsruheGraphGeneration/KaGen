#include "kagen/context.h"

#include <iomanip>
#include <ios>
#include <mpi.h>

namespace kagen {
std::unordered_map<std::string, OutputFormat> GetOutputFormatMap() {
    return {
        {"none", OutputFormat::NONE},
        {"edge-list", OutputFormat::EDGE_LIST},
        {"binary-edge-list", OutputFormat::BINARY_EDGE_LIST},
        {"metis", OutputFormat::METIS},
        {"hmetis", OutputFormat::HMETIS},
        {"dot", OutputFormat::DOT},
        {"dot-directed", OutputFormat::DOT_DIRECTED},
        {"coordinates", OutputFormat::COORDINATES},
        {"binary-parhip", OutputFormat::BINARY_PARHIP},
    };
}

std::ostream& operator<<(std::ostream& out, OutputFormat output_format) {
    switch (output_format) {
        case OutputFormat::NONE:
            return out << "none";

        case OutputFormat::EDGE_LIST:
            return out << "edge-list";

        case OutputFormat::BINARY_EDGE_LIST:
            return out << "binary-edge-list";

        case OutputFormat::METIS:
            return out << "metis";

        case OutputFormat::HMETIS:
            return out << "hmetis";

        case OutputFormat::DOT:
            return out << "dot";

        case OutputFormat::DOT_DIRECTED:
            return out << "dot-directed";

        case OutputFormat::COORDINATES:
            return out << "coordinates";

        case OutputFormat::BINARY_PARHIP:
            return out << "binary-parhip";
    }

    return out << "<invalid>";
}

std::unordered_map<std::string, OutputHeader> GetOutputHeaderMap() {
    return {
        {"always", OutputHeader::ALWAYS},
        {"root", OutputHeader::ROOT},
        {"never", OutputHeader::NEVER},
    };
}

std::ostream& operator<<(std::ostream& out, OutputHeader output_header) {
    switch (output_header) {
        case OutputHeader::ALWAYS:
            return out << "always";

        case OutputHeader::ROOT:
            return out << "root";

        case OutputHeader::NEVER:
            return out << "never";
    }

    return out << "<invalid>";
}

std::unordered_map<std::string, GeneratorType> GetGeneratorTypeMap() {
    return {
        {"gnm-directed", GeneratorType::GNM_DIRECTED},
        {"gnm-undirected", GeneratorType::GNM_UNDIRECTED},
        {"gnp-directed", GeneratorType::GNP_DIRECTED},
        {"gnp-undirected", GeneratorType::GNP_UNDIRECTED},
        {"rgg2d", GeneratorType::RGG_2D},
        {"rgg3d", GeneratorType::RGG_3D},
#ifdef KAGEN_CGAL_FOUND
        {"rdg2d", GeneratorType::RDG_2D},
        {"rdg3d", GeneratorType::RDG_3D},
#endif // KAGEN_CGAL_FOUND
        {"grid2d", GeneratorType::GRID_2D},
        {"grid3d", GeneratorType::GRID_3D},
        {"ba", GeneratorType::BA},
        {"kronecker", GeneratorType::KRONECKER},
        {"rhg", GeneratorType::RHG},
        {"rmat", GeneratorType::RMAT},
        {"image-mesh", GeneratorType::IMAGE_MESH},
    };
}

std::ostream& operator<<(std::ostream& out, GeneratorType generator_type) {
    switch (generator_type) {
        case GeneratorType::GNM_DIRECTED:
            return out << "gnm-directed";

        case GeneratorType::GNM_UNDIRECTED:
            return out << "gnm-undirected";

        case GeneratorType::GNP_DIRECTED:
            return out << "gnp-directed";

        case GeneratorType::GNP_UNDIRECTED:
            return out << "gnp-undirected";

        case GeneratorType::RGG_2D:
            return out << "rgg2d";

        case GeneratorType::RGG_3D:
            return out << "rgg3d";

#ifdef KAGEN_CGAL_FOUND
        case GeneratorType::RDG_2D:
            return out << "rdg2d";

        case GeneratorType::RDG_3D:
            return out << "rdg3d";
#endif // KAGEN_CGAL_FOUND

        case GeneratorType::GRID_2D:
            return out << "grid2d";

        case GeneratorType::GRID_3D:
            return out << "grid3d";

        case GeneratorType::BA:
            return out << "ba";

        case GeneratorType::KRONECKER:
            return out << "kronecker";

        case GeneratorType::RHG:
            return out << "rhg";

        case GeneratorType::RMAT:
            return out << "rmat";

        case GeneratorType::IMAGE_MESH:
            return out << "image-mesh";
    }

    return out << "<invalid>";
}

bool operator<=(StatisticsLevel a, StatisticsLevel b) {
    return static_cast<std::uint8_t>(a) <= static_cast<std::uint8_t>(b);
}

std::unordered_map<std::string, StatisticsLevel> GetStatisticsLevelMap() {
    return {
        {"none", StatisticsLevel::NONE},
        {"basic", StatisticsLevel::BASIC},
        {"advanced", StatisticsLevel::ADVANCED},
    };
}

std::ostream& operator<<(std::ostream& out, StatisticsLevel statistics_level) {
    switch (statistics_level) {
        case StatisticsLevel::NONE:
            return out << "none";

        case StatisticsLevel::BASIC:
            return out << "basic";

        case StatisticsLevel::ADVANCED:
            return out << "advanced";
    }

    return out << "<invalid>";
}

std::unordered_map<std::string, ImageMeshWeightModel> GetImageMeshWeightModelMap() {
    return {
        {"l2", ImageMeshWeightModel::L2},
    };
}

std::ostream& operator<<(std::ostream& out, ImageMeshWeightModel weight_model) {
    switch (weight_model) {
        case ImageMeshWeightModel::L2:
            return out << "l2";
    }
}

std::ostream& operator<<(std::ostream& out, const PGeneratorConfig& config) {
    out << "General Parameters:\n";
    out << "  Validate generated graph:           " << (config.validate_simple_graph ? "yes" : "no") << "\n";
    out << "  Statistics level:                   " << config.statistics_level << "\n";
    out << "  Generate coordinates:               " << (config.coordinates ? "yes" : "no") << "\n";
    out << "-------------------------------------------------------------------------------\n";

    out << "Generator Parameters:\n";
    out << "  Generator:                          " << config.generator << "\n";
    switch (config.generator) {
        case GeneratorType::GNM_DIRECTED:
        case GeneratorType::GNM_UNDIRECTED:
            out << "  Number of vertices:                 " << config.n << "\n";
            out << "  Number of edges:                    " << config.m << "\n";
            out << "  Self loops:                         " << (config.self_loops ? "yes" : "no") << "\n";
            out << "  Graph type:                         "
                << (config.generator == GeneratorType::GNM_DIRECTED ? "directed" : "undirected") << "\n";
            break;

        case GeneratorType::GNP_DIRECTED:
        case GeneratorType::GNP_UNDIRECTED:
            out << "  Number of vertices:                 " << config.n << "\n";
            out << "  Edge probability:                   " << config.p << "\n";
            out << "  Self loops:                         " << (config.self_loops ? "yes" : "no") << "\n";
            out << "  Graph type:                         "
                << (config.generator == GeneratorType::GNP_DIRECTED ? "directed" : "undirected") << "\n";
            break;

        case GeneratorType::RGG_2D:
        case GeneratorType::RGG_3D:
            out << "  Number of vertices:                 " << (config.n == 0 ? "auto" : std::to_string(config.n))
                << "\n";
            out << "  Number of edges:                    " << (config.m == 0 ? "auto" : std::to_string(config.m))
                << "\n";
            out << "  Edge radius:                        " << (config.r == 0.0 ? "auto" : std::to_string(config.r))
                << "\n";
            break;

#ifdef KAGEN_CGAL_FOUND
        case GeneratorType::RDG_2D:
        case GeneratorType::RDG_3D:
            out << "  Number of vertices:                 " << config.n << "\n";
            out << "  Periodic boundary condition:        " << (config.periodic ? "yes" : "no") << "\n";
            break;
#endif // KAGEN_CGAL_FOUND

        case GeneratorType::GRID_3D:
            out << "  Grid z:                             " << config.grid_z << "\n";
            // intentional fall through

        case GeneratorType::GRID_2D:
            out << "  Grid y:                             " << config.grid_y << "\n";
            out << "  Grid x:                             " << config.grid_x << "\n";
            out << "  Number of edges:                    " << (config.m == 0 ? "auto" : std::to_string(config.m))
                << "\n";
            out << "  Edge probability:                   " << (config.p == 0.0 ? "auto" : std::to_string(config.p))
                << "\n";
            out << "  Periodic boundary:                  " << (config.periodic ? "yes" : "no") << "\n";
            break;

        case GeneratorType::BA:
            out << "  Number of vertices:                 " << config.n << "\n";
            out << "  Minimum degree:                     " << config.min_degree << "\n";
            out << "  Self loops:                         " << (config.self_loops ? "yes" : "no") << "\n";
            out << "  Graph type:                         " << (config.directed ? "directed" : "undirected") << "\n";
            break;

        case GeneratorType::KRONECKER:
            out << "  Number of vertices:                 " << config.n << "\n";
            out << "  Number of edges:                    " << config.m << "\n";
            out << "  Self loops:                         " << (config.self_loops ? "yes" : "no") << "\n";
            break;

        case GeneratorType::RHG:
            out << "  Power-law exponent:                 " << config.plexp << "\n";
            out << "  Number of vertices:                 " << (config.n == 0 ? "auto" : std::to_string(config.n))
                << "\n";
            out << "  Number of edges:                    " << (config.m == 0 ? "auto" : std::to_string(config.m))
                << "\n";
            out << "  Average degree:                     "
                << (config.avg_degree == 0.0 ? "auto" : std::to_string(config.avg_degree)) << "\n";
            out << "  High-resolution floating points:    ";
            if (config.hp_floats == 0) {
                out << "auto";
            } else if (config.hp_floats == -1) {
                out << "no";
            } else {
                out << "yes";
            }
            out << "\n";
            break;

        case GeneratorType::RMAT:
            out << "  Number of vertices:                 " << config.n << "\n";
            out << "  Number of edges:                    " << config.m << "\n";
            out << "  Probabilities:                      " << std::setprecision(3) << std::fixed << config.rmat_a
                << " / " << config.rmat_b << " / " << config.rmat_c << " / "
                << 1.0 - config.rmat_a - config.rmat_b - config.rmat_c << "\n";
            out << "  Self loops:                         " << (config.self_loops ? "yes" : "no") << "\n";
            break;

        case GeneratorType::IMAGE_MESH:
            out << "  Input image:                        " << config.image_mesh.filename << "\n";
            out << "  Weight model:                       " << config.image_mesh.weight_model << "\n";
            out << "  Grid size:                          "
                << (config.image_mesh.max_grid_x == 0 ? "auto" : std::to_string(config.image_mesh.max_grid_x)) << " x "
                << (config.image_mesh.max_grid_y == 0 ? "auto" : std::to_string(config.image_mesh.max_grid_y)) << "\n";
            if (config.image_mesh.grid_x != 0 || config.image_mesh.grid_y != 0) {
                out << "    Only use top-left subgrid:        "
                    << (config.image_mesh.grid_x == 0 ? "auto" : std::to_string(config.image_mesh.grid_x)) << " x "
                    << (config.image_mesh.grid_y == 0 ? "auto" : std::to_string(config.image_mesh.grid_y)) << "\n";
            }
            out << "  Subgrid assigned to each PE:        "
                << (config.image_mesh.cols_per_pe == 0 ? "auto" : std::to_string(config.image_mesh.cols_per_pe))
                << " x "
                << (config.image_mesh.rows_per_pe == 0 ? "auto" : std::to_string(config.image_mesh.rows_per_pe))
                << "\n";
            break;
    }

    // RMAT does not use chunks
    if (config.generator != GeneratorType::RMAT && config.generator != GeneratorType::IMAGE_MESH) {
        if (config.k == 0) {
            out << "  Number of chunks:                   auto\n";
        } else {
            out << "  Number of chunks:                   " << config.k << "\n";
        }
    }
    out << "-------------------------------------------------------------------------------\n";

    out << "IO Parameters:\n";
    out << "  Filename:                           " << config.output_file << "\n";
    out << "  Output format:                      " << config.output_format << "\n";
    out << "  Output header:                      " << config.output_header << "\n";
    out << "  Distributed output:                 " << (config.output_single_file ? "no" : "yes") << "\n";
    out << "-------------------------------------------------------------------------------\n";

    return out << std::flush;
}

namespace {
using Options = std::unordered_map<std::string, std::string>;

// Parse a string such as 'key1=value1;key2=value2;key3;key4'
Options ParseOptionString(const std::string& options) {
    std::istringstream toker(options);
    std::string        pair;
    Options            result;

    while (std::getline(toker, pair, ';')) {
        const auto eq_pos = pair.find('=');
        if (eq_pos == std::string::npos) { // No equal sign -> Option is a flag
            result[pair] = "1";
        } else {
            const std::string key   = pair.substr(0, eq_pos);
            const std::string value = pair.substr(eq_pos + 1);
            result[key]             = value;
        }
    }

    return result;
}
} // namespace

PGeneratorConfig CreateConfigFromString(const std::string& options_str, PGeneratorConfig config) {
    const auto options = ParseOptionString(options_str);

    // Used to decide whether a boolean value is true
    auto is_true_value = [](const std::string& value) {
        return value == "1" || value == "true" || value == "yes";
    };

    // Select the generator type, throw an exception if the type not specified
    const GeneratorType type = [&] {
        const auto type_map = GetGeneratorTypeMap();

        // If there is a type key, use it to select the graph model
        const auto type_key_it = options.find("type");
        if (type_key_it != options.end()) {
            const std::string type_name = type_key_it->second;
            const auto        type_it   = type_map.find(type_name);
            if (type_it == type_map.end()) {
                throw std::runtime_error("invalid generator type");
            } else {
                return type_it->second;
            }
        }

        // Otherwise, scan for flags that match a graph model
        for (const auto& [key, value]: options) {
            // Only look for flags / boolean entries
            if (!is_true_value(value)) {
                continue;
            }

            // Check if the key names a graph model
            const auto type_it = type_map.find(key);
            if (type_it != type_map.end()) {
                return type_it->second;
            }
        }

        throw std::runtime_error("no generator type specified");
    }();

    auto get_sint_or_default = [&](const std::string& option, const SInt default_value = 0) {
        const auto it = options.find(option);
        return (it == options.end() ? default_value : std::stoll(it->second));
    };
    auto get_hpfloat_or_default = [&](const std::string& option, const HPFloat default_value = 0.0) {
        const auto it = options.find(option);
        return (it == options.end() ? default_value : std::stod(it->second));
    };
    auto get_bool_or_default = [&](const std::string& option, const bool default_value = false) {
        const auto it = options.find(option);
        return (it == options.end() ? default_value : is_true_value(it->second));
    };

    config.generator   = type;
    config.n           = get_sint_or_default("n", 1ull << get_sint_or_default("N"));
    config.m           = get_sint_or_default("m", 1ull << get_sint_or_default("M"));
    config.k           = get_sint_or_default("k");
    config.p           = get_hpfloat_or_default("prob");
    config.r           = get_hpfloat_or_default("radius");
    config.plexp       = get_hpfloat_or_default("gamma");
    config.periodic    = get_bool_or_default("periodic");
    config.avg_degree  = get_hpfloat_or_default("avg_degree");
    config.min_degree  = get_sint_or_default("min_degree");
    config.grid_x      = get_sint_or_default("grid_x");
    config.grid_y      = get_sint_or_default("grid_y");
    config.grid_z      = get_sint_or_default("grid_z");
    config.rmat_a      = get_sint_or_default("rmat_a");
    config.rmat_b      = get_sint_or_default("rmat_b");
    config.rmat_c      = get_sint_or_default("rmat_c");
    config.coordinates = get_bool_or_default("coordinates");
    return config;
}
} // namespace kagen
