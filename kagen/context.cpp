#include "kagen/context.h"

#include <iomanip>
#include <ios>
#include <sstream>
#include <unordered_set>
#include <mpi.h>

namespace kagen {
std::unordered_map<std::string, FileFormat> GetOutputFormatMap() {
    return {
        {"noop", FileFormat::NOOP},
        {"edgelist", FileFormat::EDGE_LIST},
        {"edgelist-undirected", FileFormat::EDGE_LIST_UNDIRECTED},
        {"binary-edgelist", FileFormat::BINARY_EDGE_LIST},
        {"binary-edgelist-undirected", FileFormat::BINARY_EDGE_LIST_UNDIRECTED},
        {"plain-edgelist", FileFormat::PLAIN_EDGE_LIST},
        {"metis", FileFormat::METIS},
        {"hmetis", FileFormat::HMETIS},
        {"dot", FileFormat::DOT},
        {"dot-directed", FileFormat::DOT_DIRECTED},
        {"coordinates", FileFormat::COORDINATES},
        {"parhip", FileFormat::PARHIP},
        {"xtrapulp", FileFormat::XTRAPULP},

        {"none", FileFormat::NOOP},                                               // @deprecated
        {"edge-list", FileFormat::EDGE_LIST},                                     // @deprecated
        {"binary-parhip", FileFormat::PARHIP},                                    // @deprecated
        {"edge-list-undirected", FileFormat::EDGE_LIST_UNDIRECTED},               // @deprecated
        {"binary-edge-list", FileFormat::BINARY_EDGE_LIST},                       // @deprecated
        {"binary-edge-list-undirected", FileFormat::BINARY_EDGE_LIST_UNDIRECTED}, // @deprecated
    };
}

std::unordered_map<std::string, FileFormat> GetInputFormatMap() {
    return {
        {"extension", FileFormat::EXTENSION},
        {"metis", FileFormat::METIS},
        {"parhip", FileFormat::PARHIP},
        {"plain-edgelist", FileFormat::PLAIN_EDGE_LIST},
    };
}

std::ostream& operator<<(std::ostream& out, FileFormat output_format) {
    switch (output_format) {
        case FileFormat::NOOP:
            return out << "noop";

        case FileFormat::EXTENSION:
            return out << "extension";

        case FileFormat::EDGE_LIST:
            return out << "edgelist";

        case FileFormat::EDGE_LIST_UNDIRECTED:
            return out << "edgelist-undirected";

        case FileFormat::BINARY_EDGE_LIST:
            return out << "binary-edgelist";

        case FileFormat::BINARY_EDGE_LIST_UNDIRECTED:
            return out << "binary-edgelist-undirected";

        case FileFormat::PLAIN_EDGE_LIST:
            return out << "plain-edgelist";

        case FileFormat::METIS:
            return out << "metis";

        case FileFormat::HMETIS:
            return out << "hmetis";

        case FileFormat::HMETIS_DIRECTED:
            return out << "hmetis-directed";

        case FileFormat::DOT:
            return out << "dot";

        case FileFormat::DOT_DIRECTED:
            return out << "dot-directed";

        case FileFormat::COORDINATES:
            return out << "coordinates";

        case FileFormat::PARHIP:
            return out << "parhip";

        case FileFormat::XTRAPULP:
            return out << "xtrapulp";
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
        {"path", GeneratorType::PATH_DIRECTED},
        {"ba", GeneratorType::BA},
        {"kronecker", GeneratorType::KRONECKER},
        {"rhg", GeneratorType::RHG},
        {"rmat", GeneratorType::RMAT},
        {"image", GeneratorType::IMAGE_MESH},
        {"imagemesh", GeneratorType::IMAGE_MESH},
        {"image-mesh", GeneratorType::IMAGE_MESH},
        {"file", GeneratorType::FILE},
        {"static", GeneratorType::FILE}, // @deprecated
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

        case GeneratorType::RDG_2D:
            return out << "rdg2d";

        case GeneratorType::RDG_3D:
            return out << "rdg3d";

        case GeneratorType::GRID_2D:
            return out << "grid2d";

        case GeneratorType::GRID_3D:
            return out << "grid3d";

        case GeneratorType::PATH_DIRECTED:
            return out << "path-directed";

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

        case GeneratorType::FILE:
            return out << "file";
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
        {"l2", ImageMeshWeightModel::L2},          {"inv-l2", ImageMeshWeightModel::INV_L2},
        {"ratio", ImageMeshWeightModel::RATIO},    {"inv-ratio", ImageMeshWeightModel::INV_RATIO},
        {"sim", ImageMeshWeightModel::SIMILARITY}, {"similarity", ImageMeshWeightModel::SIMILARITY},
    };
}

std::ostream& operator<<(std::ostream& out, ImageMeshWeightModel weight_model) {
    switch (weight_model) {
        case ImageMeshWeightModel::L2:
            return out << "l2";
        case ImageMeshWeightModel::INV_L2:
            return out << "inv-l2";
        case ImageMeshWeightModel::RATIO:
            return out << "ratio";
        case ImageMeshWeightModel::INV_RATIO:
            return out << "inv-ratio";
        case ImageMeshWeightModel::SIMILARITY:
            return out << "similarity";
    }

    return out << "<invalid>";
}

std::unordered_map<std::string, GraphDistribution> GetGraphDistributionMap() {
    return {
        {"balance-vertices", GraphDistribution::BALANCE_VERTICES},
        {"balance-edges", GraphDistribution::BALANCE_EDGES},
    };
}

std::ostream& operator<<(std::ostream& out, GraphDistribution distribution) {
    switch (distribution) {
        case GraphDistribution::BALANCE_VERTICES:
            return out << "balance-vertices";
        case GraphDistribution::BALANCE_EDGES:
            return out << "balance-edges";
    }

    return out << "<invalid>";
}

std::unordered_map<std::string, GraphRepresentation> GetGraphRepresentationMap() {
    return {
        {"edge-list", GraphRepresentation::EDGE_LIST},
        {"csr", GraphRepresentation::CSR},
    };
}

std::ostream& operator<<(std::ostream& out, const GraphRepresentation representation) {
    switch (representation) {
        case GraphRepresentation::EDGE_LIST:
            return out << "edge-list";

        case GraphRepresentation::CSR:
            return out << "csr";
    }

    return out << "<invalid>";
}

std::ostream& operator<<(std::ostream& out, const PGeneratorConfig& config) {
    out << "General Parameters:\n";
    out << "  Seed:                               " << config.seed << "\n";
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

        case GeneratorType::RDG_2D:
        case GeneratorType::RDG_3D:
            out << "  Number of vertices:                 " << config.n << "\n";
            out << "  Periodic boundary condition:        " << (config.periodic ? "yes" : "no") << "\n";
            break;

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

        case GeneratorType::PATH_DIRECTED:
            out << "  Number of vertices:                 " << config.n << "\n";
            out << "  Periodic boundary condition:        " << (config.periodic ? "yes" : "no") << "\n";
            out << "  Permute vertices:                   " << (config.permute ? "yes" : "no") << "\n";
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
            out << "  Weight model:                       " << config.image_mesh.weight_model << " (x "
                << config.image_mesh.weight_multiplier << ", + " << config.image_mesh.weight_offset << ")\n";
            out << "  Weight threshold:                   " << config.image_mesh.weight_min_threshold
                << " <= <weight> <= " << config.image_mesh.weight_max_threshold << "\n";
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

        case GeneratorType::FILE:
            out << "  Input file:                         " << config.input_graph.filename << "\n";
            out << "  File format:                        " << config.input_graph.format << "\n";
            out << "  Distribution:                       " << config.input_graph.distribution << "\n";
            break;
    }

    // RMAT does not use chunks
    if (config.generator != GeneratorType::RMAT && config.generator != GeneratorType::IMAGE_MESH
        && config.generator != GeneratorType::FILE) {
        if (config.k == 0) {
            out << "  Number of chunks:                   auto\n";
        } else {
            out << "  Number of chunks:                   " << config.k << "\n";
        }
    }
    out << "-------------------------------------------------------------------------------\n";

    if (!config.output_graph.formats.empty()) {
        out << "IO Parameters:\n";
        out << "  Filename:                           " << config.output_graph.filename
            << (config.output_graph.extension ? ".*" : "") << "\n";

        out << "  Output formats:                     " << config.output_graph.formats.front();
        for (std::size_t i = 1; i < config.output_graph.formats.size(); ++i) {
            out << ", " << config.output_graph.formats[i];
        }
        out << "\n";
        out << "  Output header:                      " << config.output_graph.header << "\n";
        out << "  Distributed output:                 " << (config.output_graph.distributed ? "yes" : "no") << "\n";
        out << "  Data type width:                    "
            << (config.output_graph.width == 0 ? "format dependent"
                                               : std::to_string(config.output_graph.width) + " bits")
            << "\n";
        out << "-------------------------------------------------------------------------------\n";
    }

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

    // store all option keys
    std::unordered_set<std::string> option_keys;
    for (const auto& [key, value]: options) {
        option_keys.insert(key);
    }

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
                option_keys.erase("type");
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
                option_keys.erase(key);
                return type_it->second;
            }
        }

        throw std::runtime_error("no generator type specified");
    }();

    auto generic_get_or_default = [&](const std::string& option, const auto& default_value, auto&& parser) {
        const auto it = options.find(option);
        if (it == options.end()) {
            return default_value;
        } else {
            option_keys.erase(option);
            return parser(it->second);
        }
    };

    auto get_sint_or_default = [&](const std::string& option, const SInt default_value = 0) {
        return generic_get_or_default(
            option, default_value, [](const std::string& value) { return std::stoull(value); });
    };

    auto get_hpfloat_or_default = [&](const std::string& option, const HPFloat default_value = 0.0) {
        return generic_get_or_default(
            option, default_value, [](const std::string& value) { return std::stold(value); });
    };

    auto get_bool_or_default = [&](const std::string& option, const bool default_value = false) {
        return generic_get_or_default(option, default_value, is_true_value);
    };

    auto get_string_or_default = [&](const std::string& option, const std::string& default_value = "") {
        return generic_get_or_default(option, default_value, [](const std::string& value) { return value; });
    };

    config.generator   = type;
    config.n           = get_sint_or_default("n", 1ull << get_sint_or_default("N"));
    config.m           = get_sint_or_default("m", 1ull << get_sint_or_default("M"));
    config.k           = get_sint_or_default("k");
    config.p           = get_hpfloat_or_default("prob", get_hpfloat_or_default("p"));
    config.r           = get_hpfloat_or_default("radius", get_hpfloat_or_default("r"));
    config.plexp       = get_hpfloat_or_default("gamma", get_hpfloat_or_default("g"));
    config.periodic    = get_bool_or_default("periodic");
    config.avg_degree  = get_hpfloat_or_default("avg_degree", get_hpfloat_or_default("d"));
    config.min_degree  = get_sint_or_default("min_degree", get_sint_or_default("d"));
    config.grid_x      = get_sint_or_default("grid_x", get_sint_or_default("x"));
    config.grid_y      = get_sint_or_default("grid_y", get_sint_or_default("y"));
    config.grid_z      = get_sint_or_default("grid_z", get_sint_or_default("z"));
    config.rmat_a      = get_hpfloat_or_default("rmat_a", get_hpfloat_or_default("a"));
    config.rmat_b      = get_hpfloat_or_default("rmat_b", get_hpfloat_or_default("b"));
    config.rmat_c      = get_hpfloat_or_default("rmat_c", get_hpfloat_or_default("c"));
    config.coordinates = get_bool_or_default("coordinates");
    config.permute     = get_bool_or_default("permute");

    if (config.generator == GeneratorType::IMAGE_MESH) {
        const std::string filename = get_string_or_default("filename");
        if (filename.empty()) {
            throw std::runtime_error("missing filename");
        }
        config.image_mesh.filename = filename;

        config.image_mesh.weight_multiplier =
            get_hpfloat_or_default("weight_multiplier", config.image_mesh.weight_multiplier);
        config.image_mesh.weight_offset = get_hpfloat_or_default("weight_offset", config.image_mesh.weight_offset);
        config.image_mesh.weight_min_threshold =
            get_hpfloat_or_default("min_weight_threshold", config.image_mesh.weight_min_threshold);
        config.image_mesh.weight_max_threshold =
            get_hpfloat_or_default("max_weight_threshold", config.image_mesh.weight_max_threshold);
        config.image_mesh.neighborhood = get_sint_or_default("neighborhood", config.image_mesh.neighborhood);
        config.image_mesh.max_grid_x   = get_sint_or_default("max_grid_x", config.image_mesh.max_grid_x);
        config.image_mesh.max_grid_y   = get_sint_or_default("max_grid_y", config.image_mesh.max_grid_y);
        config.image_mesh.grid_x       = get_sint_or_default("grid_x", config.image_mesh.grid_x);
        config.image_mesh.grid_y       = get_sint_or_default("grid_y", config.image_mesh.grid_y);
        config.image_mesh.cols_per_pe  = get_sint_or_default("cols_per_pe", config.image_mesh.cols_per_pe);
        config.image_mesh.rows_per_pe  = get_sint_or_default("rows_per_pe", config.image_mesh.rows_per_pe);

        const auto        weight_models = GetImageMeshWeightModelMap();
        const std::string weight_model_name =
            get_string_or_default("weight_model", StringifyEnum(config.image_mesh.weight_model));
        const auto weight_model_it = weight_models.find(weight_model_name);
        if (weight_model_it == weight_models.end()) {
            throw std::runtime_error("invalid weight model name");
        }
        config.image_mesh.weight_model = weight_model_it->second;
    } else if (config.generator == GeneratorType::FILE) {
        const std::string filename = get_string_or_default("filename");
        if (filename.empty()) {
            throw std::runtime_error("missing filename");
        }
        config.input_graph.filename = filename;

        const auto        distributions = GetGraphDistributionMap();
        const std::string distribution_name =
            get_string_or_default("distribution", StringifyEnum(config.input_graph.distribution));
        const auto distribution_it = distributions.find(distribution_name);
        if (distribution_it == distributions.end()) {
            throw std::runtime_error("invalid graph distribution");
        }
        config.input_graph.distribution = distribution_it->second;

        const auto        formats     = GetInputFormatMap();
        const std::string format_name = get_string_or_default("input_format", StringifyEnum(config.input_graph.format));
        const auto        format_it   = formats.find(format_name);
        if (format_it == formats.end()) {
            throw std::runtime_error("invalid graph input format");
        }
        config.input_graph.format = format_it->second;
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == ROOT && !option_keys.empty()) {
        std::stringstream sstr;
        sstr << "WARNING: unused options: ";
        for (const auto& key: option_keys) {
            sstr << key << " ";
        }
        sstr << "\n";
        std::cout << sstr.str();
    }

    return config;
}
} // namespace kagen
