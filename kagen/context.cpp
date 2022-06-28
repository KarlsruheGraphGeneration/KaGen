#include "kagen/context.h"

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

std::ostream& operator<<(std::ostream& out, const PGeneratorConfig& config) {
    out << "General Parameters:\n";
    out << "  Validate simple graph:              " << (config.validate_simple_graph ? "yes" : "no") << "\n";
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
            break;

        case GeneratorType::GNP_DIRECTED:
        case GeneratorType::GNP_UNDIRECTED:
            out << "  Number of vertices:                 " << config.n << "\n";
            out << "  Edge probability:                   " << config.p << "\n";
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
            break;

        case GeneratorType::BA:
            out << "  Number of vertices:                 " << config.n << "\n";
            out << "  Minimum degree:                     " << config.min_degree << "\n";
            break;

        case GeneratorType::KRONECKER:
            out << "  Number of vertices:                 " << config.n << "\n";
            out << "  Number of edges:                    " << config.m << "\n";
            break;

        case GeneratorType::RHG:
            out << "  Power-law exponent:                 " << config.plexp << "\n";
            out << "  Number of vertices:                 " << (config.n == 0 ? "auto" : std::to_string(config.n))
                << "\n";
            out << "  Number of edges:                    " << (config.m == 0 ? "auto" : std::to_string(config.m))
                << "\n";
            out << "  Average degree:                     "
                << (config.avg_degree == 0.0 ? "auto" : std::to_string(config.avg_degree)) << "\n";
            out << "  High-resolution floating points:    " << (config.hi_res ? "yes" : "no") << "\n";
            break;
    }
    if (config.k == 0) {
        out << "  Number of chunks:                   auto\n";
    } else {
        out << "  Number of chunks:                   " << config.k << "\n";
    }
    out << "-------------------------------------------------------------------------------\n";

    out << "IO Parameters:\n";
    out << "  Filename:                           " << config.output_file << "\n";
    out << "  Output format:                      " << config.output_format << "\n";
    out << "  Output header:                      " << config.output_header << "\n";
    out << "  Single file:                        " << (config.output_single_file ? "yes" : "no") << "\n";
    out << "-------------------------------------------------------------------------------\n";

    return out << std::flush;
}
} // namespace kagen
