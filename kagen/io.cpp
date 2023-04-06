#include "kagen/io.h"

#include <algorithm>
#include <cassert>
#include <memory>
#include <sstream>

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/io/buffered_writer.h"
#include "kagen/io/coordinates.h"
#include "kagen/io/dot.h"
#include "kagen/io/edgelist.h"
#include "kagen/io/hmetis.h"
#include "kagen/io/metis.h"
#include "kagen/io/parhip.h"
#include "kagen/tools/statistics.h"

namespace kagen {
const std::unordered_map<FileFormat, std::unique_ptr<FileFormatFactory>>& GetGraphFormatFactories() {
    static std::unordered_map<FileFormat, std::unique_ptr<FileFormatFactory>> factories;
    if (factories.empty()) {
        factories[FileFormat::NOOP]                        = std::make_unique<NoopFactory>();
        factories[FileFormat::EDGE_LIST]                   = std::make_unique<EdgelistFactory>();
        factories[FileFormat::EDGE_LIST_UNDIRECTED]        = std::make_unique<UndirectedEdgelistFactory>();
        factories[FileFormat::BINARY_EDGE_LIST]            = std::make_unique<BinaryEdgelistFactory>();
        factories[FileFormat::BINARY_EDGE_LIST_UNDIRECTED] = std::make_unique<UndirectedBinaryEdgelistFactory>();
        factories[FileFormat::XTRAPULP]                    = std::make_unique<XtrapulpFactory>();
        factories[FileFormat::METIS]                       = std::make_unique<MetisFactory>();
        factories[FileFormat::HMETIS]                      = std::make_unique<HmetisFactory>();
        factories[FileFormat::HMETIS_DIRECTED]             = std::make_unique<DirectedHmetisFactory>();
        factories[FileFormat::DOT]                         = std::make_unique<DotFactory>();
        factories[FileFormat::DOT_DIRECTED]                = std::make_unique<DirectedDotFactory>();
        factories[FileFormat::COORDINATES]                 = std::make_unique<CoordinatesFactory>();
        factories[FileFormat::PARHIP]                      = std::make_unique<ParhipFactory>();
    }
    return factories;
}

const std::unique_ptr<FileFormatFactory>& GetGraphFormatFactory(const FileFormat format) {
    const auto& factories = GetGraphFormatFactories();
    const auto  it        = factories.find(format);
    if (it != factories.end()) {
        return (*it).second;
    }

    std::stringstream error_msg;
    error_msg << "file format " << format << " not available for writing";
    throw IOError(error_msg.str());
}

namespace {
std::string GetExtension(const std::string& filename) {
    const auto last_dot_pos = filename.find_last_of('.');
    if (last_dot_pos != std::string::npos) {
        return filename.substr(last_dot_pos + 1);
    }
    return filename;
}
} // namespace

std::unique_ptr<GraphReader> CreateGraphReader(const std::string& filename, const InputGraphConfig& config) {
    const std::string extension = GetExtension(filename);

    const auto& factories = GetGraphFormatFactories();
    for (const auto& [format, factory]: factories) {
        if (factory->DefaultExtension() == extension) {
            auto reader = factory->CreateReader(config);
            if (reader != nullptr) {
                return reader;
            }
        }
    }

    std::stringstream error_msg;
    error_msg << "no file format found for filename " << filename << " with file extension " << extension;
    throw IOError(error_msg.str());
}

std::unique_ptr<GraphReader> CreateGraphReader(const FileFormat format, const InputGraphConfig& config) {
    if (format == FileFormat::EXTENSION) {
        return CreateGraphReader(config.filename, config);
    }

    const auto& factories = GetGraphFormatFactories();
    const auto  it        = factories.find(format);
    if (it != factories.end()) {
        auto reader = (*it).second->CreateReader(config);
        if (reader != nullptr) {
            return reader;
        }
    }

    std::stringstream error_msg;
    error_msg << "file format " << format << " not available for reading";
    throw IOError(error_msg.str());
}
} // namespace kagen
