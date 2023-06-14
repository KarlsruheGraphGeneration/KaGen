#pragma once

#include <string>
#include <unordered_map>

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/io/graph_format.h"

namespace kagen {
const std::unordered_map<FileFormat, std::unique_ptr<FileFormatFactory>>& GetGraphFormatFactories();

const std::unique_ptr<FileFormatFactory>& GetGraphFormatFactory(FileFormat format);

std::unique_ptr<GraphReader>
CreateGraphReader(const std::string& filename, const InputGraphConfig& config, PEID rank, PEID size);

std::unique_ptr<GraphReader>
CreateGraphReader(const FileFormat format, const InputGraphConfig& config, PEID rank, PEID size);
} // namespace kagen
