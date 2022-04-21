#pragma once

#include <string>

#include "kagen/definitions.h"
#include "kagen/generator_config.h"

namespace kagen {
void WriteGraph(const PGeneratorConfig& config, EdgeList& edges, VertexRange vertex_range);

void WriteEdgeList(
    const std::string& filename, bool omit_header, bool single_file, const EdgeList& edges, VertexRange vertex_range);

void WriteBinaryEdgeList(
    const std::string& filename, bool omit_header, bool single_file, const EdgeList& edges, VertexRange vertex_range);

void WriteMetis(const std::string& filename, EdgeList& edges, VertexRange vertex_range);

void WriteHMetis(const std::string& filename, EdgeList& edges, VertexRange vertex_range);
} // namespace kagen
