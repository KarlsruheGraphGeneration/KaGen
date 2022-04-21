#pragma once

#include <memory>

#include "kagen/definitions.h"
#include "kagen/context.h"
#include "kagen/generators/generator.h"

namespace kagen {
std::unique_ptr<Generator> CreateGenerator(const PGeneratorConfig& config, PEID rank, PEID size);

std::pair<EdgeList, VertexRange> Generate(const PGeneratorConfig& config);

std::pair<EdgeList, VertexRange> Generate(const PGeneratorConfig& config, PEID rank, PEID size);
} // namespace kagen
