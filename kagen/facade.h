#pragma once

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/kagen.h"

#include <mpi.h>

#include <memory>

namespace kagen {
std::unique_ptr<GeneratorFactory> CreateGeneratorFactory(GeneratorType type);

Graph Generate(const PGeneratorConfig& config, GraphRepresentation representation, MPI_Comm comm);
} // namespace kagen
