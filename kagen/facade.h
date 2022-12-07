#pragma once

#include <memory>

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/generator.h"

namespace kagen {
std::unique_ptr<GeneratorFactory> CreateGeneratorFactory(GeneratorType type);

Graph Generate(const PGeneratorConfig& config, MPI_Comm comm);
} // namespace kagen
