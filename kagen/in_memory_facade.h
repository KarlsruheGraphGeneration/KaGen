#pragma once

#include "kagen/context.h"
#include "kagen/kagen.h"

#include <mpi.h>

namespace kagen {
void  GenerateInMemoryToDisk(PGeneratorConfig config, MPI_Comm comm);
Graph GenerateInMemory(const PGeneratorConfig& config, GraphRepresentation representation, MPI_Comm comm);
} // namespace kagen
