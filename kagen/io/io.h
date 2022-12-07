#pragma once

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/definitions.h"

namespace kagen {
void WriteGraph(const PGeneratorConfig& config, Graph& graph, MPI_Comm comm);
} // namespace kagen
