#pragma once

#include "kagen/context.h"
#include "kagen/kagen.h"

#include <mpi.h>

namespace kagen {
void GenerateExternalMemoryToDisk(PGeneratorConfig config, MPI_Comm comm);
}
