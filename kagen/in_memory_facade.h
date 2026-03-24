#pragma once

#include "kagen/comm/comm.h"
#include "kagen/context.h"
#include "kagen/kagen.h"

namespace kagen {
void  GenerateInMemoryToDisk(PGeneratorConfig config, Comm& comm);
Graph GenerateInMemory(const PGeneratorConfig& config, GraphRepresentation representation, Comm& comm);
} // namespace kagen
