#pragma once

#include "kagen/comm/comm.h"
#include "kagen/context.h"
#include "kagen/kagen.h"

namespace kagen {
void GenerateExternalMemoryToDisk(PGeneratorConfig config, Comm& comm);
}
