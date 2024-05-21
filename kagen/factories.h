#pragma once

#include "kagen/generators/generator.h"
#include "kagen/kagen.h"

#include <memory>

namespace kagen {
std::unique_ptr<GeneratorFactory> CreateGeneratorFactory(GeneratorType type);
}
