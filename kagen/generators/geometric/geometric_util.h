#pragma once

#include "kagen/kagen.h"

namespace kagen {

// n, x_off, y_off, generated, offset
using Chunk = std::tuple<SInt, LPFloat, LPFloat, bool, SInt>;
// n, x_off, y_off, generated, offset
using Cell = std::tuple<SInt, LPFloat, LPFloat, bool, SInt>;
// x, y, id
using Vertex = std::tuple<LPFloat, LPFloat, SInt>;

} // namespace kagen