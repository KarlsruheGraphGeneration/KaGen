#pragma once

#include "kagen/kagen.h"

namespace kagen {
constexpr PEID        ROOT             = 0;
constexpr std::size_t NEWTON_MAX_ITERS = 10000;
constexpr double      NEWTON_EPS       = 0.001;

enum Direction {
    Up,
    Down,
    Left,
    Right,
    Front,
    Back,
};
} // namespace kagen
