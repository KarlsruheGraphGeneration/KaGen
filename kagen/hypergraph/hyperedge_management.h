#pragma once

#include "kagen/kagen.h"

using SInt = unsigned long long;

namespace kagen {

std::pair<std::vector<SInt>, std::vector<PinRange>>
NormalizeCurrentHyperedge(std::vector<SInt> pins, std::vector<PinRange> ranges, SInt min_run_length = 4);
}
