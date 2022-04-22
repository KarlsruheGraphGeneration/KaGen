#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "hash.hpp"
#include "kagen/generators/barabassi/barabassi.h"

namespace kagen {
Barabassi::Barabassi(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      min_degree_(config_.min_degree),
      total_degree_(2 * config_.min_degree) {
    from_ = rank * std::ceil(config_.n / (LPFloat)size);
    to_   = std::min((SInt)((rank + 1) * ceil(config_.n / (LPFloat)size) - 1), config_.n - 1);
}

void Barabassi::GenerateImpl() {
    GenerateEdges();
    SetVertexRange(from_, to_ + 1);
}

void Barabassi::GenerateEdges() {
    for (SInt v = from_; v <= to_; v++) {
        for (SInt i = 0; i < min_degree_; i++) {
            SInt r = 2 * (v * min_degree_ + i) + 1;
            do {
                // compute hash h(r)
                SInt hash = sampling::Spooky::hash(config_.seed + r);
                r         = hash % r;
            } while (r % 2 == 1);
            SInt w = r / total_degree_;
            PushEdge(v, w);
        }
    }
}
} // namespace kagen
