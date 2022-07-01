#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "hash.hpp"

#include "kagen/generators/barabassi/barabassi.h"
#include "kagen/generators/generator.h"

namespace kagen {
PGeneratorConfig BarabassiFactory::NormalizeParameters(PGeneratorConfig config, bool output) const {
    if (config.min_degree == 0) {
        if (config.n == 0 || config.m == 0) {
            throw ConfigurationError("at least two parameters out of {n, m, d} must be nonzero");
        }

        config.min_degree = 1.0 * config.m / config.n;
        if (output) {
            std::cout << "Setting minimum degree to " << config.min_degree << std::endl;
        }
    } else if (config.n == 0) {
        if (config.m == 0 || config.min_degree == 0) {
            throw ConfigurationError("at least two parameters out of {n, m, d} must be nonzero");
        }

        config.n = 1.0 * config.m / config.min_degree;
        if (output) {
            std::cout << "Setting number of vertices to " << config.n << std::endl;
        }
    }

    return config;
}

std::unique_ptr<Generator>
BarabassiFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<Barabassi>(config, rank, size);
}

Barabassi::Barabassi(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      min_degree_(config_.min_degree),
      total_degree_(2 * config_.min_degree) {
    from_ = rank * std::ceil(config_.n / (LPFloat)size);
    to_   = std::min((SInt)((rank + 1) * std::ceil(config_.n / (LPFloat)size) - 1), config_.n - 1);
}

bool Barabassi::IsAlmostUndirected() const {
    // If we do not want a directed graph, we want an undirected graph
    // Since the generator does not always generate reverse edges across PEs, we are 'almost' undirected
    return !config_.directed;
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
            if (config_.self_loops || v != w) {
                PushEdge(v, w);
                if (v != w && from_ <= w && w <= to_) {
                    PushEdge(w, v);
                }
            }
        }
    }

    FilterDuplicateEdges();
}
} // namespace kagen
