#pragma once

#include <cstdlib>

namespace kagen {
template <typename Fn, typename De>
double FindRoot(Fn&& fn, De&& de, const double x0, const double eps, const std::size_t max_iter) {
    double x = x0;
    for (std::size_t iter = 0; iter < max_iter && std::abs(fn(x)) > eps; ++iter) {
        x = x - 1.0 * fn(x) / de(x);
    }
    return x;
}
} // namespace kagen
