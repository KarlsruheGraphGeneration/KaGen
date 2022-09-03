#include "kagen/generators/geometric/rgg.h"

#include "kagen/tools/newton.h"

namespace kagen {
namespace {
// According to https://mathworld.wolfram.com/SquareLinePicking.html,
// when placing two points at a random position within a unit square, the
// probability that their distance is less than r (where 0 <= r <= 1)
// is given by
//
// prob. edge = 1/2 l^4 - 8/3 l^3 + pi l^2
//
// Note that the case there r > 1 is not implemented.
inline double ComputeEdgeProability2D(const LPFloat r) {
    const double r2 = r * r;
    return r2 * (r2 / 2.0 - 8.0 / 3.0 * r + M_PI);
}

// = d/dx ComputeEdgeProability(x)
inline double ComputeDerivedEdgeProbability2D(const LPFloat r) {
    const double r2 = r * r;
    return (2 * r * (0.5 * r2 - 8.0 / 3.0 * r + M_PI) + r2 * (r - 8.0 / 3.0));
}

// Approximate r such that n(n - 1) * ComputeEdgeProability(r) - m = 0 using a
// naive implementation of Newton's method
double ApproxRadius2D(const SInt n, const SInt m) {
    const HPFloat max_m = 1.0l * n * (n - 1);
    return FindRoot(
        [max_m, m](const double r) { return max_m * ComputeEdgeProability2D(r) - m; },
        [max_m](const double r) { return max_m * ComputeDerivedEdgeProbability2D(r); }, 0.5, NEWTON_EPS,
        NEWTON_MAX_ITERS);
}

SInt ApproxNumNodes2D(const SInt m, const double r) {
    const double p = ComputeEdgeProability2D(r);
    return (p + std::sqrt(p * p + 4 * p * m)) / (2.0 * p);
}

// According to https://mathworld.wolfram.com/CubeLinePicking.html,
// when placing two points at a random position in a unit cube, the probability
// that the distance between these points is equal to l (with 0 <= l <= 1) is
// given by
//
// P(l) = -l^2 [(l - 8) l^2 + pi (6l - 4)]
//
// Thus, the probability for there to be an edge is given by
//
// prob. edge = int_0^r P(l) dl
//            = 1/30 [(48 - 5 r) r^5 - 5 pi r^3 (9r - 8)]
//
// Note that the case where r > 1 is not implemented.
inline double ComputeEdgeProability3D(const LPFloat r) {
    const double r3 = r * r * r;
    const double r5 = r3 * r * r;
    return (1.0 / 30.0) * ((48.0 - 5.0 * r) * r5 - 5 * M_PI * r3 * (9.0 * r - 8.0));
}

// = d/dx ComputeEdgeProability(r)
inline double ComputeDerivedEdgeProbability3D(const LPFloat r) {
    return -1.0 * r * r * ((r - 8.0) * r * r + M_PI * (6 * r - 4));
}

// Approximate r such that n(n - 1) * ComputeEdgeProability(r) - m = 0 using a naive implementation of
// Newton's method
double ApproxRadius3D(const SInt n, const SInt m) {
    const HPFloat max_m = 1.0l * n * (n - 1);
    return FindRoot(
        [max_m, m](const double r) { return max_m * ComputeEdgeProability3D(r) - m; },
        [max_m](const double r) { return max_m * ComputeDerivedEdgeProbability3D(r); }, 0.5, NEWTON_EPS,
        NEWTON_MAX_ITERS);
}

SInt ApproxNumNodes3D(const SInt m, const double r) {
    const double p = ComputeEdgeProability3D(r);
    return (p + std::sqrt(p * p + 4 * p * m)) / (2.0 * p);
}

template <typename ApproxRadius, typename ApproxNumNodes>
PGeneratorConfig NormalizeParametersCommon(
    PGeneratorConfig config, ApproxRadius&& approx_radius, ApproxNumNodes&& approx_num_nodes, const bool output) {
    // We have three parameters, from which two must be given:
    // - Number of nodes
    // - Number of edges
    // - Radius
    // If number of nodes or radius is missing, compute it from the other two values
    if (config.r == 0.0) {
        if (config.n == 0 || config.m == 0) {
            throw ConfigurationError("at least two parameters out of {n, m, r} must be nonzero");
        }

        config.r = approx_radius(config.n, config.m * 2); // x2 because we want undirected edges
        if (output) {
            std::cout << "Setting radius to " << config.r << std::endl;
        }
    } else if (config.n == 0) {
        if (config.r == 0.0 || config.m == 0) {
            throw ConfigurationError("at least two parameters out of {n, m, r} must be nonzero");
        }
        if (config.r > 1.0 && output) {
            std::cerr
                << "Warning: deducing the number of nodes with a radius larger than 1 might give unexpected results\n";
        }

        config.n = approx_num_nodes(config.m * 2, config.r); // x2 because we want undirecte edges
        if (output) {
            std::cout << "Setting number of nodes to " << config.n << std::endl;
        }
    }

    if (config.r < 0) {
        throw ConfigurationError("generator configuration infeasible (negative radius)");
    }

    return config;
}
} // namespace

int RGG2DFactory::Requirements() const {
    return GeneratorRequirement::SQUARE_CHUNKS | GeneratorRequirement::POWER_OF_TWO_COMMUNICATOR_SIZE;
}

PGeneratorConfig RGG2DFactory::NormalizeParameters(PGeneratorConfig config, const bool output) const {
    return NormalizeParametersCommon(config, &ApproxRadius2D, &ApproxNumNodes2D, output);
}

std::unique_ptr<Generator>
RGG2DFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<RGG2D>(config, rank, size);
}

int RGG3DFactory::Requirements() const {
    return GeneratorRequirement::CUBIC_CHUNKS | GeneratorRequirement::POWER_OF_TWO_COMMUNICATOR_SIZE;
}

PGeneratorConfig RGG3DFactory::NormalizeParameters(PGeneratorConfig config, const bool output) const {
    return NormalizeParametersCommon(config, &ApproxRadius3D, &ApproxNumNodes3D, output);
}

std::unique_ptr<Generator>
RGG3DFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<RGG3D>(config, rank, size);
}
} // namespace kagen
