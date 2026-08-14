
#include "kagen/hypergraph/hypergraph_utils.h"

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/kagen.h"
#include "kagen/sampling/hash.hpp"
#include "kagen/tools/geometry.h"
#include "kagen/tools/rng_wrapper.h"

#include <algorithm>
#include <cstdlib>
#include <format>

namespace kagen {

double SampleHyperedgeRadius(
    SInt identifier, const PGeneratorConfig& config, double lower_bound, double upper_bound, RNGWrapper<>& rng) {
    const SInt   seed          = sampling::Spooky::hash(config.seed + (7919 * config.n) + identifier);
    const double uniformRandom = rng.GenerateUniformDouble(seed);

    return SampleHyperedgeRadiusFromUniform(config, uniformRandom, lower_bound, upper_bound);
}

double SampleHyperedgeRadius(
    const PGeneratorConfig& config, const double lower_bound, const double upper_bound, Mersenne& mersenne) {
    const double uniformRandom = mersenne.Random();

    return SampleHyperedgeRadiusFromUniform(config, uniformRandom, lower_bound, upper_bound);
}

double getSampledOrConstantRadius(
    const PGeneratorConfig& config, const SInt identifier, const double actual_lower_bound,
    const double actual_upper_bound, RNGWrapper<>& rng) {
    if (!config.random_radius) {
        return config.r;
    }

    return SampleHyperedgeRadius(identifier, config, actual_lower_bound, actual_upper_bound, rng);
}

double SampleHyperedgeRadiusFromUniform(
    const PGeneratorConfig& config, const double uniform_random, const double lower_bound, const double upper_bound) {
    if (lower_bound <= 0.0 || upper_bound <= 0.0 || lower_bound > upper_bound) {
        throw ConfigurationError(
            std::format("Invalid hyperedge radius bounds: lower bound:{} upper bound:{}", lower_bound, upper_bound));
    }

    const double alpha = config.hyperedge_radius_exponent;

    const double log_lower = -alpha * std::log(lower_bound);
    const double log_upper = -alpha * std::log(upper_bound);

    const double max_log = std::max(log_lower, log_upper);

    const double a = std::exp(log_lower - max_log);
    const double b = std::exp(log_upper - max_log);

    const double mixed = ((1.0 - uniform_random) * a) + (uniform_random * b);

    double sampled = std::exp(-(std::log(mixed) + max_log) / alpha);
    sampled        = std::clamp(sampled, lower_bound, upper_bound);

    if (!std::isfinite(sampled) || sampled <= 0.0) {
        throw ConfigurationError("Invalid sampled hyperedge radius");
    }

    return sampled;
}

/**
 * Performs checks to make sure that config is valid in terms of radius data.
 */
bool RandomRadiusChecks(PGeneratorConfig& config) {
    if (config.random_radius) {
        if (config.min_hyperedge_radius <= 0.0 || config.max_hyperedge_radius <= 0.0
            || config.hyperedge_radius_exponent <= 0.0) {
            throw ConfigurationError("non-constant hyperball radii require min radius, max radius and exponent > 0");
        }
        if (config.min_hyperedge_radius > config.max_hyperedge_radius) {
            throw ConfigurationError("min hyperedge radius must not exceed max hyperedge radius");
        }
    } else if (config.r <= 0) {
        throw ConfigurationError("hyperbolic hypergraph generator requires hyperedge radius > 0");
    }
    return true;
}

PinRange getRandomPinRange(
    const SInt target_cell_size, const SInt range_size, const SInt target_cell_offset, Mersenne& mersenne) {
    if (range_size < 0 || range_size > target_cell_size) {
        throw ConfigurationError("Cannot sample a range larger than the target cell");
    }

    if (range_size == target_cell_size) {
        return {
            .begin = target_cell_offset,
            .end   = target_cell_offset + target_cell_size,
        };
    }

    const LPFloat u                = mersenne.Random();
    const SInt    number_of_starts = target_cell_size - range_size + 1;

    const SInt interval_start = static_cast<SInt>(u * static_cast<LPFloat>(number_of_starts));

    return {
        .begin = target_cell_offset + interval_start,
        .end   = target_cell_offset + interval_start + range_size,
    };
}

double QuantileOrConstantHyperedgeRadius(const PGeneratorConfig& config) {
    if (!config.random_radius) {
        if (config.r <= 0.0 || !std::isfinite(config.r)) {
            throw ConfigurationError("Constant hyperedge radius must be positive and finite");
        }

        return config.r;
    }

    const double alpha = config.hyperedge_radius_exponent;
    if (alpha <= 0.0 || !std::isfinite(alpha)) {
        throw ConfigurationError("Random hyperedge radius requires a positive finite exponent");
    }

    const double default_lower = std::sqrt(2.0 / (M_PI * static_cast<double>(config.n)));

    const double lower = config.min_hyperedge_radius == -1.0 ? default_lower : config.min_hyperedge_radius;

    const double upper = config.max_hyperedge_radius == -1.0 ? 1.0 : config.max_hyperedge_radius;

    if (lower <= 0.0 || upper <= 0.0 || !std::isfinite(lower) || !std::isfinite(upper)) {
        throw ConfigurationError("Hyperedge radius bounds must be positive and finite");
    }

    if (lower > upper) {
        throw ConfigurationError("Minimum hyperedge radius must not exceed maximum hyperedge radius");
    }

    const double q = config.quantile;

    if (!std::isfinite(q) || q < 0.0 || q > 1.0) {
        throw ConfigurationError("Hyperedge radius quantile must be finite and in [0, 1]");
    }

    return SampleHyperedgeRadiusFromUniform(config, q, lower, upper);
}

double ExpectedSquaredHyperedgeRadius(const double lower, const double upper, const double alpha) {
    if (lower <= 0.0 || upper <= 0.0 || lower > upper || alpha <= 0.0 || !std::isfinite(alpha)) {
        throw ConfigurationError("invalid radius distribution parameters");
    }

    if (lower == upper) {
        return lower * lower;
    }

    constexpr double eps = 1e-12;

    double numerator;
    if (std::abs(alpha - 2.0) < eps) {
        numerator = std::log(upper / lower);
    } else {
        numerator = (std::pow(upper, 2.0 - alpha) - std::pow(lower, 2.0 - alpha)) / (2.0 - alpha);
    }

    const double denominator = (std::pow(lower, -alpha) - std::pow(upper, -alpha)) / alpha;

    return numerator / denominator;
}

double SolveRadiusExponentForExpectedPins(const double target_e_r2, const double lower, const double upper) {
    if (lower <= 0.0 || upper <= 0.0 || lower > upper) {
        throw ConfigurationError("invalid radius bounds");
    }

    const double min_e_r2 = lower * lower;
    const double max_e_r2 = upper * upper;

    if (target_e_r2 < min_e_r2 || target_e_r2 > max_e_r2) {
        throw ConfigurationError(
            std::format(
                "expected total pins target is incompatible with radius bounds: target E[r^2]={}, allowed=[{}, {}]",
                target_e_r2, min_e_r2, max_e_r2));
    }

    if (std::abs(target_e_r2 - min_e_r2) < 1e-15) {
        return 1024.0;
    }

    double lo = 1e-12;
    double hi = 1.0;

    while (ExpectedSquaredHyperedgeRadius(lower, upper, hi) > target_e_r2) {
        hi *= 2.0;

        if (hi > 1024.0) {
            throw ConfigurationError("could not solve hyperedge radius exponent");
        }
    }

    for (int iter = 0; iter < 80; ++iter) {
        const double mid = 0.5 * (lo + hi);
        const double cur = ExpectedSquaredHyperedgeRadius(lower, upper, mid);

        if (cur > target_e_r2) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    return 0.5 * (lo + hi);
}

double HyperbolicRadialQuantile(double q, double alpha, double min_r, double max_r) {
    const double min_cdf = std::cosh(alpha * min_r);
    const double max_cdf = std::cosh(alpha * max_r);

    return std::acosh((q * (max_cdf - min_cdf)) + min_cdf) / alpha;
}

struct RadialSample {
    double r;
    double sinh_r;
    double cosh_r;
    double weight;
};

std::vector<RadialSample> BuildRadialSamples(int samples, double alpha, double target_r) {
    std::vector<RadialSample> out;
    out.reserve(samples);

    for (int i = 0; i < samples; ++i) {
        const double q = (i + 0.5) / samples;
        const double r = HyperbolicRadialQuantile(q, alpha, 0.0, target_r);

        out.push_back({
            .r      = r,
            .sinh_r = std::sinh(r),
            .cosh_r = std::cosh(r),
            .weight = 1.0 / samples,
        });
    }

    return out;
}

double ExpectedPinsForCenterRadius(
    double center_r, double radius, const PGeneratorConfig& config, double target_r,
    const std::vector<RadialSample>& samples) {
    if (radius >= center_r + target_r) {
        return static_cast<double>(config.n);
    }

    if (radius <= 0.0) {
        return 0.0;
    }

    double expected = 0.0;

    // Region 1: all sample radii r with r <= radius - center_r are fully inside.
    const double full_inside_until = radius - center_r;

    auto full_inside_end = std::upper_bound(
        samples.begin(), samples.end(), full_inside_until,
        [](double value, const RadialSample& sample) { return value < sample.r; });

    for (auto it = samples.begin(); it != full_inside_end; ++it) {
        expected += it->weight;
    }

    // Region 2: only samples with |center_r - r| < radius < center_r + r need acos.
    const double boundary_begin_r = std::max(0.0, center_r - radius);
    const double boundary_end_r   = std::min(target_r, center_r + radius);

    auto boundary_begin = std::upper_bound(
        samples.begin(), samples.end(), boundary_begin_r,
        [](double value, const RadialSample& sample) { return value < sample.r; });

    auto boundary_end =
        std::lower_bound(samples.begin(), samples.end(), boundary_end_r, [](const RadialSample& sample, double value) {
            return sample.r < value;
        });

    boundary_begin = std::max(boundary_begin, full_inside_end);

    const double center_cosh = std::cosh(center_r);
    const double center_sinh = std::sinh(center_r);
    const double radius_cosh = std::cosh(radius);

    for (auto it = boundary_begin; it != boundary_end; ++it) {
        const auto& s = *it;

        if (center_sinh <= 0.0 || s.sinh_r <= 0.0) {
            if (std::abs(center_r - s.r) <= radius) {
                expected += s.weight;
            }
            continue;
        }

        const double x = ((center_cosh * s.cosh_r) - radius_cosh) / (center_sinh * s.sinh_r);

        if (x <= -1.0) {
            expected += s.weight;
        } else if (x < 1.0) {
            expected += s.weight * (std::acos(x) / M_PI);
        }
    }

    return static_cast<double>(config.n) * expected;
}

static inline double
SampleRadiusQuantile(const double q, const double lower, const double upper, const double exponent) {
    if (lower == upper) {
        return lower;
    }

    const double log_lower = -exponent * std::log(lower);
    const double log_upper = -exponent * std::log(upper);
    const double max_log   = std::max(log_lower, log_upper);

    const double scaled_lower = std::exp(log_lower - max_log);
    const double scaled_upper = std::exp(log_upper - max_log);

    const double mixed = std::lerp(scaled_lower, scaled_upper, q);

    return std::clamp(std::exp(-(std::log(mixed) + max_log) / exponent), lower, upper);
}

double ExpectedPinsForCenterRadiusGeneratorApprox(
    double center_r, double radius, const PGeneratorConfig& config, double alpha, double target_r) {
    const int total_annuli = std::max(1, static_cast<int>(std::floor(alpha * target_r / std::numbers::ln2)));

    const double total_area = PGGeometry<double>::RadiusToHyperbolicArea(alpha * target_r);

    double expected = 0.0;

    for (int a = 0; a < total_annuli; ++a) {
        const double min_r = a * target_r / total_annuli;
        const double max_r = (a + 1) * target_r / total_annuli;

        const double ring_area = PGGeometry<double>::RadiusToHyperbolicArea(alpha * max_r)
                                 - PGGeometry<double>::RadiusToHyperbolicArea(alpha * min_r);

        const double expected_annulus_n = static_cast<double>(config.n) * ring_area / total_area;

        const double query_r = std::clamp(center_r, min_r, max_r);
        const double denom   = std::sinh(center_r) * std::sinh(query_r);

        double angular_reach;

        if (denom <= 0.0) {
            angular_reach = M_PI;
        } else {
            const double x = ((std::cosh(center_r) * std::cosh(query_r)) - std::cosh(radius)) / denom;

            if (x <= -1.0) {
                angular_reach = M_PI;
            } else if (x >= 1.0) {
                angular_reach = 0.0;
            } else {
                angular_reach = std::acos(x);
            }
        }

        const double angular_fraction = std::min(1.0, angular_reach / M_PI);

        expected += expected_annulus_n * angular_fraction;
    }

    return expected;
}

namespace {

struct PreparedRadiusDistribution {
    double lower;
    double upper;
    double exponent;
    double max_log;
    double scaled_lower;
    double scaled_upper;

    PreparedRadiusDistribution(const double lower_, const double upper_, const double exponent_)
        : lower(lower_),
          upper(upper_),
          exponent(exponent_) {
        const double log_lower = -exponent * std::log(lower);
        const double log_upper = -exponent * std::log(upper);

        max_log      = std::max(log_lower, log_upper);
        scaled_lower = std::exp(log_lower - max_log);
        scaled_upper = std::exp(log_upper - max_log);
    }

    double Quantile(const double q) const {
        if (lower == upper) {
            return lower;
        }

        const double mixed = std::lerp(scaled_lower, scaled_upper, q);

        return std::clamp(std::exp(-(std::log(mixed) + max_log) / exponent), lower, upper);
    }
};

double ClampOpen01(double q) {
    constexpr double eps = 1e-14;
    return std::clamp(q, eps, 1.0 - eps);
}

template <typename F>
double Simpson(double a, double b, double fa, double fm, double fb) {
    return (b - a) * (fa + 4.0 * fm + fb) / 6.0;
}

template <typename F>
double AdaptiveSimpsonRecursive(
    const F& f, double a, double b, double eps, double fa, double fm, double fb, double whole, int depth) {
    const double m  = 0.5 * (a + b);
    const double lm = 0.5 * (a + m);
    const double rm = 0.5 * (m + b);

    const double flm = f(lm);
    const double frm = f(rm);

    const double left  = Simpson<F>(a, m, fa, flm, fm);
    const double right = Simpson<F>(m, b, fm, frm, fb);

    const double refined = left + right;
    const double error   = refined - whole;

    if (depth <= 0 || std::abs(error) <= 15.0 * eps) {
        return refined + error / 15.0;
    }

    return AdaptiveSimpsonRecursive(f, a, m, 0.5 * eps, fa, flm, fm, left, depth - 1)
           + AdaptiveSimpsonRecursive(f, m, b, 0.5 * eps, fm, frm, fb, right, depth - 1);
}

template <typename F>
double AdaptiveSimpson(const F& f, double a, double b, double eps, int max_depth = 24) {
    const double m  = 0.5 * (a + b);
    const double fa = f(a);
    const double fm = f(m);
    const double fb = f(b);

    const double whole = Simpson<F>(a, b, fa, fm, fb);

    return AdaptiveSimpsonRecursive(f, a, b, eps, fa, fm, fb, whole, max_depth);
}

template <typename F>
double IntegrateRadiusQuantileTailAware(const F& f) {
    static constexpr std::array<double, 11> cuts{
        0.0, 0.50, 0.90, 0.99, 0.999, 0.9999, 0.99999, 0.999999, 0.9999999, 0.99999999, 1.0,
    };

    double result = 0.0;

    for (std::size_t i = 0; i + 1 < cuts.size(); ++i) {
        const double a = cuts[i];
        const double b = cuts[i + 1];

        if (a < b) {
            result += AdaptiveSimpson(f, a, b, 1e-4);
        }
    }

    return result;
}

double ExpectedPinsOverRadius(
    double center_r, double lower, double upper, double exponent, const PGeneratorConfig& config, double target_r,
    const std::vector<RadialSample>& radial_samples) {
    if (lower <= 0.0 || upper <= 0.0 || lower > upper) {
        throw ConfigurationError("invalid hyperedge radius bounds in expected pin solver");
    }

    if (lower == upper) {
        return ExpectedPinsForCenterRadius(center_r, lower, config, target_r, radial_samples);
    }

    const PreparedRadiusDistribution distribution{lower, upper, exponent};
    auto                             f = [&](double q) {
        q = ClampOpen01(q);

        const double radius = distribution.Quantile(q);

        return ExpectedPinsForCenterRadius(center_r, radius, config, target_r, radial_samples);
    };

    return IntegrateRadiusQuantileTailAware(f);
}

} // namespace

double SolveHyperbolicRadiusExponentForExpectedPins(const PGeneratorConfig& config) {
    const double alpha = (config.plexp - 1.0) / 2.0;

    const double target_r = PGGeometry<double>::GetTargetRadius(config.n, config.n * config.avg_degree / 2.0, alpha);

    const double target_per_edge = static_cast<double>(config.size_dist_pin_budget) / static_cast<double>(config.m);

    if (config.size_dist_upper_bound > 0 && config.size_dist_lower_bound > config.size_dist_upper_bound) {
        throw ConfigurationError(
            "lower hyperedge size bound "
            "must not exceed upper hyperedge size bound");
    }

    const int total_annuli = std::max(1, static_cast<int>(std::floor(alpha * target_r / std::numbers::ln2)));

    const auto radial_samples = BuildRadialSamples(64, alpha, target_r);

    //
    // Find the hyperbolic radius whose hyperball has the requested
    // expected number of pins for a center at center_r.
    //
    auto radius_for_expected_pins = [&](const double center_r, const double desired_pins) {
        if (desired_pins <= 0.0) {
            throw ConfigurationError(
                "expected hyperedge size "
                "must be positive");
        }

        if (desired_pins >= static_cast<double>(config.n)) {
            return center_r + target_r;
        }

        double lo = 0.0;
        double hi = center_r + target_r;

        for (int i = 0; i < 48; ++i) {
            const double mid = 0.5 * (lo + hi);

            const double expected = ExpectedPinsForCenterRadiusGeneratorApprox(center_r, mid, config, alpha, target_r);

            if (expected >= desired_pins) {
                hi = mid;
            } else {
                lo = mid;
            }
        }

        return hi;
    };

    //
    // Generation approximates center-dependent size bounds once per
    // center annulus, using the annulus midpoint. Reproduce that here
    // so the pin-budget solver uses the same radius distribution.
    //
    std::vector<double> lower_by_annulus(total_annuli);

    std::vector<double> upper_by_annulus(total_annuli);

    for (int a = 0; a < total_annuli; ++a) {
        const double min_r = a * target_r / total_annuli;

        const double max_r = (a + 1) * target_r / total_annuli;

        const double mid_r = 0.5 * (min_r + max_r);

        if (config.min_hyperedge_radius != -1.0) {
            lower_by_annulus[a] = config.min_hyperedge_radius;
        } else {
            lower_by_annulus[a] = radius_for_expected_pins(mid_r, static_cast<double>(config.size_dist_lower_bound));
        }

        if (config.max_hyperedge_radius != -1.0) {
            upper_by_annulus[a] = config.max_hyperedge_radius;
        } else if (config.size_dist_upper_bound > 0) {
            upper_by_annulus[a] = radius_for_expected_pins(mid_r, static_cast<double>(config.size_dist_upper_bound));
        } else {
            //
            // Sentinel: no expected-size upper bound.
            // The actual upper radius will depend on the sampled
            // center position.
            //
            upper_by_annulus[a] = -1.0;
        }
    }

    auto expected_per_edge = [&](double exponent) {
        const double total_area = PGGeometry<double>::RadiusToHyperbolicArea(alpha * target_r);

        double weighted_total = 0.0;
        double total_weight   = 0.0;

        constexpr int CENTER_SAMPLES_PER_ANNULUS = 4;

        for (int a = 0; a < total_annuli; ++a) {
            const double ann_min_r = a * target_r / total_annuli;

            const double ann_max_r = (a + 1) * target_r / total_annuli;

            const double ring_area = PGGeometry<double>::RadiusToHyperbolicArea(alpha * ann_max_r)
                                     - PGGeometry<double>::RadiusToHyperbolicArea(alpha * ann_min_r);

            const double annulus_weight = ring_area / total_area;

            const double lower = lower_by_annulus[a];

            for (int i = 0; i < CENTER_SAMPLES_PER_ANNULUS; ++i) {
                const double q0 = static_cast<double>(i) / CENTER_SAMPLES_PER_ANNULUS;

                const double q1 = static_cast<double>(i + 1) / CENTER_SAMPLES_PER_ANNULUS;

                const double qc = 0.5 * (q0 + q1);

                const double center_r = HyperbolicRadialQuantile(qc, alpha, ann_min_r, ann_max_r);

                double upper;

                if (upper_by_annulus[a] > 0.0) {
                    upper = upper_by_annulus[a];
                } else {
                    upper = center_r + target_r;
                }

                upper = std::max(lower, upper);

                const double expected =
                    ExpectedPinsOverRadius(center_r, lower, upper, exponent, config, target_r, radial_samples);

                const double weight = annulus_weight / CENTER_SAMPLES_PER_ANNULUS;

                weighted_total += weight * expected;

                total_weight += weight;
            }
        }

        return weighted_total / total_weight;
    };

    double lo = 1e-6;
    double hi = 1.0;

    while (expected_per_edge(hi) > target_per_edge) {
        hi *= 2.0;

        if (hi > 1024.0) {
            throw ConfigurationError("could not bracket hyperedge radius exponent");
        }
    }

    for (int iter = 0; iter < 32; ++iter) {
        const double mid = 0.5 * (lo + hi);

        const double cur = expected_per_edge(mid);

        const double relative_error = std::abs(cur - target_per_edge) / std::max(1.0, target_per_edge);

        if (relative_error <= 1e-4) {
            lo = mid;
            hi = mid;
            break;
        }

        if ((hi - lo) <= 1e-8 * std::max(1.0, mid)) {
            break;
        }

        if (cur > target_per_edge) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    const double exponent = 0.5 * (lo + hi);

    return exponent;
}

double EuclideanRadiusForExpectedHyperedgeSize(const SInt expected_size, const SInt num_vertices) {
    if (expected_size <= 0) {
        throw ConfigurationError("expected hyperedge size must be positive");
    }

    if (num_vertices <= 0) {
        throw ConfigurationError("number of vertices must be positive");
    }

    return std::sqrt(static_cast<double>(expected_size) / (M_PI * static_cast<double>(num_vertices)));
}

double ExpectedHyperedgeSizeForEuclideanRadius(const double radius, const SInt num_vertices) {
    return (static_cast<double>(num_vertices) * M_PI * radius * radius);
}

} // namespace kagen