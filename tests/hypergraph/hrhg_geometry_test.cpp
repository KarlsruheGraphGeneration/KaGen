#include "kagen/generators/hyper/h_hyperbolic/hyper_hyperbolic.h"
#include "kagen/generators/hyper/h_hyperbolic/hyper_hyperbolic_policy.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <vector>

using namespace kagen;

namespace {

using Generator = LowPrecisionHyperHyperbolic;
using Policy    = HyperbolicGeometryPolicy<LPFloat>;
using Center    = HyperbolicHyperedgeCenter<LPFloat>;
using Cell      = Policy::Cell;

PGeneratorConfig GeometryConfig(const double radius) {
    PGeneratorConfig config;

    config.generator = GeneratorType::H_RHG;

    config.n = 5000;
    config.m = 1000;
    config.k = 4;

    config.avg_degree = 16.0;
    config.plexp      = 2.8;

    config.r = radius;

    config.seed  = 42;
    config.quiet = true;
    config.debug = false;

    config.random_radius = false;

    config.partial_cell_mode = PartialCellMode::GenerateAndCheck;

    Hyper_HyperbolicFactory factory;

    return factory.NormalizeParameters(config, 0, 1, false);
}

struct GeometryFixture {
    explicit GeometryFixture(const double radius)
        : config(GeometryConfig(radius)),
          generator(config, 0, 1),
          policy(generator, 0, 0) {}

    PGeneratorConfig config;
    Generator        generator;
    Policy           policy;
};

double NormalizePhi(double phi) {
    const double two_pi = 2.0 * M_PI;

    while (phi < 0.0) {
        phi += two_pi;
    }

    while (phi >= two_pi) {
        phi -= two_pi;
    }

    return phi;
}

double AngularDistance(double a, double b) {
    a = NormalizePhi(a);
    b = NormalizePhi(b);

    const double distance = std::abs(a - b);

    return std::min(distance, (2.0 * M_PI) - distance);
}

bool AngleInCell(const double phi, const Cell& cell) {
    const double normalized_phi = NormalizePhi(phi);

    const double min_phi = NormalizePhi(cell.min_phi);

    const double max_phi = NormalizePhi(cell.max_phi);

    if (min_phi <= max_phi) {
        return min_phi <= normalized_phi && normalized_phi <= max_phi;
    }

    return normalized_phi >= min_phi || normalized_phi <= max_phi;
}

double MinAngularDistance(const double phi, const Cell& cell) {
    if (AngleInCell(phi, cell)) {
        return 0.0;
    }

    return std::min(AngularDistance(phi, cell.min_phi), AngularDistance(phi, cell.max_phi));
}

double MaxAngularDistance(const double phi, const Cell& cell) {
    const double antipode = NormalizePhi(phi + M_PI);

    if (AngleInCell(antipode, cell)) {
        return M_PI;
    }

    return std::max(AngularDistance(phi, cell.min_phi), AngularDistance(phi, cell.max_phi));
}

double CoshHyperbolicDistance(const double center_r, const double query_r, const double delta_phi) {
    return (std::cosh(center_r) * std::cosh(query_r))
           - (std::sinh(center_r) * std::sinh(query_r) * std::cos(delta_phi));
}

//
// Exact minimum cosh-distance from the hyperedge center to
// any point of the cell.
//
// For fixed angular distance delta:
//
//   cosh(d)
//     = cosh(rc) cosh(r)
//       - sinh(rc) sinh(r) cos(delta)
//
// The radial minimum may occur at an endpoint or at
//
//   tanh(r) = tanh(rc) cos(delta).
//
double MinimumCoshDistance(const Center& center, const Cell& cell) {
    const double delta = MinAngularDistance(center.phi, cell);

    double best = std::min(
        CoshHyperbolicDistance(center.r, cell.min_r, delta), CoshHyperbolicDistance(center.r, cell.max_r, delta));

    const double stationary_tanh = std::tanh(center.r) * std::cos(delta);

    if (stationary_tanh > 0.0 && stationary_tanh < 1.0) {
        const double stationary_r = std::atanh(stationary_tanh);

        if (stationary_r >= cell.min_r && stationary_r <= cell.max_r) {
            best = std::min(best, CoshHyperbolicDistance(center.r, stationary_r, delta));
        }
    }

    return best;
}

//
// For fixed angular distance, cosh-distance is convex in r.
// Therefore its maximum over the radial interval occurs at
// one of the two radial endpoints.
//
// The maximum angular distance independently maximizes the
// distance because cos(delta) decreases as delta approaches pi.
//
double MaximumCoshDistance(const Center& center, const Cell& cell) {
    const double delta = MaxAngularDistance(center.phi, cell);

    return std::max(
        CoshHyperbolicDistance(center.r, cell.min_r, delta), CoshHyperbolicDistance(center.r, cell.max_r, delta));
}

CellBallRelation ExactRelation(const Center& center, const double radius, const Cell& cell) {
    const double cosh_radius = std::cosh(radius);

    if (MinimumCoshDistance(center, cell) > cosh_radius) {
        return CellBallRelation::OUTSIDE;
    }

    if (MaximumCoshDistance(center, cell) <= cosh_radius) {
        return CellBallRelation::INSIDE;
    }

    return CellBallRelation::PARTIAL;
}

//
// Obtain the complete cell universe through the same public
// geometry interface.
//
// Centering at r = 0 with a radius much larger than the
// generated disk gives full radial and angular reach.
//
std::vector<Cell> AllCells(Policy& policy) {
    const Center center{
        .phi        = 0.0,
        .r          = 0.0,
        .sampled_id = 0,
        .annulus_id = 0,
    };

    std::vector<Cell> cells;

    policy.CandidateCells(center, 100.0, cells);

    return cells;
}

std::set<SInt> CellIds(const std::vector<Cell>& cells) {
    std::set<SInt> ids;

    for (const Cell& cell: cells) {
        ids.insert(cell.global_cell_id);
    }

    return ids;
}

double IntervalOverlap(const double a_begin, const double a_end, const double b_begin, const double b_end) {
    return std::max(0.0, std::min(a_end, b_end) - std::max(a_begin, b_begin));
}

double ExactCircularOverlap(double a_begin, double a_end, double b_begin, double b_end) {
    const double two_pi = 2.0 * M_PI;

    auto split = [two_pi](double begin, double end) {
        std::vector<std::pair<double, double>> parts;

        const double width = end - begin;

        if (width >= two_pi) {
            parts.emplace_back(0.0, two_pi);

            return parts;
        }

        begin = NormalizePhi(begin);
        end   = NormalizePhi(end);

        if (begin <= end) {
            parts.emplace_back(begin, end);
        } else {
            parts.emplace_back(begin, two_pi);

            parts.emplace_back(0.0, end);
        }

        return parts;
    };

    const auto a_parts = split(a_begin, a_end);

    const auto b_parts = split(b_begin, b_end);

    double overlap = 0.0;

    for (const auto& a: a_parts) {
        for (const auto& b: b_parts) {
            overlap += IntervalOverlap(a.first, a.second, b.first, b.second);
        }
    }

    return overlap;
}

double ExactAngularCoverageAtRadius(
    const Center& center, const double hyperball_radius, const Cell& cell, const double query_r) {
    const double cell_width = cell.max_phi - cell.min_phi;

    if (cell_width <= 0.0) {
        return 0.0;
    }

    const double denominator = std::sinh(center.r) * std::sinh(query_r);

    double half_angle;

    if (denominator <= std::numeric_limits<double>::epsilon()) {
        const double radial_distance = std::abs(center.r - query_r);

        half_angle = radial_distance <= hyperball_radius ? M_PI : 0.0;
    } else {
        const double argument =
            ((std::cosh(center.r) * std::cosh(query_r)) - std::cosh(hyperball_radius)) / denominator;

        if (argument <= -1.0) {
            half_angle = M_PI;
        } else if (argument >= 1.0) {
            half_angle = 0.0;
        } else {
            half_angle = std::acos(argument);
        }
    }

    if (half_angle >= M_PI) {
        return 1.0;
    }

    if (half_angle <= 0.0) {
        return 0.0;
    }

    const double overlap =
        ExactCircularOverlap(cell.min_phi, cell.max_phi, center.phi - half_angle, center.phi + half_angle);

    return std::clamp(overlap / cell_width, 0.0, 1.0);
}

double Simpson(const double a, const double b, const double fa, const double fm, const double fb) {
    return (b - a) * (fa + (4.0 * fm) + fb) / 6.0;
}

template <typename Function>
double AdaptiveSimpsonRecursive(
    const Function& function, const double a, const double b, const double fa, const double fm, const double fb,
    const double whole, const double tolerance, const int depth) {
    const double midpoint = 0.5 * (a + b);

    const double left_midpoint = 0.5 * (a + midpoint);

    const double right_midpoint = 0.5 * (midpoint + b);

    const double f_left_midpoint = function(left_midpoint);

    const double f_right_midpoint = function(right_midpoint);

    const double left = Simpson(a, midpoint, fa, f_left_midpoint, fm);

    const double right = Simpson(midpoint, b, fm, f_right_midpoint, fb);

    const double combined = left + right;

    if (depth <= 0 || std::abs(combined - whole) <= 15.0 * tolerance) {
        return combined + ((combined - whole) / 15.0);
    }

    return AdaptiveSimpsonRecursive(function, a, midpoint, fa, f_left_midpoint, fm, left, tolerance / 2.0, depth - 1)
           + AdaptiveSimpsonRecursive(
               function, midpoint, b, fm, f_right_midpoint, fb, right, tolerance / 2.0, depth - 1);
}

template <typename Function>
double AdaptiveSimpson(const Function& function, const double a, const double b, const double tolerance = 1e-10) {
    const double midpoint = 0.5 * (a + b);

    const double fa = function(a);

    const double fm = function(midpoint);

    const double fb = function(b);

    const double whole = Simpson(a, b, fa, fm, fb);

    return AdaptiveSimpsonRecursive(function, a, b, fa, fm, fb, whole, tolerance, 20);
}

double ReferenceCellCoverage(
    const PGeneratorConfig& config, const Center& center, const double hyperball_radius, const Cell& cell) {
    const double alpha = (config.plexp - 1.0) / 2.0;

    const double u_min = std::cosh(alpha * cell.min_r);

    const double u_max = std::cosh(alpha * cell.max_r);

    if (u_max <= u_min) {
        return 0.0;
    }

    auto integrand = [&](const double u) {
        const double query_r = std::acosh(u) / alpha;

        return ExactAngularCoverageAtRadius(center, hyperball_radius, cell, query_r);
    };

    const double integral = AdaptiveSimpson(integrand, u_min, u_max, 1e-10);

    return std::clamp(integral / (u_max - u_min), 0.0, 1.0);
}

} // namespace

//
// Candidate collection
//

TEST(HRHGGeometry, CandidateCollectionMissesNoIntersectingCells) {
    constexpr double radius = 1.0;

    GeometryFixture fixture(radius);

    const auto all_cells = AllCells(fixture.policy);

    ASSERT_FALSE(all_cells.empty());

    const Center centers[] = {
        {
            .phi        = 1.0,
            .r          = 2.0,
            .sampled_id = 0,
            .annulus_id = 0,
        },

        //
        // Close to 0: angular search wraps around.
        //
        {
            .phi        = 0.01,
            .r          = 2.0,
            .sampled_id = 1,
            .annulus_id = 0,
        },

        //
        // Same edge case from the other side of 2pi.
        //
        {
            .phi        = (2.0 * M_PI) - 0.01,
            .r          = 3.0,
            .sampled_id = 2,
            .annulus_id = 0,
        },

        {
            .phi        = M_PI,
            .r          = 1.0,
            .sampled_id = 3,
            .annulus_id = 0,
        },
    };

    for (const Center& center: centers) {
        std::vector<Cell> candidates;

        fixture.policy.CandidateCells(center, radius, candidates);

        const auto candidate_ids = CellIds(candidates);

        for (const Cell& cell: all_cells) {
            const CellBallRelation exact = ExactRelation(center, radius, cell);

            if (exact == CellBallRelation::OUTSIDE) {
                continue;
            }

            EXPECT_TRUE(candidate_ids.contains(cell.global_cell_id))
                << "candidate collection missed intersecting cell"
                << " annulus=" << cell.annulus_id << " chunk=" << cell.chunk_id << " cell=" << cell.cell_id
                << " min_r=" << cell.min_r << " max_r=" << cell.max_r << " min_phi=" << cell.min_phi
                << " max_phi=" << cell.max_phi << " center_phi=" << center.phi << " center_r=" << center.r
                << " radius=" << radius;
        }
    }
}

TEST(HRHGGeometry, CandidateCellsRespectRadialSearchRange) {
    constexpr double radius = 0.5;

    GeometryFixture fixture(radius);

    const Center center{
        .phi        = 1.0,
        .r          = 2.0,
        .sampled_id = 0,
        .annulus_id = 0,
    };

    std::vector<Cell> cells;

    fixture.policy.CandidateCells(center, radius, cells);

    ASSERT_FALSE(cells.empty());

    for (const Cell& cell: cells) {
        EXPECT_GE(cell.max_r, center.r - radius);

        EXPECT_LE(cell.min_r, center.r + radius);
    }
}

TEST(HRHGGeometry, CandidateCellsHandleAngularWraparound) {
    constexpr double radius = 1.0;

    GeometryFixture fixture(radius);

    const Center center{
        .phi        = 0.01,
        .r          = 2.0,
        .sampled_id = 0,
        .annulus_id = 0,
    };

    std::vector<Cell> cells;

    fixture.policy.CandidateCells(center, radius, cells);

    ASSERT_FALSE(cells.empty());

    bool has_low_phi  = false;
    bool has_high_phi = false;

    for (const Cell& cell: cells) {
        if (cell.min_phi < 0.5) {
            has_low_phi = true;
        }

        if (cell.max_phi > (2.0 * M_PI) - 0.5) {
            has_high_phi = true;
        }
    }

    EXPECT_TRUE(has_low_phi);

    EXPECT_TRUE(has_high_phi);
}

//
// Classification
//

TEST(HRHGGeometry, ClassificationMatchesExactHyperbolicRelation) {
    constexpr double radius = 1.0;

    GeometryFixture fixture(radius);

    const Center centers[] = {
        {
            .phi        = 1.0,
            .r          = 2.0,
            .sampled_id = 0,
            .annulus_id = 0,
        },

        {
            .phi        = 0.01,
            .r          = 2.0,
            .sampled_id = 1,
            .annulus_id = 0,
        },

        {
            .phi        = (2.0 * M_PI) - 0.01,
            .r          = 3.0,
            .sampled_id = 2,
            .annulus_id = 0,
        },

        {
            .phi        = 2.3,
            .r          = 4.0,
            .sampled_id = 3,
            .annulus_id = 0,
        },
    };

    for (const Center& center: centers) {
        std::vector<Cell> candidates;

        fixture.policy.CandidateCells(center, radius, candidates);

        ASSERT_FALSE(candidates.empty());

        for (const Cell& cell: candidates) {
            const CellBallRelation expected = ExactRelation(center, radius, cell);

            const CellBallRelation actual = fixture.policy.ClassifyCell(center, radius, cell);

            EXPECT_EQ(actual, expected) << "incorrect classification"
                                        << " annulus=" << cell.annulus_id << " chunk=" << cell.chunk_id
                                        << " cell=" << cell.cell_id << " min_r=" << cell.min_r
                                        << " max_r=" << cell.max_r << " min_phi=" << cell.min_phi
                                        << " max_phi=" << cell.max_phi << " center_phi=" << center.phi
                                        << " center_r=" << center.r << " radius=" << radius;
        }
    }
}

//
// Coverage smoke / edge-case tests.
//
// The detailed precision test comes next; that requires an
// independent radial integration oracle rather than simply
// reimplementing CellCoverage().
//

TEST(HRHGGeometry, CoverageAlwaysLiesInUnitInterval) {
    constexpr double radius = 1.0;

    GeometryFixture fixture(radius);

    const Center center{
        .phi        = 1.0,
        .r          = 2.0,
        .sampled_id = 0,
        .annulus_id = 0,
    };

    std::vector<Cell> cells;

    fixture.policy.CandidateCells(center, radius, cells);

    ASSERT_FALSE(cells.empty());

    for (const Cell& cell: cells) {
        const double coverage = fixture.policy.CellCoverage(center, radius, cell);

        EXPECT_GE(coverage, 0.0);

        EXPECT_LE(coverage, 1.0);
    }
}

TEST(HRHGGeometry, InsideCellsHaveFullCoverage) {
    constexpr double radius = 2.0;

    GeometryFixture fixture(radius);

    const Center center{
        .phi        = 1.0,
        .r          = 2.0,
        .sampled_id = 0,
        .annulus_id = 0,
    };

    std::vector<Cell> cells;

    fixture.policy.CandidateCells(center, radius, cells);

    ASSERT_FALSE(cells.empty());

    bool found_inside = false;

    for (const Cell& cell: cells) {
        if (ExactRelation(center, radius, cell) != CellBallRelation::INSIDE) {
            continue;
        }

        found_inside = true;

        EXPECT_NEAR(fixture.policy.CellCoverage(center, radius, cell), 1.0, 1e-12)
            << "inside cell did not receive full coverage"
            << " annulus=" << cell.annulus_id << " chunk=" << cell.chunk_id << " cell=" << cell.cell_id;
    }

    EXPECT_TRUE(found_inside);
}

TEST(HRHGGeometry, PartialCellCoverageMatchesHighPrecisionReference) {
    constexpr double radius = 1.0;

    GeometryFixture fixture(radius);

    const Center centers[] = {
        {
            .phi        = 1.0,
            .r          = 2.0,
            .sampled_id = 0,
            .annulus_id = 0,
        },

        {
            .phi        = 0.01,
            .r          = 2.0,
            .sampled_id = 1,
            .annulus_id = 0,
        },

        {
            .phi        = 2.3,
            .r          = 4.0,
            .sampled_id = 2,
            .annulus_id = 0,
        },
    };

    std::size_t tested_partial_cells = 0;

    for (const Center& center: centers) {
        std::vector<Cell> cells;

        fixture.policy.CandidateCells(center, radius, cells);

        for (const Cell& cell: cells) {
            if (ExactRelation(center, radius, cell) != CellBallRelation::PARTIAL) {
                continue;
            }

            const double reference = ReferenceCellCoverage(fixture.config, center, radius, cell);

            const double estimated = fixture.policy.CellCoverage(center, radius, cell);

            EXPECT_NEAR(estimated, reference, 0.01)
                << "coverage error"
                << " annulus=" << cell.annulus_id << " chunk=" << cell.chunk_id << " cell=" << cell.cell_id
                << " min_r=" << cell.min_r << " max_r=" << cell.max_r << " min_phi=" << cell.min_phi
                << " max_phi=" << cell.max_phi << " center_phi=" << center.phi << " center_r=" << center.r
                << " estimated=" << estimated << " reference=" << reference
                << " absolute_error=" << std::abs(estimated - reference);

            ++tested_partial_cells;
        }
    }

    EXPECT_GT(tested_partial_cells, std::size_t{0});
}