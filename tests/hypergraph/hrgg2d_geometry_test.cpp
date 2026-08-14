#include "kagen/context.h"
#include "kagen/generators/hyper/h_geometric/h_rgg.h"
#include "kagen/generators/hyper/h_geometric/hyper_rgg2d.h"
#include "kagen/generators/hyper/h_geometric/hyper_rgg2d_policy.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <set>
#include <utility>
#include <vector>

using namespace kagen;

namespace {

constexpr SInt kNumVertices = 1024;
constexpr SInt kNumChunks   = 16;

struct CellBounds {
    double min_x;
    double max_x;
    double min_y;
    double max_y;
};

PGeneratorConfig GeometryConfig() {
    PGeneratorConfig config;

    config.n = kNumVertices;
    config.m = 1;
    config.k = kNumChunks;

    config.seed  = 1;
    config.quiet = true;
    config.debug = false;

    config.random_radius = false;
    config.r             = 0.1;

    config.partial_cell_mode = PartialCellMode::GenerateAndCheck;

    HyperRGG2DFactory factory;

    return factory.NormalizeParameters(config, 0, 1, false);
}

struct GeometryHarness {
    GeometryHarness() : config(GeometryConfig()), generator(config, 0, 1), policy(generator) {}

    PGeneratorConfig config;
    HyperRGG2D       generator;
    HyperRGG2DPolicy policy;
};

SInt TotalCellsPerDimension(HyperRGG2DPolicy& policy) {
    const HyperRGG2DPolicy::Center center{
        .x          = 0.5,
        .y          = 0.5,
        .radius     = 2.0,
        .sampled_id = 0,
        .chunk_id   = 0,
        .cell_id    = 0,
    };

    std::vector<HyperRGG2DPolicy::Cell> cells;

    policy.CandidateCells(center, 2.0, cells);

    const auto cells_per_dim = static_cast<SInt>(std::llround(std::sqrt(static_cast<double>(cells.size()))));

    EXPECT_EQ(cells_per_dim * cells_per_dim, cells.size());

    return cells_per_dim;
}

CellBounds BoundsOf(const SInt x, const SInt y, const SInt cells_per_dim) {
    const double h = 1.0 / static_cast<double>(cells_per_dim);

    return {
        .min_x = static_cast<double>(x) * h,
        .max_x = static_cast<double>(x + 1) * h,
        .min_y = static_cast<double>(y) * h,
        .max_y = static_cast<double>(y + 1) * h,
    };
}

CellBounds BoundsOf(const HyperRGG2DPolicy::Cell& cell, const SInt cells_per_dim) {
    return BoundsOf(cell.global_cell_x, cell.global_cell_y, cells_per_dim);
}

CellBallRelation
ExactRelation(const double center_x, const double center_y, const double radius, const CellBounds& cell) {
    const double radius_sq = radius * radius;

    //
    // Minimum squared distance from center to rectangle.
    //
    double dx = 0.0;

    if (center_x < cell.min_x) {
        dx = cell.min_x - center_x;
    } else if (center_x > cell.max_x) {
        dx = center_x - cell.max_x;
    }

    double dy = 0.0;

    if (center_y < cell.min_y) {
        dy = cell.min_y - center_y;
    } else if (center_y > cell.max_y) {
        dy = center_y - cell.max_y;
    }

    const double min_distance_sq = dx * dx + dy * dy;

    if (min_distance_sq > radius_sq) {
        return CellBallRelation::OUTSIDE;
    }

    //
    // Maximum squared distance occurs at one of the four corners.
    //
    const double corner_x[] = {
        cell.min_x,
        cell.max_x,
    };

    const double corner_y[] = {
        cell.min_y,
        cell.max_y,
    };

    double max_distance_sq = 0.0;

    for (const double x: corner_x) {
        for (const double y: corner_y) {
            const double cx = x - center_x;

            const double cy = y - center_y;

            max_distance_sq = std::max(max_distance_sq, cx * cx + cy * cy);
        }
    }

    if (max_distance_sq <= radius_sq) {
        return CellBallRelation::INSIDE;
    }

    return CellBallRelation::PARTIAL;
}

std::set<std::pair<SInt, SInt>> CandidateCoordinates(const std::vector<HyperRGG2DPolicy::Cell>& cells) {
    std::set<std::pair<SInt, SInt>> result;

    for (const auto& cell: cells) {
        result.emplace(cell.global_cell_x, cell.global_cell_y);
    }

    return result;
}

double CirclePrimitive(const double x, const double radius) {
    const double clamped_x = std::clamp(x, -radius, radius);

    const double inside = std::max(0.0, (radius * radius) - (clamped_x * clamped_x));

    return 0.5 * ((clamped_x * std::sqrt(inside)) + (radius * radius * std::asin(clamped_x / radius)));
}

double IntegralUpperCircle(const double x0, const double x1, const double radius) {
    if (x1 <= x0) {
        return 0.0;
    }

    const double begin = std::max(x0, -radius);

    const double end = std::min(x1, radius);

    if (end <= begin) {
        return 0.0;
    }

    return CirclePrimitive(end, radius) - CirclePrimitive(begin, radius);
}

double CircleAreaBelowY(const double x0, const double x1, const double y, const double radius) {
    if (x1 <= x0 || radius <= 0.0) {
        return 0.0;
    }

    const double begin = std::max(x0, -radius);

    const double end = std::min(x1, radius);

    if (end <= begin) {
        return 0.0;
    }

    //
    // Entire vertical circle section lies below y.
    //
    if (y >= radius) {
        return 2.0 * IntegralUpperCircle(begin, end, radius);
    }

    //
    // Entire vertical circle section lies above y.
    //
    if (y <= -radius) {
        return 0.0;
    }

    //
    // Solve
    //
    //     sqrt(r² - x²) = |y|
    //
    // for x.
    //
    const double crossing = std::sqrt(std::max(0.0, (radius * radius) - (y * y)));

    double area = 0.0;

    //
    // Handle the interval in pieces because the horizontal
    // line y intersects the circle at x = +/- crossing.
    //
    const double cuts[] = {
        begin,
        -crossing,
        crossing,
        end,
    };

    for (int i = 0; i < 3; ++i) {
        const double a = std::max(begin, cuts[i]);

        const double b = std::min(end, cuts[i + 1]);

        if (b <= a) {
            continue;
        }

        const double mid = 0.5 * (a + b);

        const double circle_height = std::sqrt(std::max(0.0, (radius * radius) - (mid * mid)));

        if (y <= -circle_height) {
            //
            // Rectangle threshold is below the circle here:
            // no circle area lies below y.
            //
            continue;
        }

        if (y >= circle_height) {
            //
            // Entire vertical circle section lies below y.
            //
            area += 2.0 * IntegralUpperCircle(a, b, radius);

            continue;
        }

        //
        // y cuts through the vertical circle section.
        //
        // Area below y is:
        //
        //     y - (-sqrt(...))
        //   = y + sqrt(...)
        //
        area += (y * (b - a)) + IntegralUpperCircle(a, b, radius);
    }

    return area;
}

double ExactCircleRectIntersectionArea(
    const double center_x, const double center_y, const double radius, const CellBounds& cell) {
    if (radius <= 0.0) {
        return 0.0;
    }

    //
    // Translate rectangle so that the circle center is at
    // the origin.
    //
    const double x0 = cell.min_x - center_x;

    const double x1 = cell.max_x - center_x;

    const double y0 = cell.min_y - center_y;

    const double y1 = cell.max_y - center_y;

    const double below_top = CircleAreaBelowY(x0, x1, y1, radius);

    const double below_bottom = CircleAreaBelowY(x0, x1, y0, radius);

    return std::max(0.0, below_top - below_bottom);
}

double
ExactCircleRectCoverage(const double center_x, const double center_y, const double radius, const CellBounds& cell) {
    const double cell_area = (cell.max_x - cell.min_x) * (cell.max_y - cell.min_y);

    if (cell_area <= 0.0) {
        return 0.0;
    }

    const double intersection = ExactCircleRectIntersectionArea(center_x, center_y, radius, cell);

    return std::clamp(intersection / cell_area, 0.0, 1.0);
}

const HyperRGG2DPolicy::Cell* FindCell(const std::vector<HyperRGG2DPolicy::Cell>& cells, const SInt x, const SInt y) {
    const auto it = std::find_if(cells.begin(), cells.end(), [=](const auto& cell) {
        return cell.global_cell_x == x && cell.global_cell_y == y;
    });

    if (it == cells.end()) {
        return nullptr;
    }

    return &*it;
}

} // namespace

//
// Candidate completeness
//

TEST(HRGG2DGeometry, CandidateCollectionMissesNoIntersectingCells) {
    GeometryHarness harness;

    const SInt cells_per_dim = TotalCellsPerDimension(harness.policy);

    const double test_cases[][3] = {
        {0.50, 0.50, 0.10}, {0.31, 0.67, 0.08}, {0.02, 0.02, 0.12},
        {0.98, 0.50, 0.07}, {0.50, 0.98, 0.09}, {0.123, 0.456, 0.035},
    };

    for (const auto& test: test_cases) {
        const double center_x = test[0];
        const double center_y = test[1];
        const double radius   = test[2];

        const HyperRGG2DPolicy::Center center{
            .x          = center_x,
            .y          = center_y,
            .radius     = radius,
            .sampled_id = 0,
            .chunk_id   = 0,
            .cell_id    = 0,
        };

        std::vector<HyperRGG2DPolicy::Cell> candidates;

        harness.policy.CandidateCells(center, radius, candidates);

        const auto candidate_coordinates = CandidateCoordinates(candidates);

        for (SInt x = 0; x < cells_per_dim; ++x) {
            for (SInt y = 0; y < cells_per_dim; ++y) {
                const CellBounds bounds = BoundsOf(x, y, cells_per_dim);

                const CellBallRelation exact = ExactRelation(center_x, center_y, radius, bounds);

                if (exact == CellBallRelation::OUTSIDE) {
                    continue;
                }

                EXPECT_TRUE(candidate_coordinates.contains({x, y}))
                    << "Candidate collection missed intersecting cell "
                    << "(" << x << ", " << y << ")"
                    << " for center=(" << center_x << ", " << center_y << ")"
                    << " radius=" << radius;
            }
        }
    }
}

//
// Classification correctness
//

TEST(HRGG2DGeometry, ClassificationMatchesExactCircleRectangleRelation) {
    GeometryHarness harness;

    const SInt cells_per_dim = TotalCellsPerDimension(harness.policy);

    const double test_cases[][3] = {
        {0.50, 0.50, 0.10}, {0.31, 0.67, 0.08}, {0.02, 0.02, 0.12},
        {0.98, 0.50, 0.07}, {0.50, 0.98, 0.09}, {0.123, 0.456, 0.035},
    };

    for (const auto& test: test_cases) {
        const double center_x = test[0];
        const double center_y = test[1];
        const double radius   = test[2];

        const HyperRGG2DPolicy::Center center{
            .x          = center_x,
            .y          = center_y,
            .radius     = radius,
            .sampled_id = 0,
            .chunk_id   = 0,
            .cell_id    = 0,
        };

        std::vector<HyperRGG2DPolicy::Cell> candidates;

        harness.policy.CandidateCells(center, radius, candidates);

        for (const auto& cell: candidates) {
            const CellBounds bounds = BoundsOf(cell, cells_per_dim);

            const CellBallRelation expected = ExactRelation(center_x, center_y, radius, bounds);

            const CellBallRelation actual = harness.policy.ClassifyCell(center, radius, cell);

            EXPECT_EQ(actual, expected) << "Incorrect classification for cell "
                                        << "(" << cell.global_cell_x << ", " << cell.global_cell_y << ")"
                                        << " center=(" << center_x << ", " << center_y << ")"
                                        << " radius=" << radius;
        }
    }
}

TEST(HRGG2DGeometry, TangentCellIsNotClassifiedInside) {
    GeometryHarness harness;

    const SInt cells_per_dim = TotalCellsPerDimension(harness.policy);

    const double h = 1.0 / static_cast<double>(cells_per_dim);

    const SInt x = cells_per_dim / 2;

    const SInt y = cells_per_dim / 2;

    //
    // Center exactly one radius left of the neighboring cell edge.
    //
    const double center_x = static_cast<double>(x + 1) * h - 0.5 * h;

    const double center_y = (static_cast<double>(y) + 0.5) * h;

    const double radius = 0.5 * h;

    const HyperRGG2DPolicy::Center center{
        .x          = center_x,
        .y          = center_y,
        .radius     = radius,
        .sampled_id = 0,
        .chunk_id   = 0,
        .cell_id    = 0,
    };

    std::vector<HyperRGG2DPolicy::Cell> candidates;

    harness.policy.CandidateCells(center, radius, candidates);

    const auto* neighboring_cell = FindCell(candidates, x + 1, y);

    ASSERT_NE(neighboring_cell, nullptr);

    EXPECT_NE(harness.policy.ClassifyCell(center, radius, *neighboring_cell), CellBallRelation::INSIDE);
}

//
// Coverage
//

TEST(HRGG2DGeometry, CoverageMatchesExactCircleRectangleIntersection) {
    GeometryHarness harness;

    const SInt cells_per_dim = TotalCellsPerDimension(harness.policy);

    const double h = 1.0 / static_cast<double>(cells_per_dim);

    struct CoverageCase {
        double center_x;
        double center_y;
        double radius;

        SInt cell_x;
        SInt cell_y;
    };

    const SInt center_cell_x = cells_per_dim / 2;

    const SInt center_cell_y = cells_per_dim / 2;

    const double cell_center_x = (static_cast<double>(center_cell_x) + 0.5) * h;

    const double cell_center_y = (static_cast<double>(center_cell_y) + 0.5) * h;

    const CoverageCase cases[] = {
        //
        // Circle centered inside the cell, partial coverage.
        //
        {
            cell_center_x,
            cell_center_y,
            0.6 * h,
            center_cell_x,
            center_cell_y,
        },

        //
        // Circle center exactly on a grid corner.
        //
        {
            0.5,
            0.5,
            h,
            center_cell_x,
            center_cell_y,
        },

        //
        // Off-center boundary crossing.
        //
        {
            cell_center_x + (0.23 * h),
            cell_center_y - (0.17 * h),
            0.73 * h,
            center_cell_x,
            center_cell_y,
        },

        //
        // Very small intersection with the right-hand neighbor.
        //
        {
            cell_center_x,
            cell_center_y,
            0.52 * h,
            center_cell_x + 1,
            center_cell_y,
        },

        //
        // Almost all of the center cell covered.
        //
        {
            cell_center_x,
            cell_center_y,
            0.70 * h,
            center_cell_x,
            center_cell_y,
        },
    };

    for (const auto& test: cases) {
        const HyperRGG2DPolicy::Center center{
            .x          = test.center_x,
            .y          = test.center_y,
            .radius     = test.radius,
            .sampled_id = 0,
            .chunk_id   = 0,
            .cell_id    = 0,
        };

        std::vector<HyperRGG2DPolicy::Cell> candidates;

        harness.policy.CandidateCells(center, test.radius, candidates);

        const auto* cell = FindCell(candidates, test.cell_x, test.cell_y);

        ASSERT_NE(cell, nullptr) << "cell=(" << test.cell_x << ", " << test.cell_y << ")"
                                 << " center=(" << test.center_x << ", " << test.center_y << ")"
                                 << " radius=" << test.radius;

        const CellBounds bounds = BoundsOf(*cell, cells_per_dim);

        const double reference = ExactCircleRectCoverage(test.center_x, test.center_y, test.radius, bounds);

        const double estimated = harness.policy.CellCoverage(center, test.radius, *cell);

        EXPECT_NEAR(estimated, reference, 0.01)
            << "Coverage mismatch"
            << " cell=(" << test.cell_x << ", " << test.cell_y << ")"
            << " center=(" << test.center_x << ", " << test.center_y << ")"
            << " radius=" << test.radius << " reference=" << reference << " estimated=" << estimated;
    }
}

TEST(HRGG2DGeometry, CoverageHandlesExactOutsideAndInsideCases) {
    GeometryHarness harness;

    const SInt cells_per_dim = TotalCellsPerDimension(harness.policy);

    const double h = 1.0 / static_cast<double>(cells_per_dim);

    const SInt x = cells_per_dim / 2;

    const SInt y = cells_per_dim / 2;

    //
    // Fully inside.
    //
    {
        const double center_x = (static_cast<double>(x) + 0.5) * h;

        const double center_y = (static_cast<double>(y) + 0.5) * h;

        const double radius = 0.8 * h;

        const HyperRGG2DPolicy::Center center{
            .x          = center_x,
            .y          = center_y,
            .radius     = radius,
            .sampled_id = 0,
            .chunk_id   = 0,
            .cell_id    = 0,
        };

        std::vector<HyperRGG2DPolicy::Cell> cells;

        harness.policy.CandidateCells(center, radius, cells);

        const auto* cell = FindCell(cells, x, y);

        ASSERT_NE(cell, nullptr);

        EXPECT_DOUBLE_EQ(harness.policy.CellCoverage(center, radius, *cell), 1.0);
    }

    //
    // Fully outside but still present in candidate bounding box.
    //
    {
        const double center_x = (static_cast<double>(x) + 0.5) * h;

        const double center_y = (static_cast<double>(y) + 0.5) * h;

        const double radius = 0.6 * h;

        const HyperRGG2DPolicy::Center center{
            .x          = center_x,
            .y          = center_y,
            .radius     = radius,
            .sampled_id = 0,
            .chunk_id   = 0,
            .cell_id    = 0,
        };

        std::vector<HyperRGG2DPolicy::Cell> cells;

        harness.policy.CandidateCells(center, radius, cells);

        const auto* diagonal = FindCell(cells, x + 1, y + 1);

        ASSERT_NE(diagonal, nullptr);

        EXPECT_DOUBLE_EQ(harness.policy.CellCoverage(center, radius, *diagonal), 0.0);
    }
}

TEST(HRGG2DGeometry, QuarterCircleCoverageHasExpectedValue) {
    GeometryHarness harness;

    const SInt cells_per_dim = TotalCellsPerDimension(harness.policy);

    const double h = 1.0 / static_cast<double>(cells_per_dim);

    //
    // Cell (0,0) is [0,h] x [0,h].
    //
    // A radius-h circle centered at (0,0) covers exactly
    // one quarter of a circle inside that cell.
    //
    const HyperRGG2DPolicy::Center center{
        .x          = 0.0,
        .y          = 0.0,
        .radius     = h,
        .sampled_id = 0,
        .chunk_id   = 0,
        .cell_id    = 0,
    };

    std::vector<HyperRGG2DPolicy::Cell> cells;

    harness.policy.CandidateCells(center, h, cells);

    const auto* cell = FindCell(cells, 0, 0);

    ASSERT_NE(cell, nullptr);

    const double estimated = harness.policy.CellCoverage(center, h, *cell);

    EXPECT_NEAR(estimated, M_PI / 4.0, 0.01);
}