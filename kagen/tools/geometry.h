/*******************************************************************************
 * include/tools/geometry.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
#pragma once

#include "kagen/kagen.h"

#include <cmath>

#ifdef KAGEN_CGAL_FOUND
    #include <CGAL/Dimension.h>
#endif // KAGEN_CGAL_FOUND

namespace kagen {
template <typename Double>
class PGGeometry {
public:
    static inline Double HyperbolicAreaToRadius(const Double area) {
        return std::acosh(area / (2 * M_PI) + 1);
    }

    static bool TestTargetRadius(const Double n, const Double m, const Double alpha, const Double epsilon = 0.01) {
        const Double r = GetTargetRadius(n, m, alpha, epsilon);
        return std::abs(GetExpectedDegree(n, alpha, r) - (m / n) * 2) <= epsilon;
    }

    static Double GetTargetRadius(const Double n, const Double m, const Double alpha, const Double epsilon = 0.01) {
        const Double pl_exp            = 2 * alpha + 1;
        const Double target_avg_degree = (m / n) * 2;
        const Double xi_inv            = ((pl_exp - 2) / (pl_exp - 1));
        const Double v                 = target_avg_degree * (M_PI / 2) * xi_inv * xi_inv;
        Double       current_r         = 2 * std::log(n / v);
        Double       lower_bound       = current_r / 2;
        Double       upper_bound       = current_r * 2;
        do {
            current_r        = (lower_bound + upper_bound) / 2;
            Double current_k = GetExpectedDegree(n, alpha, current_r);
            if (current_k < target_avg_degree) {
                upper_bound = current_r;
            } else {
                lower_bound = current_r;
            }
        } while (std::abs(upper_bound - lower_bound) > epsilon // converge if target_avg_degree is infeasible
                 && std::abs(GetExpectedDegree(n, alpha, current_r) - target_avg_degree) > epsilon);
        return current_r;
    }

    static Double GetExpectedDegree(const Double n, const Double alpha, const Double R) {
        const Double gamma          = 2 * alpha + 1;
        const Double xi             = (gamma - 1) / (gamma - 2);
        const Double first_sum_term = std::exp(-R / 2);
        const Double second_sum_term =
            std::exp(-alpha * R)
            * (alpha * (R / 2) * ((M_PI / 4) * std::pow((1 / alpha), 2) - (M_PI - 1) * (1 / alpha) + (M_PI - 2)) - 1);
        const Double exp_degree = (2 / M_PI) * xi * xi * n * (first_sum_term + second_sum_term);
        return exp_degree;
    }

    static Double EuclideanRadiusToHyperbolic(const Double euclidean_radius) {
        const Double eusq   = euclidean_radius * euclidean_radius;
        const Double result = std::acosh(1 + 2 * eusq / ((1 - eusq)));
        return result;
    }

    static Double HyperbolicRadiusToEuclidean(const Double hyperbolic_radius) {
        const Double ch = std::cosh(hyperbolic_radius);
        return std::sqrt((ch - 1) / (ch + 1));
    }

    static std::pair<Double, Double> PolarToCartesian(const Double phi, const Double r) {
        return std::make_pair(r * std::cos(phi), r * std::sin(phi));
    }

    template <typename HybVertex>
    static Double HyperbolicDistance(const HybVertex& v1, const HybVertex& v2) {
        const Double x1      = std::get<2>(v1);
        const Double x2      = std::get<2>(v2);
        const Double y1      = std::get<3>(v1);
        const Double y2      = std::get<3>(v2);
        const Double gamma1  = std::get<4>(v1);
        const Double gamma2  = std::get<4>(v2);
        const Double delta_x = x1 - x2;
        const Double delta_y = y1 - y2;
        return (delta_x * delta_x + delta_y * delta_y) * gamma1 * gamma2;
    }

    static Double HyperbolicDistance(const Double r1, const Double r2, const Double phi1, const Double phi2) {
        Double result;
        if (phi1 == phi2) {
            result = std::abs(r1 - r2);
        } else {
            const Double delta_phi = M_PI - std::abs(M_PI - std::abs(phi1 - phi2));
            const Double cosh_dist =
                std::cosh(r1) * std::cosh(r2) - std::sinh(r1) * std::sinh(r2) * std::cos(delta_phi);
            result = (cosh_dist >= 1) ? std::acosh(cosh_dist) : 0;
        }
        return result;
    }

    template <typename EucVertex>
    static inline Double SquaredEuclideanDistance(const EucVertex& v1, const EucVertex& v2) {
        const Double x = std::get<0>(v1) - std::get<0>(v2);
        const Double y = std::get<1>(v1) - std::get<1>(v2);
        return x * x + y * y;
    }

    static void CartesianToPolar(const std::pair<Double, Double> a, Double& phi, Double& r) {
        r = sqrt(a.first * a.first + a.second * a.second);
        if (r == 0) {
            phi = 0;
        } else if (a.second >= 0) {
            phi = std::acos(a.first / r);
        } else {
            phi = -std::acos(a.first / r);
        }

        if (phi < 0) {
            phi += 2 * M_PI;
        }
    }

    static void GetEuclideanCircle(
        const std::pair<Double, Double> hyperbolic_center, const Double hyperbolic_radius,
        std::pair<Double, Double>& euclidean_center, Double& euclidean_radius) {
        Double phi_h, r_h;
        PGGeometry::CartesianToPolar(hyperbolic_center, phi_h, r_h);
        Double r_c;
        PGGeometry::GetEuclideanCircle(r_h, hyperbolic_radius, r_c, euclidean_radius);
        euclidean_center = PGGeometry::PolarToCartesian(phi_h, r_c);
    }

    static void GetEuclideanCircle(
        const Double r_h, const Double hyperbolic_radius, Double& euclidean_center_r, Double& euclidean_radius) {
        const Double a     = std::cosh(hyperbolic_radius) - 1;
        const Double b     = 1 - (r_h * r_h);
        euclidean_center_r = (2 * r_h) / (b * a + 2);
        euclidean_radius   = std::sqrt(euclidean_center_r * euclidean_center_r - (2 * r_h * r_h - b * a) / (b * a + 2));
    }

    static inline Double RadiusToHyperbolicArea(const Double radius) {
        return 2 * M_PI * (std::cosh(radius) - 1);
    }
};

#ifdef KAGEN_CGAL_FOUND
/* Tests whehter a sphere intersects with the box */
template <typename BOX, typename SPHERE>
static bool boxIntersects(const BOX& box, const SPHERE& sphere) {
    const auto D  = CGAL::Ambient_dimension<BOX>::value;
    auto       r2 = sphere.squared_radius();
    auto       r  = std::sqrt(r2);

    double dist = 0;
    for (SInt d = 0; d < D; ++d) {
        double e =
            std::max(box.min(d) - sphere.center()[d], (double)0) + std::max(sphere.center()[d] - box.max(d), (double)0);

        if (r < e)
            return false;
        dist += e * e;
    }

    if (dist <= r2)
        return true;

    return false;
}

/* tests whether sphere is FULLY contained in box */
template <typename BOX, typename SPHERE>
static bool boxContains(const BOX& box, const SPHERE& sphere) {
    const auto D  = CGAL::Ambient_dimension<BOX>::value;
    auto       r2 = sphere.squared_radius();

    // check whether center is in box
    for (SInt d = 0; d < D; ++d) {
        if (!(box.min(d) <= sphere.center()[d] && sphere.center()[d] <= box.max(d))) {
            return false;
        }
    }

    //  return boxIntersects(box, sphere);

    // the center of the sphere is within the box
    for (SInt d = 0; d < D; ++d) {
        // project p to the boundary of box in dimension d closest to center of
        // the sphere
        double dist =
            sphere.center()[d] - (sphere.center()[d] < (box.max(d) + box.min(d)) / 2 ? box.min(d) : box.max(d));
        dist = dist * dist;

        if (dist < r2)
            return false;
    }

    return true;
}
#endif // KAGEN_CGAL_FOUND
} // namespace kagen
