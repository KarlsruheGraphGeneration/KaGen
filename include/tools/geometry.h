/*******************************************************************************
 * include/tools/geometry.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include <cmath>
#include <CGAL/Dimension.h>

namespace kagen {

class PGGeometry {
 public:
  PGGeometry();
  virtual ~PGGeometry();

  static inline double HyperbolicAreaToRadius(double area) {
    return acosh(area / (2 * M_PI) + 1);
  }

  static double GetTargetRadius(double n, double m, double alpha = 1,
                                double /* T */ = 0, double epsilon = 0.01) {
    double pl_exp = 2 * alpha + 1;
    double target_avg_degree = (m / n) * 2;
    double xi_inv = ((pl_exp - 2) / (pl_exp - 1));
    double v = target_avg_degree * (M_PI / 2) * xi_inv * xi_inv;
    double current_r = 2 * log(n / v);
    double lower_bound = current_r / 2;
    double upper_bound = current_r * 2;
    do {
      current_r = (lower_bound + upper_bound) / 2;
      double current_k = getExpectedDegree(n, alpha, current_r);
      if (current_k < target_avg_degree) {
        upper_bound = current_r;
      } else {
        lower_bound = current_r;
      }
    } while (std::abs(getExpectedDegree(n, alpha, current_r) -
                      target_avg_degree) > epsilon);
    return current_r;
  }

  static double getExpectedDegree(double n, double alpha, double R) {
    double gamma = 2 * alpha + 1;
    double xi = (gamma - 1) / (gamma - 2);
    double first_sum_term = exp(-R / 2);
    double second_sum_term =
        exp(-alpha * R) * (alpha * (R / 2) *
                               ((M_PI / 4) * pow((1 / alpha), 2) -
                                (M_PI - 1) * (1 / alpha) + (M_PI - 2)) -
                           1);
    double exp_degree =
        (2 / M_PI) * xi * xi * n * (first_sum_term + second_sum_term);
    return exp_degree;
  }

  static double EuclideanRadiusToHyperbolic(double euclidean_radius) {
    double eusq = euclidean_radius * euclidean_radius;
    double result = acosh(1 + 2 * eusq / ((1 - eusq)));
    return result;
  }

  static double HyperbolicRadiusToEuclidean(double hyperbolic_radius) {
    double ch = cosh(hyperbolic_radius);
    return sqrt((ch - 1) / (ch + 1));
  }

  static std::pair<double, double> PolarToCartesian(double phi, double r) {
    return std::pair<double, double>(r * cos(phi), r * sin(phi));
  }

  template <typename HypVertex>
  static double HyperbolicDistance(const HypVertex &v1, const HypVertex &v2) {
    double x1 = std::get<2>(v1);
    double x2 = std::get<2>(v2);
    double y1 = std::get<3>(v1);
    double y2 = std::get<3>(v2);
    double gamma1 = std::get<4>(v1);
    double gamma2 = std::get<4>(v2);
    double delta_x = x1 - x2;
    double delta_y = y1 - y2;
    // return hyperbolicDistance(std::get<1>(v1), std::get<1>(v2),
    // std::get<0>(v1), std::get<0>(v2));
    return (delta_x * delta_x + delta_y * delta_y) * gamma1 * gamma2;
  }

  static double HyperbolicDistance(double r1, double r2, double phi1,
                                   double phi2) {
    double result;
    if (phi1 == phi2) {
      result = std::abs(r1 - r2);
    } else {
      double delta_phi = M_PI - std::abs(M_PI - std::abs(phi1 - phi2));
      double cosh_dist =
          cosh(r1) * cosh(r2) - sinh(r1) * sinh(r2) * cos(delta_phi);
      if (cosh_dist >= 1)
        result = acosh(cosh_dist);
      else
        result = 0;
    }
    return result;
  }

  template <typename EucVertex>
  static inline double SquaredEuclideanDistance(const EucVertex &v1,
                                                const EucVertex &v2) {
    LPFloat x = std::get<0>(v1) - std::get<0>(v2);
    LPFloat y = std::get<1>(v1) - std::get<1>(v2);
    return x * x + y * y;
  }

  static void CartesianToPolar(std::pair<double, double> a, double &phi,
                               double &r) {
    r = sqrt(a.first * a.first + a.second * a.second);
    if (r == 0)
      phi = 0;
    else if (a.second >= 0) {
      phi = acos(a.first / r);
    } else {
      phi = -acos(a.first / r);
    }
    if (phi < 0) phi += 2 * M_PI;
  }

  static void GetEuclideanCircle(std::pair<double, double> hyperbolic_center,
                                 double hyperbolic_radius,
                                 std::pair<double, double> &euclidean_center,
                                 double &euclidean_radius) {
    double phi_h, r_h;
    PGGeometry::CartesianToPolar(hyperbolic_center, phi_h, r_h);
    double r_c;
    PGGeometry::GetEuclideanCircle(r_h, hyperbolic_radius, r_c, euclidean_radius);
    euclidean_center = PGGeometry::PolarToCartesian(phi_h, r_c);
  }

  static void GetEuclideanCircle(double r_h, double hyperbolic_radius,
                                 double &euclidean_center_r,
                                 double &euclidean_radius) {
    double a = cosh(hyperbolic_radius) - 1;
    double b = 1 - (r_h * r_h);
    euclidean_center_r = (2 * r_h) / (b * a + 2);
    euclidean_radius = sqrt(euclidean_center_r * euclidean_center_r -
                            (2 * r_h * r_h - b * a) / (b * a + 2));
  }

  static inline double RadiusToHyperbolicArea(double radius) {
    return 2 * M_PI * (cosh(radius) - 1);
  }
};

/* Tests whehter a sphere intersects with the box */
template <typename BOX, typename SPHERE>
static bool boxIntersects(const BOX &box, const SPHERE &sphere) {
  const auto D = CGAL::Ambient_dimension<BOX>::value;
  auto r2 = sphere.squared_radius();
  auto r = std::sqrt(r2);

  double dist = 0;
  for (SInt d = 0; d < D; ++d) {
    double e = std::max(box.min(d) - sphere.center()[d], (double)0) +
               std::max(sphere.center()[d] - box.max(d), (double)0);

    if (r < e) return false;
    dist += e * e;
  }

  if (dist <= r2) return true;

  return false;
}

/* tests whether sphere is FULLY contained in box */
template <typename BOX, typename SPHERE>
static bool boxContains(const BOX &box, const SPHERE &sphere) {
  const auto D = CGAL::Ambient_dimension<BOX>::value;
  auto r2 = sphere.squared_radius();

  // check whether center is in box
  for (SInt d = 0; d < D; ++d) {
    if (!(box.min(d) <= sphere.center()[d] &&
          sphere.center()[d] <= box.max(d))) {
      return false;
    }
  }

  //  return boxIntersects(box, sphere);

  // the center of the sphere is within the box
  for (SInt d = 0; d < D; ++d) {
    // project p to the boundary of box in dimension d closest to center of
    // the sphere
    double dist =
        sphere.center()[d] - (sphere.center()[d] < (box.max(d) + box.min(d)) / 2
                                  ? box.min(d)
                                  : box.max(d));
    dist = dist * dist;

    if (dist < r2) return false;
  }

  return true;
}

}
#endif
