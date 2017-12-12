#ifndef _BENCHMARK_H_
#define _BENCHMARK_H_

#include <cassert>
#include <cmath>

// template-based loop unrolling
template <size_t N>
struct FauxUnroll {
  template <typename F>
  static void call(F &&f) {
    FauxUnroll<N - 1>::call(f);
    f(N - 1);
  }
};

template <>
struct FauxUnroll<0> {
  template <typename F>
  static void call(F &&) {}
};

struct Statistics {
  // Single-pass standard deviation calculation as described in Donald Knuth:
  // The Art of Computer Programming, Volume 2, Chapter 4.2.2, Equations 15&16
  double mean_;
  double nvar_;  // approx n * variance; stddev = sqrt(nvar_ / (count_-1))
  size_t count_;

  Statistics() : mean_(0.0), nvar_(0.0), count_(0) {}

  void Push(double t) {
    ++count_;
    if (count_ == 1) {
      mean_ = t;
    } else {
      double oldmean = mean_;
      mean_ += (t - oldmean) / count_;
      nvar_ += (t - oldmean) * (t - mean_);
    }
  }

  double Avg() const { return mean_; }
  double Stddev() const { return count_ > 1 ? sqrt(nvar_ / (count_ - 1)) : 0; }
};

#endif
