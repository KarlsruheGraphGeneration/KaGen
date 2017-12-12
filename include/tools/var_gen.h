#ifndef _VAR_GEN_H_
#define _VAR_GEN_H_

#include <cstddef>

#include <gmp.h>
#include <mpfr.h>

#include "definitions.h"

template <typename Float = HPFloat>
class VarGen {
 public:
  VarGen(SInt seed) {
    initialized_ = false;

    bino_n_last_ = -1.0;
    bino_p_last_ = -1.0;

    bino_h_ = -1.0;
    bino_a_ = -1.0;
    bino_bound_ = -1.0;

    bino_mode_ = -1.0;
    bino_r1_ = -1.0;
    bino_g_ = -1.0;

    hyp_n_last_ = -1.0;
    hyp_m_last_ = -1.0;
    hyp_N_last_ = -1.0;

    hyp_h_ = -1.0;
    hyp_a_ = -1.0;
    hyp_fm_ = -1.0;
    hyp_bound_ = -1.0;

    srand(seed);
    // gmp_randinit_mt(rng);
    // gmp_randseed_ui(rng, seed);
  }

  void RandomInit(SInt seed) {
    srand(seed);
    // gmp_randseed_ui(rng, seed);
  }

  Float Binomial(Float n, double p) {
    SInt inv = 0;  // invert
    Float x;       // result

    // TODO: p==1 float inaccuracy might lead to deadlock
    Float eps = 0.000001;
    if (n == 0 || std::abs(p - 1) < eps) return n;
    if (p == 0) return 0.0;

    if (p > 0.5) {  // faster calculation by inversion
      p = 1. - p;
      inv = 1;
    }

    Float np = n * p;
    if (np < 35.) {
      if (np < 1.E-6) {
        // Poisson approximation for extremely low np
        PoissonLow(np, x);
      } else {
        // inversion method, using chop-down search from 0
        BinomialInver(n, p, x);
      }
    } else {
      // ratio of uniforms method
      BinomialRatioOfUniforms(n, p, x);
    }

    if (inv) {
      x = n - x;  // undo inversion
    }
    // TODO: hotfix
    if (x > n) x = n;
    if (x < 0) x = 0;
    return x;
  }

  Float Hypergeometric(Float n, Float m, Float N) {
    Float fak, addd;
    Float x;

    fak = 1;
    addd = 0;
    if (m > N / 2) {
      // invert m
      m = N - m;
      fak = -1;
      addd = n;
    }
    if (n > N / 2) {
      // invert n
      n = N - n;
      addd += fak * m;
      fak = -fak;
    }
    if (n > m) {
      // swap n and m
      x = n;
      n = m;
      m = x;
    }

    if (N > 680 || n > 70) {
      // use ratio-of-uniforms method
      HypRatioOfUnifoms(n, m, N, x);
    } else {
      // inversion method, using chop-down search from mode
      HypInversionMod(n, m, N, x);
    }

    return x * fak + addd;
  }

 private:
  Float bino_n_last_, bino_p_last_;
  Float bino_h_, bino_a_, bino_bound_;
  Float bino_mode_, bino_r1_, bino_g_;

  Float hyp_n_last_, hyp_m_last_, hyp_N_last_;
  Float hyp_h_, hyp_a_, hyp_fm_, hyp_bound_;
  Float hyp_mode_, hyp_mp_;

  // gmp_randstate_t rng;

  bool initialized_;
  static const SInt kFakLen = 1024;
  const double C0 = 0.918938533204672722;
  const double C1 = 1. / 12.;
  const double C3 = -1. / 360.;
  Float fac_table[kFakLen];

  const double SHAT1 = 2.943035529371538573;   // 8/e
  const double SHAT2 = 0.8989161620588987408;  // 3-sqrt(12/e)

  Float LnFac(Float n) {
    // log factorial function. gives natural logarithm of n!

    if (n < (Float)kFakLen) {
      if (n <= 1) {
        return 0;
      }
      if (!initialized_) {  // first time. Must initialize table
                            // make table of ln(n!)
        Float sum = fac_table[0] = 0.;
        for (SInt i = 1; i < kFakLen; i++) {
          sum += log(Float(i));
          fac_table[i] = sum;
        }
        initialized_ = 1;
      }
      return fac_table[(SInt)n];
    }
    // not found in table. use Stirling approximation
    // float128  n1, r;
    Float n1 = n;
    Float r = 1. / n1;
    return (n1 + 0.5) * log(n1) - n1 + C0 + r * (C1 + r * r * C3);
  }

  void fc_lnpk(Float& k, Float& L, Float& m, Float& n, Float& result) {
    result = LnFac(k) + LnFac(m - k) + LnFac(n - k) + LnFac(L + k);
  }

  void PoissonLow(double L, Float& result) {
    double d, r;
    d = sqrt(L);
    if (rand() / double(RAND_MAX) >= d) {
      result = 0;
      return;
    }
    r = (rand() / double(RAND_MAX)) * d;
    if (r > L * (1. - L)) {
      result = 0;
      return;
    }
    if (r > 0.5 * L * L * (1. - L)) {
      result = 1;
      return;
    }
    result = 2;
    return;
  }

  void BinomialInver(Float& n, double p, Float& result) {
    double f0, f, q;
    SInt bound;
    double pn, r, rc;
    SInt x, n1, i;

    // f(0) = probability of x=0 is (1-p)^n
    // fast calculation of (1-p)^n
    f0 = 1.;
    pn = 1. - p;
    n1 = n;
    while (n1) {
      if (n1 & 1) f0 *= pn;
      pn *= pn;
      n1 >>= 1;
    }
    // calculate safety bound
    rc = (n + 1) * p;
    bound = (SInt)(rc + 11.0 * (sqrt(rc) + 1.0));
    if (bound > n) bound = n;
    q = p / (1. - p);

    while (1) {
      r = rand() / double(RAND_MAX);
      // recursive calculation: f(x) = f(x-1) * (n-x+1)/x*p/(1-p)
      f = f0;
      x = 0;
      i = n;
      do {
        r -= f;
        if (r <= 0) {
          result = x;
          return;
        }
        x++;
        f *= q * i;
        r *= x;  // it is faster to multiply r by x than dividing f by x
        i--;
      } while (x <= bound);
    }
  }

  void BinomialRatioOfUniforms(Float& n, double p, Float& result) {
    /*
    Subfunction for Binomial distribution. Assumes p < 0.5.

    Uses the Ratio-of-Uniforms rejection method.

    The computation time hardly depends on the parameters, except that it
    matters a lot whether parameters are within the range where the LnFac
    function is tabulated.

    Reference: E. Stadlober: "The ratio of uniforms approach for generating
    discrete random variates". Journal of Computational and Applied Mathematics,
    vol. 31, no. 1, 1990, pp. 181-189.
    */
    Float u;    // uniform random
    Float q1;   // 1-p
    Float np;   // n*p
    Float var;  // variance
    Float lf;   // ln(f(x))
    Float x;    // real sample
    Float k;    // integer sample

    if (bino_n_last_ != n || bino_p_last_ != p) {  // Set_up
      bino_n_last_ = n;
      bino_p_last_ = p;
      q1 = 1.0 - p;
      np = n * p;
      bino_mode_ = floor(np + p);  // mode
      bino_a_ = np + 0.5;          // hat center
      bino_r1_ = log(p / q1);
      bino_g_ = LnFac(bino_mode_) + LnFac(n - bino_mode_);
      var = np * q1;                                 // variance
      bino_h_ = sqrt(SHAT1 * (var + 0.5)) + SHAT2;   // hat width
      bino_bound_ = floor(bino_a_ + 6.0 * bino_h_);  // safety-bound
      if (bino_bound_ > n) bino_bound_ = n;          // safety-bound
    }

    while (1) {  // rejection loop
      // mpfr_urandomb(u.mpfr_ptr(), rng);
      u = rand() / double(RAND_MAX);
      if (u == 0) continue;  // avoid division by 0
      // mpfr_urandomb(x.mpfr_ptr(), rng);
      x = rand() / double(RAND_MAX);
      x = bino_a_ + bino_h_ * (x - 0.5) / u;
      if (x < 0.) continue;  // reject, avoid overflow
      k = floor(x);          // truncate
      lf = (k - bino_mode_) * bino_r1_ + bino_g_ - LnFac(k) -
           LnFac(n - k);                     // ln(f(k))
      if (u * (4.0 - u) - 3.0 <= lf) break;  // lower squeeze accept
      if (u * (u - lf) > 1.0) continue;      // upper squeeze reject
      if (2.0 * log(u) <= lf) break;         // final acceptance
    }

    // return k;
    result = k;
  }

  void HypInversionMod(SInt n, SInt m, SInt N, Float& result) {
    /*
    Subfunction for Hypergeometric distribution. Assumes 0 <= n <= m <= N/2.
    Overflow protection is needed when N > 680 or n > 75.

    Hypergeometric distribution by inversion method, using down-up
    search starting at the mode using the chop-down technique.

    This method is faster than the rejection method when the variance is low.
    */

    // Sampling
    SInt I;              // Loop counter
    SInt L = N - m - n;  // Parameter
    double modef;        // mode, float
    double Mp, np;       // m + 1, n + 1
    double p;            // temporary
    double U;            // uniform random
    double c, d;         // factors in iteration
    double divisor;      // divisor, eliminated by scaling
    double k1, k2;       // float version of loop counter
    double L1 = L;       // float version of L

    Mp = (double)(m + 1);
    np = (double)(n + 1);

    if (N != hyp_N_last_ || m != hyp_m_last_ || n != hyp_n_last_) {
      // set-up when parameters have changed
      hyp_N_last_ = N;
      hyp_m_last_ = m;
      hyp_n_last_ = n;

      p = Mp / (N + 2.);
      modef = np * p;           // mode, real
      hyp_mode_ = (SInt)modef;  // mode, integer
      if (hyp_mode_ == modef && p == 0.5) {
        hyp_mp_ = hyp_mode_--;
      } else {
        hyp_mp_ = hyp_mode_ + 1;
      }
      // mode probability, using log factorial function
      // (may read directly from fac_table if N < FAK_LEN)
      hyp_fm_ = exp(LnFac(N - m) - LnFac(L + hyp_mode_) - LnFac(n - hyp_mode_) +
                    LnFac(m) - LnFac(m - hyp_mode_) - LnFac(hyp_mode_) -
                    LnFac(N) + LnFac(N - n) + LnFac(n));

      // safety bound - guarantees at least 17 significant decimal digits
      // bound = min(n, (SInt)(modef + k*c'))
      hyp_bound_ = (SInt)(
          modef + 11. * sqrt(modef * (1. - p) * (1. - n / (double)N) + 1.));
      if (hyp_bound_ > n) hyp_bound_ = n;
    }

    // loop until accepted
    while (1) {
      U = rand() / double(RAND_MAX);  // uniform random number to be converted

      // start chop-down search at mode
      if ((U -= hyp_fm_) <= 0.) {
        result = hyp_mode_;
        return;
      }
      c = d = hyp_fm_;

      // alternating down- and upward search from the mode
      k1 = hyp_mp_ - 1;
      k2 = hyp_mode_ + 1;
      for (I = 1; I <= hyp_mode_; I++, k1--, k2++) {
        // Downward search from k1 = hyp_mp_ - 1
        divisor = (np - k1) * (Mp - k1);
        // Instead of dividing c with divisor, we multiply U and d because
        // multiplication is faster. This will give overflow if N > 800
        U *= divisor;
        d *= divisor;
        c *= k1 * (L1 + k1);
        if ((U -= c) <= 0.) {
          result = hyp_mp_ - I - 1;  // = k1 - 1
          return;
        }

        // Upward search from k2 = hyp_mode_ + 1
        divisor = k2 * (L1 + k2);
        // re-scale parameters to avoid time-consuming division
        U *= divisor;
        c *= divisor;
        d *= (np - k2) * (Mp - k2);
        if ((U -= d) <= 0.) {
          result = hyp_mode_ + I;
          return;
        }
        // Values of n > 75 or N > 680 may give overflow if you leave out this..
        // overflow protection
        // if (U > 1.E100) {U *= 1.E-100; c *= 1.E-100; d *= 1.E-100;}
      }

      // Upward search from k2 = 2*mode + 1 to bound
      for (k2 = I = hyp_mp_ + hyp_mode_; I <= hyp_bound_; I++, k2++) {
        divisor = k2 * (L1 + k2);
        U *= divisor;
        d *= (np - k2) * (Mp - k2);
        if ((U -= d) <= 0.) {
          result = I;
          return;
        }
        // more overflow protection
        // if (U > 1.E100) {U *= 1.E-100; d *= 1.E-100;}
      }
    }
  }

  void HypRatioOfUnifoms(Float& n, Float& m, Float& N, Float& result) {
    /*
    Subfunction for Hypergeometric distribution using the ratio-of-uniforms
    rejection method.

    This code is valid for 0 < n <= m <= N/2.

    The computation time hardly depends on the parameters, except that it
    matters a lot whether parameters are within the range where the LnFac
    function is tabulated.

    Reference: E. Stadlober: "The ratio of uniforms approach for generating
    discrete random variates". Journal of Computational and Applied Mathematics,
    vol. 31, no. 1, 1990, pp. 181-189.
    */

    // convert inputs
    Float L;     // N-m-n
    Float mode;  // mode
    Float k;     // integer sample
    Float x;     // real sample
    Float rNN;   // 1/(N*(N+2))
    Float my;    // mean
    Float var;   // variance
    Float u;     // uniform random
    Float lf;    // ln(f(x))
    Float lnfac;

    L = N - m - n;
    if (hyp_N_last_ != N || hyp_m_last_ != m || hyp_n_last_ != n) {
      hyp_N_last_ = N;
      hyp_m_last_ = m;
      hyp_n_last_ = n;             // Set-up
      rNN = 1. / (N * (N + 2));    // make two divisions in one
      my = n * m * rNN * (N + 2);  // mean = n*m/N
      mode = floor((n + 1) * (m + 1) * rNN *
                   N);  // mode = floor((n+1)*(m+1)/(N+2))
      var = n * m * (N - m) * (N - n) / (N * N * (N - 1));  // variance
      hyp_h_ = sqrt(SHAT1 * (var + 0.5)) + SHAT2;           // hat width
      hyp_a_ = my + 0.5;                                    // hat center
      fc_lnpk(mode, L, m, n, lnfac);
      hyp_fm_ = lnfac;
      hyp_bound_ = floor(hyp_a_ + 4.0 * hyp_h_);  // safety-bound
      if (hyp_bound_ > n) hyp_bound_ = n;
    }

    while (1) {
      // mpfr_urandomb(u.mpfr_ptr(), rng);
      u = rand() / double(RAND_MAX);
      if (u == 0) continue;  // avoid division by 0
      // mpfr_urandomb(x.mpfr_ptr(), rng);
      x = rand() / double(RAND_MAX);
      x = hyp_a_ + hyp_h_ * (x - 0.5) / u;  // generate hat distribution
      if (x < 0.) continue;                 // reject, avoid overflow
      k = floor(x);
      if (k > hyp_bound_) continue;  // reject if outside range
      // lf = hyp_fm_ - fc_lnpk(k,L,m,n);          // ln(f(k))
      fc_lnpk(k, L, m, n, lnfac);
      lf = hyp_fm_ - lnfac;                  // ln(f(k))
      if (u * (4.0 - u) - 3.0 <= lf) break;  // lower squeeze accept
      if (u * (u - lf) > 1.0) continue;      // upper squeeze reject
      if (2.0 * log(u) <= lf) break;         // final acceptance
    }

    // return k;
    result = k;
  }
};

#endif
