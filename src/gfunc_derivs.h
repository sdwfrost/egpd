#ifndef EGPD_GFUNC_DERIVS_H
#define EGPD_GFUNC_DERIVS_H
// Shared helpers for the native EGPD likelihoods: the stabilised continuous-EGPD
// objective (egpd_cont_nll) and its exact Jet version used where the generated
// closed-form derivatives are unreliable, predictor-scale stencils for d34,
// the truncated-beta (carrier 6) normaliser, Gauss-Kronrod constants shared with
// the DEGPD evaluators, and truncated-normal G-function building blocks.
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <functional>
#include <limits>
#include "egpd_jet.h"

// Continuous GPD log survival, including its removable singularity at xi=0.
// The series also keeps numerical differentiation smooth through zero.
static inline double egpd_tail_exponent(double z, double xi) {
  double v = xi * z;
  if (std::abs(v) < 1e-4)
    return z * (1.0 + v * (-0.5 + v * (1.0/3.0 + v * (-0.25 +
      v * (0.2 + v * (-1.0/6.0 + v / 7.0))))));
  return std::log1p(v) / xi;
}

static inline double egpd_logpower(double coefficient, double logx) {
  return coefficient == 0.0 ? 0.0 : coefficient * logx;
}

static inline double egpd_expm1_over_z(double z) {
  if (std::abs(z) < 0.1) {
    double r = 1.0;
    for (int k = 16; k >= 1; --k)
      r = 1.0 + (z / static_cast<double>(k + 1)) * r;
    return r;
  }
  return std::expm1(z) / z;
}

// Exact log G3(1 - exp(-A); exp(ldelta)) = log(1 - exp(-A)*(1 + (1 - exp(-d*A))/d))
// across A -> 0+, d -> 0+, and large d*A (avoiding pbeta underflow to 1.0).
static inline double egpd_log_G3(double A, double logA, double ldelta) {
  if (A <= 0.0) return R_NegInf;
  if (A >= 1.0) {
    double d = std::exp(std::min(ldelta, 700.0));
    double z = -d * A;
    double ratio;
    if (z < -50.0) ratio = -std::exp(-ldelta);
    else if (std::abs(z) < 0.1) ratio = -A * egpd_expm1_over_z(z);
    else ratio = std::expm1(z) / d;
    double logsf = -A + std::log1p(-ratio);
    return logsf < -0.6931471805599453 ? std::log1p(-std::exp(logsf)) : std::log(-std::expm1(logsf));
  }
  double ldu = ldelta + logA;
  if (ldelta <= 40.0) {
    double dp1 = 1.0 + std::exp(ldelta);
    if (dp1 * A < 0.2) {
      double P[17];
      P[0] = 1.0;
      for (int k = 1; k < 17; ++k) P[k] = 1.0 + dp1 * P[k - 1];
      double s = P[16];
      for (int k = 15; k >= 0; --k)
        s = P[k] - (A / static_cast<double>(k + 3)) * s;
      return std::log(dp1) - std::log(2.0) + 2.0 * logA + std::log(s);
    }
  } else if (ldu < -1.6) {
    double du = std::exp(ldu);
    double s = 1.0;
    for (int k = 16; k >= 1; --k)
      s = 1.0 - (du / static_cast<double>(k + 2)) * s;
    return logA + ldu - std::log(2.0) + std::log(s);
  }
  double term2 = (ldu > 50.0) ? std::exp(-ldu) : egpd_expm1_over_z(-std::exp(ldu));
  return -A + logA + std::log(egpd_expm1_over_z(A) - term2);
}

static inline bool egpd_cont_oos(double y, double lpsi, double xi) {
  if (y < 0.0 || !std::isfinite(y) || !std::isfinite(lpsi) || !std::isfinite(xi))
    return true;
  double z = y * std::exp(-lpsi);
  return (xi < 0.0 && 1.0 + xi * z <= 0.0);
}

// Carrier 5 (truncated normal) normaliser Phi(sqrt(kappa)) - 1/2. As kappa -> 0
// the direct 1/2 - pnorm(-sqrt(kappa)) cancels (relative error ~eps/sqrt(kappa);
// no digits are left below kappa ~ 1e-32), so small kappa uses the equal
// erf(sqrt(kappa/2))/2, which involves no subtraction.
static inline double egpd_tn_norm(double kappa) {
  if (kappa < 1e-2) return 0.5 * std::erf(std::sqrt(0.5 * kappa));
  return 0.5 - R::pnorm(-std::sqrt(kappa), 0.0, 1.0, 1, 0);
}

// Carrier 5 with kappa = exp(lkappa) below DBL_MIN (lkappa < -708.4): kappa is
// then subnormal or 0 and keeps only a few significant bits, while 0.5*lkappa
// is exact, so 0.5*lkappa - log(egpd_tn_norm(kappa)) no longer cancelled: the
// objective varied by up to 0.2 per observation (likelihoods far above the
// maximum near lkappa = -743), and was the 1e20 penalty below -744.1. Writing
// Phi(sqrt(kappa)) - 1/2 = sqrt(kappa/(2 pi)) M(kappa), M(kappa) =
// sum_n (-kappa/2)^n / (n! (2n+1)) = 1 - kappa/6 + O(kappa^2), the cancellation
// is exact: log g(v) = -kappa (1 - v)^2 / 2 - log M(kappa), whose kappa terms
// (-log M = kappa/6 here) are themselves below DBL_MIN: the GPD limit.
static inline bool egpd_tn_subnormal(double kappa) {
  return kappa < std::numeric_limits<double>::min();
}

// Carrier 6 log-normaliser log Z(kappa) (see tb_normaliser below).
static inline double egpd_tb_logZ(double lkappa);

// Predictor-scale objective shared by the exact likelihood and the local
// derivative stencils. Models use the native (not public carrier-type) ordering;
// model == -2 uses direct lkappa2 for egpd2d34.
static inline double egpd_cont_nll(double y, const std::vector<double>& t, int model) {
  for (size_t i = 0; i < t.size(); ++i) {
    if (!std::isfinite(t[i])) return 1e20;
  }
  if (egpd_cont_oos(y, t[0], t[1]))
    return 1e20;
  if (y == 0.0) {
    if (model == 1) return (t[2] == 0.0) ? t[0] : 1e20;
    if (model == 2) {
      if (t[2] != 0.0) return 1e20;
      double lp = -std::max(-t[4], 0.0) - std::log1p(std::exp(-std::abs(t[4])));
      return t[0] - lp;
    }
    if (model == -2) {
      if (t[2] == 0.0 && t[3] == 0.0) return t[0];
      if (t[2] == 0.0 && t[3] > 0.0) {
        double lp = -std::max(-t[4], 0.0) - std::log1p(std::exp(-std::abs(t[4])));
        return t[0] - lp;
      }
      if (t[3] == 0.0 && t[2] > 0.0) {
        double lq = -std::max(t[4], 0.0) - std::log1p(std::exp(-std::abs(t[4])));
        return t[0] - lq;
      }
      return 1e20;
    }
    if (model == 3) return 1e20;
    if (model == 4) {
      if (t[3] != 0.0) return 1e20;
      double log1pd = t[2] > 0.0 ? t[2] + std::log1p(std::exp(-t[2])) : std::log1p(std::exp(t[2]));
      return t[0] - 0.5 * (log1pd - std::log(2.0));
    }
    double k = std::exp(t[2]);
    // kappa = 0 (lkappa < -745.1) is the kappa -> 0 limit of carriers 5 and 6
    if (!std::isfinite(k)) return 1e20;
    if (model == 5) {
      if (egpd_tn_subnormal(k)) return t[0] - (-0.5 * k + k / 6.0);
      double denom = egpd_tn_norm(k);
      double logg = 0.5 * t[2] - 0.5 * std::log(2.0 * M_PI) - 0.5 * k - std::log(denom);
      return std::isfinite(logg) ? (t[0] - logg) : 1e20;
    }
    // Carrier 6: the incomplete-beta difference below cancels as kappa -> 0
    // (relative error ~eps/kappa). There use the equal normalised form
    // g = (1 - C s^2)^(kappa - 1) / Z(kappa), C = (15/16)^2, here at s = 1.
    if (k < 1e-2) {
      double logg = (k - 1.0) * std::log1p(-225.0 / 256.0) - egpd_tb_logZ(t[2]);
      return std::isfinite(logg) ? (t[0] - logg) : 1e20;
    }
    double lo = 1.0 / 32.0, width = 0.5 - lo;
    double denom = R::pbeta(0.5, k, k, 1, 0) - R::pbeta(lo, k, k, 1, 0);
    double logg = std::log(width) + R::dbeta(lo, k, k, 1) - std::log(denom);
    return std::isfinite(logg) ? (t[0] - logg) : 1e20;
  }
  double z = y * std::exp(-t[0]), xi = t[1], v = xi * z;
  double r_over_z;
  if (std::abs(v) < 1e-4) {
    r_over_z = 1.0 + v * (-0.5 + v * (1.0/3.0 + v * (-0.25 +
      v * (0.2 + v * (-1.0/6.0 + v / 7.0)))));
  } else {
    r_over_z = std::log1p(v) / v;
  }
  double A = z * r_over_z;
  double logA = std::log(y) - t[0] + std::log(r_over_z);
  double logH = (A < 0.1)
    ? (logA + std::log(egpd_expm1_over_z(-A)))
    : (A > 0.6931471805599453 ? std::log1p(-std::exp(-A)) : std::log(-std::expm1(-A)));
  double k = std::exp(t[2]), logg;
  if (model == 1) {
    logg = t[2] + egpd_logpower(k - 1.0, logH);
  } else if (model == 2 || model == -2) {
    double logd;
    if (model == 2) {
      double mxk = std::max(t[2], t[3]);
      logd = mxk + std::log1p(std::exp(-std::abs(t[2] - t[3])));
    } else {
      logd = t[3];
    }
    double d = std::exp(logd);
    double lp = -std::max(-t[4], 0.0) - std::log1p(std::exp(-std::abs(t[4])));
    double lq = -std::max(t[4], 0.0) - std::log1p(std::exp(-std::abs(t[4])));
    double a = lp + t[2] + egpd_logpower(k - 1.0, logH);
    double b = lq + logd + egpd_logpower(d - 1.0, logH);
    double mx = std::max(a, b);
    logg = std::isinf(mx) ? mx : mx + std::log1p(std::exp(-std::abs(a - b)));
  } else if (model == 3 || model == 4) {
    double ldelta = t[2];
    double ldu = ldelta + logA;
    if (ldu < -2.3) {
      double du = std::exp(ldu);
      double log1pd = ldelta > 0.0 ? ldelta + std::log1p(std::exp(-ldelta)) : std::log1p(std::exp(ldelta));
      logg = log1pd + logA + std::log(egpd_expm1_over_z(-du));
    } else {
      double log1p_inv_d = ldelta > 0.0 ? std::log1p(std::exp(-ldelta)) : std::log1p(std::exp(ldelta)) - ldelta;
      double log_1m_edA;
      if (ldu > 50.0) {
        log_1m_edA = 0.0;
      } else {
        double du = std::exp(ldu);
        log_1m_edA = du > 0.6931471805599453 ? std::log1p(-std::exp(-du)) : std::log(-std::expm1(-du));
      }
      logg = log1p_inv_d + log_1m_edA;
    }
    if (model == 4) {
      k = std::exp(t[3]);
      double logG = egpd_log_G3(A, logA, ldelta);
      logg += t[3] - std::log(2.0) + egpd_logpower(k / 2.0 - 1.0, logG);
    }
  } else if (model == 5) {
    if (egpd_tn_subnormal(k)) {
      logg = -0.5 * k * std::exp(-2.0 * A) + k / 6.0;
    } else {
      double denom = egpd_tn_norm(k);
      logg = 0.5 * t[2] - 0.5 * std::log(2.0 * M_PI) -
        0.5 * k * std::exp(-2.0 * A) - std::log(denom);
    }
  } else if (k < 1e-2) {
    // carrier 6 at small kappa, as at y = 0, with s = exp(-A)
    logg = (k - 1.0) * std::log1p(-(225.0 / 256.0) * std::exp(-2.0 * A)) -
      egpd_tb_logZ(t[2]);
  } else {
    double lo = 1.0 / 32.0, width = 0.5 - lo;
    double u = lo + width * (-std::expm1(-A));
    double denom = R::pbeta(0.5, k, k, 1, 0) - R::pbeta(lo, k, k, 1, 0);
    logg = std::log(width) + R::dbeta(u, k, k, 1) - std::log(denom);
  }
  double res = t[0] + (1.0 + xi) * A - logg;
  return std::isfinite(res) ? res : 1e20;
}

// The generated analytic expressions divide by xi repeatedly. In a small,
// explicit neighbourhood of zero use centred predictor-scale stencils of the
// same smooth objective, rather than perturbing or clipping the model's shape.
static inline bool egpd_cont_near_zero(double y, const std::vector<double>& t) {
  return std::abs(t[1]) * std::max(1.0, y * std::exp(-t[0])) < 1e-3;
}

static inline double egpd_cont_partial(double y, std::vector<double>& t, int model,
                                      const std::vector<int>& orders, int axis,
                                      const std::vector<double>& steps) {
  static const double weights[5][5] = {
    {0, 0, 1, 0, 0}, {1.0/12, -2.0/3, 0, 2.0/3, -1.0/12},
    {-1.0/12, 4.0/3, -2.5, 4.0/3, -1.0/12},
    {-0.5, 1, 0, -1, 0.5}, {1, -4, 6, -4, 1}};
  if (axis == static_cast<int>(t.size())) {
    // Never difference the out-of-support penalty: flag the stencil instead.
    double f = egpd_cont_nll(y, t, model);
    return f >= 1e19 ? R_NaN : f;
  }
  int order = orders[axis];
  if (order == 0) return egpd_cont_partial(y, t, model, orders, axis + 1, steps);
  double sum = 0.0, original = t[axis];
  for (int j = 0; j < 5; ++j) {
    if (weights[order][j] == 0.0) continue;
    t[axis] = original + (j - 2) * steps[axis];
    sum += weights[order][j] * egpd_cont_partial(y, t, model, orders, axis + 1, steps);
  }
  t[axis] = original;
  return sum / std::pow(steps[axis], order);
}

static inline arma::rowvec egpd_cont_derivs(double y, std::vector<double> t,
                                           int model, bool higher = false) {
  int n = t.size();
  std::vector<double> steps(n, higher ? 0.01 : 0.0002);
  double z = y * std::exp(-t[0]), v = t[1] * z, g = 1.0 + v;
  steps[1] /= std::max(1.0, z);
  // Near the xi < 0 endpoint the objective varies on the scale of its distance
  // g = 1 + xi*z to it. Keep the whole stencil (+-2 steps) inside 2% of that
  // distance; powers of two keep the shrunken offsets exactly representable.
  const double frac = 0.02;
  double h0 = std::abs(v) > 0.0 ? frac * g / (4.0 * std::abs(v)) : steps[0];
  double h1 = z > 0.0 ? frac * g / (4.0 * z) : steps[1];
  if (h0 < steps[0]) steps[0] = std::ldexp(1.0, std::ilogb(h0));
  if (h1 < steps[1]) steps[1] = std::ldexp(1.0, std::ilogb(h1));
  for (int i = 0; i < 2; ++i) {
    // Offsets below ~64 ulps of the coordinate cannot resolve derivatives.
    if (!(steps[i] > 64.0 * std::numeric_limits<double>::epsilon() *
          std::abs(t[i])) || !(steps[i] > 0.0)) {
      int nout = higher
        ? n * (n + 1) * (n + 2) / 6 + n * (n + 1) * (n + 2) * (n + 3) / 24
        : n + n * (n + 1) / 2;
      arma::rowvec bad(nout);
      bad.fill(R_NaN);
      return bad;
    }
  }
  std::vector<int> orders(n, 0);
  std::vector<double> result;
  std::function<void(int, int)> visit = [&](int start, int remaining) {
    if (remaining == 0) {
      result.push_back(egpd_cont_partial(y, t, model, orders, 0, steps));
      return;
    }
    for (int i = start; i < n; ++i) {
      ++orders[i]; visit(i, remaining - 1); --orders[i];
    }
  };
  visit(0, higher ? 3 : 1);
  visit(0, higher ? 4 : 2);
  return arma::rowvec(result);
}

// --- discrete GPD CDF F = 1 - (1+xi m/sigma)^(-1/xi) (log-xi link) -------------
// F and its first/second derivatives w.r.t. (lsigma, lxi). ok=false if support is violated.
struct FDeriv { double F, F_ls, F_lx, F_lsls, F_lslx, F_lxlx; bool ok; };

static inline FDeriv gpd_Fderiv(double m, double sigma, double xi) {
  FDeriv r;
  double q = 1.0 / xi, v = xi * m / sigma, t = 1.0 + v;
  if (t <= 0.0) { r.ok = false; r.F = r.F_ls = r.F_lx = r.F_lsls = r.F_lslx = r.F_lxlx = 0.0; return r; }
  double L = std::log(t), w = std::pow(t, -q);
  r.F = 1.0 - w;
  double e_ls = q * v / t;
  double e_lx = q * (L - v / t);
  double e_lsls = -q * (v / t - v * v / (t * t));
  double e_lslx = -q * v * v / (t * t);
  double e_lxlx = -q * (L - v / t - v * v / (t * t));
  r.F_ls = -w * e_ls;
  r.F_lx = -w * e_lx;
  r.F_lsls = -w * (e_ls * e_ls + e_lsls);
  r.F_lslx = -w * (e_ls * e_lx + e_lslx);
  r.F_lxlx = -w * (e_lx * e_lx + e_lxlx);
  r.ok = true;
  return r;
}

// G and its theta-derivatives (theta = lsigma, lxi, lkappa) at one F point.
struct GTheta { double G, ls, lx, lk, lsls, lslx, lslk, lxlx, lxlk, lklk; };

// --- truncated-normal G: kappa-side precompute (shared across both F endpoints) --
struct TNkappa { double kappa, sk, skp, skpp, Fmin, phisk, D, Dk, Dkk, Fmin_kk; };

static inline TNkappa tn_kappa(double kappa) {
  TNkappa k;
  k.kappa = kappa;
  k.sk = std::sqrt(kappa);
  k.skp = 0.5 / k.sk;
  k.skpp = -k.skp / (2.0 * kappa);
  k.Fmin = R::pnorm(-k.sk, 0.0, 1.0, 1, 0);
  k.phisk = R::dnorm(k.sk, 0.0, 1.0, 0);
  k.D = 0.5 - k.Fmin; if (k.D < 1e-300) k.D = 1e-300;
  k.Dk = k.skp * k.phisk;
  k.Dkk = k.phisk * (k.skpp - k.skp * k.skp * k.sk);
  k.Fmin_kk = k.phisk * (k.skp * k.skp * k.sk - k.skpp);
  return k;
}

// truncated-normal G and full theta-derivatives at an F point (fully closed form).
static inline GTheta tn_Gtheta(const FDeriv& f, const TNkappa& k) {
  GTheta g;
  double sk = k.sk, D = k.D, kappa = k.kappa;
  double z = sk * (f.F - 1.0);
  double phiz = R::dnorm(z, 0.0, 1.0, 0);
  double Phiz = R::pnorm(z, 0.0, 1.0, 1, 0);
  double N = Phiz - k.Fmin;
  double G = N / D;
  double G_F = phiz * sk / D;
  double G_FF = -sk * sk * z * phiz / D;
  double zk = k.skp * (f.F - 1.0);
  double zkk = k.skpp * (f.F - 1.0);
  double Nk = phiz * zk + k.skp * k.phisk;
  double G_K = (Nk * D - N * k.Dk) / (D * D);
  double Nkk = phiz * (zkk - z * zk * zk) - k.Fmin_kk;
  double G_KK = (Nkk * D - N * k.Dkk) / (D * D) - 2.0 * k.Dk * G_K / D;
  double G_FK = ((-z * phiz) * zk * sk + phiz * k.skp) / D - phiz * sk * k.Dk / (D * D);
  g.G    = G;
  g.ls   = G_F * f.F_ls;
  g.lx   = G_F * f.F_lx;
  g.lk   = G_K * kappa;
  g.lsls = G_FF * f.F_ls * f.F_ls + G_F * f.F_lsls;
  g.lslx = G_FF * f.F_ls * f.F_lx + G_F * f.F_lslx;
  g.lxlx = G_FF * f.F_lx * f.F_lx + G_F * f.F_lxlx;
  g.lslk = G_FK * kappa * f.F_ls;
  g.lxlk = G_FK * kappa * f.F_lx;
  g.lklk = G_KK * kappa * kappa + G_K * kappa;
  return g;
}

// Gauss-Kronrod 7-15 nodes/weights, shared with the DEGPD Jet quadrature.
static const double EGPD_XGK[8] = {0.991455371120813,0.949107912342759,0.864864423359769,
  0.741531185599394,0.586087235467691,0.405845151377397,0.207784955007898,0.0};
static const double EGPD_WGK[8] = {0.022935322010529,0.063092092629979,0.104790010322250,
  0.140653259715525,0.169004726639267,0.190350578064785,0.204432940075298,0.209482141084728};
static const double EGPD_WG[4] = {0.129484966168870,0.279705391489277,0.381830050505119,0.417959183673469};

// --- truncated-beta carrier (type 3, omega = 1/32): kappa-side normaliser -------
// With s = 1 - F the GPD survival, t = 1/2 - (15/32) s and t(1 - t) =
// (1 - C s^2)/4, C = (15/16)^2, so the Beta(kappa, kappa) constant and 4^(1-kappa)
// cancel and the carrier density is
//   g(F) = (1 - C s^2)^(kappa - 1) / Z(kappa),  Z(kappa) = int_0^1 (1 - C s^2)^(kappa-1) ds.
// With L(s) = log(1 - C s^2) and moments under the normalised weight,
//   d log Z / d kappa = E[L],   d^2 log Z / d kappa^2 = Var[L].
// The weight and its L-moments are integrated directly (adaptive Gauss-Kronrod,
// as in the DEGPD model-6 Jet quadrature). The former route through
// d/dkappa of incomplete-beta functions (K1 - I*b1, ...) cancelled
// catastrophically and gave wrong kappa derivatives for kappa >~ 30.
struct TBnorm { double logZ, mean, var; };

static inline void egpd_tb_qk15(double a, double b, double km1, double lo, double width,
                                double v[3], double e[3]) {
  const double C = 225.0 / 256.0;
  double ctr = 0.5 * (a + b), hl = 0.5 * (b - a);
  double rk[3] = {0.0, 0.0, 0.0}, rg[3] = {0.0, 0.0, 0.0};
  for (int j = 0; j < 8; ++j) {
    int npt = (j == 7) ? 1 : 2;
    double f[3] = {0.0, 0.0, 0.0};
    for (int side = 0; side < npt; ++side) {
      double ds = width * (ctr + (side ? -1.0 : 1.0) * hl * EGPD_XGK[j]);
      // L(lo + ds) - L(lo), even for a cell deep in the tail.
      double L = std::log1p(-C * ds * (2.0 * lo + ds) / (1.0 - C * lo * lo));
      double w = std::exp(km1 * L);
      f[0] += w; f[1] += w * L; f[2] += w * L * L;
    }
    for (int q = 0; q < 3; ++q) {
      rk[q] += EGPD_WGK[j] * f[q];
      if (j == 7) rg[q] += EGPD_WG[3] * f[q];
      else if (j % 2 == 1) rg[q] += EGPD_WG[(j - 1) / 2] * f[q];
    }
  }
  for (int q = 0; q < 3; ++q) {
    v[q] = rk[q] * hl;
    e[q] = std::fabs((rk[q] - rg[q]) * hl);
  }
}

static inline void egpd_tb_adapt(double a, double b, double km1, double lo, double width,
                                 int depth, double I[3]) {
  const double tol = 1e-12;
  double v[3], e[3];
  egpd_tb_qk15(a, b, km1, lo, width, v, e);
  // A non-finite panel cannot converge: stop rather than bisect 2^depth times.
  bool finite = std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
  bool done = !finite || (e[0] <= tol * v[0] && e[1] <= tol * std::fabs(v[1]) &&
    e[2] <= tol * v[2]);
  if (done || depth >= 30) {
    for (int q = 0; q < 3; ++q) I[q] += v[q];
    return;
  }
  double m = 0.5 * (a + b);
  egpd_tb_adapt(a, m, km1, lo, width, depth + 1, I);
  egpd_tb_adapt(m, b, km1, lo, width, depth + 1, I);
}

static inline TBnorm tb_interval(double lkappa, double lo, double hi) {
  const double C = 225.0 / 256.0;
  double km1 = std::expm1(lkappa);
  // kappa = Inf (lkappa > 709.78) or NaN: no normaliser; callers see non-finite
  // derivatives and fall back to the out-of-support row (d0 is 1e20 there).
  if (!std::isfinite(km1)) {
    TBnorm bad = {R_NaN, R_NaN, R_NaN};
    return bad;
  }
  const double L0 = std::log1p(-C * lo * lo);
  double width = hi - lo;
  if (width == 0.0) return TBnorm{R_NegInf, L0, 0.0};
  if (km1 > 100.0) {
    // Drop only the region where the weight is below exp(-100) of its peak.
    double h2 = (1.0 - C * lo * lo) * -std::expm1(-100.0 / km1) / C;
    width = std::min(width, h2 / (std::sqrt(lo * lo + h2) + lo));
  }
  double I[3] = {0.0, 0.0, 0.0};
  egpd_tb_adapt(0.0, 1.0, km1, lo, width, 0, I);
  TBnorm r;
  r.logZ = std::log(width) + km1 * L0 + std::log(I[0]);
  double dL = I[1] / I[0];
  r.mean = L0 + dL;
  r.var = I[2] / I[0] - dL * dL;
  return r;
}

static inline TBnorm tb_normaliser(double lkappa) {
  return tb_interval(lkappa, 0.0, 1.0);
}

// log Z(kappa) for egpd_cont_nll. Rows sharing kappa (e.g. an intercept-only
// kappa) ask for the same value, so the last one is kept; it depends on lkappa only.
static inline double egpd_tb_logZ(double lkappa) {
  static double last_lk = R_NaN, last_logZ = R_NaN;
  if (!(lkappa == last_lk)) {
    last_logZ = tb_normaliser(lkappa).logZ;
    last_lk = lkappa;
  }
  return last_logZ;
}

// --- exact derivatives of the continuous objective (forward-mode AD) ----------
// The same stabilised formulas as egpd_cont_nll (y > 0, in support), evaluated
// with second-order Jets. Used for d12 wherever the generated closed forms are
// unreliable (y/sigma -> 0, xi -> 0) or overflow (e.g. next to the xi < 0
// endpoint). This replaces predictor-scale finite-difference stencils, which
// carried ~eps*|NLL|/h^2 noise and could step across the endpoint into the
// 1e20 out-of-support penalty.
template<int N>
inline egpd_ad::Jet<N> egpd_softplus_jet(const egpd_ad::Jet<N>& x) {
  using namespace egpd_ad;
  return x.v > 0.0 ? x + jlog1p(jexp(-x)) : jlog1p(jexp(x));
}

// Jet version of egpd_log_G3 (same branches).
template<int N>
inline egpd_ad::Jet<N> egpd_log_G3_jet(const egpd_ad::Jet<N>& A, const egpd_ad::Jet<N>& logA,
                                       const egpd_ad::Jet<N>& ldelta) {
  using namespace egpd_ad;
  typedef Jet<N> J;
  if (A.v >= 1.0) {
    J ratio;
    if (-std::exp(std::min(ldelta.v, 700.0)) * A.v < -50.0) {
      ratio = -jexp(-ldelta);
    } else {
      J d = jexp(ldelta), z = -(d * A);
      ratio = std::abs(z.v) < 0.1 ? -(A * expm1_over_z(z)) : jexpm1(z) / d;
    }
    return log1mexp(-A + jlog1p(-ratio));
  }
  J ldu = ldelta + logA;
  if (ldelta.v <= 40.0) {
    J dp1 = J(1.0) + jexp(ldelta);
    if (dp1.v * A.v < 0.2) {
      J P[17];
      P[0] = J(1.0);
      for (int k = 1; k < 17; ++k) P[k] = J(1.0) + dp1 * P[k - 1];
      J s = P[16];
      for (int k = 15; k >= 0; --k)
        s = P[k] - (1.0 / static_cast<double>(k + 3)) * (A * s);
      return jlog(dp1) - std::log(2.0) + 2.0 * logA + jlog(s);
    }
  } else if (ldu.v < -1.6) {
    J du = jexp(ldu);
    J s(1.0);
    for (int k = 16; k >= 1; --k)
      s = J(1.0) - (1.0 / static_cast<double>(k + 2)) * (du * s);
    return logA + ldu - std::log(2.0) + jlog(s);
  }
  J term2 = (ldu.v > 50.0) ? jexp(-ldu) : expm1_over_z(-jexp(ldu));
  return -A + logA + jlog(expm1_over_z(A) - term2);
}

template<int N>
inline egpd_ad::Jet<N> egpd_cont_nll_jet(double y, const egpd_ad::Jet<N>* p, int model) {
  using namespace egpd_ad;
  typedef Jet<N> J;
  J z = y * jexp(-p[0]);
  J v = p[1] * z;
  J r;                                   // log1p(v)/v, analytic through xi = 0
  if (std::abs(v.v) < 0.1) {
    r = J(1.0 / 21.0);
    for (int k = 19; k >= 0; --k)
      r = J((k % 2 ? -1.0 : 1.0) / (k + 1.0)) + v * r;
  } else {
    r = jlog1p(v) / v;
  }
  J A = z * r;                           // -log(1 - H(y))
  J logA = J(std::log(y)) - p[0] + jlog(r);
  J logH = (A.v < 0.1) ? (logA + jlog(expm1_over_z(-A))) : log1mexp(-A);
  J logg;
  if (model == 1) {
    logg = p[2] + (jexp(p[2]) - 1.0) * logH;
  } else if (model == 2 || model == -2) {
    J logd = (model == 2) ? logadd(p[2], p[3]) : p[3];
    J lp = -egpd_softplus_jet(-p[4]), lq = -egpd_softplus_jet(p[4]);
    logg = logadd(lp + p[2] + (jexp(p[2]) - 1.0) * logH,
                  lq + logd + (jexp(logd) - 1.0) * logH);
  } else if (model == 3 || model == 4) {
    J ldu = p[2] + logA;
    if (ldu.v < -2.3) {
      logg = egpd_softplus_jet(p[2]) + logA + jlog(expm1_over_z(-jexp(ldu)));
    } else {
      J l1m = (ldu.v > 50.0) ? J(0.0) : log1mexp(-jexp(ldu));
      logg = egpd_softplus_jet(-p[2]) + l1m;
    }
    if (model == 4)
      logg = logg + (p[3] - std::log(2.0)) +
        (0.5 * jexp(p[3]) - 1.0) * egpd_log_G3_jet(A, logA, p[2]);
  } else if (model == 5) {
    double k = std::exp(p[2].v);
    if (egpd_tn_subnormal(k)) {
      // -log M(kappa) = kappa/6 + O(kappa^2), each lkappa derivative also kappa/6
      logg = unary(p[2], k / 6.0, k / 6.0, k / 6.0) - 0.5 * jexp(p[2]) * jexp(-2.0 * A);
    } else {
      double sk = std::sqrt(k);
      double denom = egpd_tn_norm(k);
      double ratio = sk * R::dnorm(sk, 0.0, 1.0, 0) / (2.0 * denom);
      J lognorm = unary(p[2], std::log(denom), ratio, 0.5 * ratio * (1.0 - k - 2.0 * ratio));
      logg = 0.5 * p[2] - 0.5 * std::log(2.0 * M_PI) -
        0.5 * jexp(p[2]) * jexp(-2.0 * A) - lognorm;
    }
  } else {
    TBnorm nz = tb_normaliser(p[2].v);
    double k = std::exp(p[2].v);
    J logZ = unary(p[2], nz.logZ, k * nz.mean, k * nz.mean + k * k * nz.var);
    J s = jexp(-A);
    logg = (jexp(p[2]) - 1.0) * jlog1p(-(225.0 / 256.0) * (s * s)) - logZ;
  }
  return p[0] + (1.0 + p[1]) * A - logg;
}

// Gradient and packed upper-triangular Hessian (d12 layout) at one row.
static inline arma::rowvec egpd_cont_jet_row(double y, const std::vector<double>& t, int model) {
  using namespace egpd_ad;
  const int n = t.size();
  Jet<5> p[5];
  for (int i = 0; i < n; ++i) p[i] = variable<5>(t[i], i);
  Jet<5> f = egpd_cont_nll_jet<5>(y, p, model);
  arma::rowvec out(n + n * (n + 1) / 2);
  for (int i = 0; i < n; ++i) out[i] = f.g[i];
  int c = n;
  for (int i = 0; i < n; ++i)
    for (int j = i; j < n; ++j) out[c++] = f.h[i][j];
  return out;
}

// Guard out-of-support, boundary y == 0, near-zero xi, and tiny y/sigma
// consistently with egpd_cont_nll across all continuous EGPD d12 and d34 routines.
static inline bool egpd_cont_try_special(double y, const std::vector<double>& t,
                                         int model, bool higher,
                                         arma::mat& out, int j) {
  for (size_t i = 0; i < t.size(); ++i) {
    if (!std::isfinite(t[i])) {
      out.row(j).zeros();
      return true;
    }
  }
  if (egpd_cont_oos(y, t[0], t[1])) {
    out.row(j).zeros();
    return true;
  }
  if (y == 0.0) {
    out.row(j).zeros();
    if (model == 1) {
      if (t[2] == 0.0 && !higher) out(j, 0) = 1.0;
      return true;
    }
    if (model == 2) {
      if (t[2] == 0.0 && !higher) {
        double p = 1.0 / (1.0 + std::exp(-t[4]));
        double q = 1.0 / (1.0 + std::exp(t[4]));
        out(j, 0) = 1.0;
        out(j, 4) = -q;
        out(j, 19) = p * q;
      }
      return true;
    }
    if (model == -2) {
      if ((t[2] == 0.0 && t[3] > 0.0) || (t[3] == 0.0 && t[2] > 0.0)) {
        double p = 1.0 / (1.0 + std::exp(-t[4]));
        double q = 1.0 / (1.0 + std::exp(t[4]));
        double pq = p * q;
        out(j, 34) = pq * (1.0 - 2.0 * p);
        out(j, 104) = pq * (1.0 - 6.0 * pq);
      }
      return true;
    }
    if (model == 3) {
      return true;
    }
    if (model == 4) {
      if (t[3] == 0.0) {
        double q = 1.0 / (1.0 + std::exp(-t[2]));
        double q1mq = q * (1.0 / (1.0 + std::exp(t[2])));
        if (!higher) {
          out(j, 0) = 1.0;
          out(j, 2) = -0.5 * q;
          out(j, 11) = -0.5 * q1mq;
        } else {
          out(j, 16) = -0.5 * q1mq * (1.0 - 2.0 * q);
          out(j, 50) = -0.5 * q1mq * (1.0 - 6.0 * q1mq);
        }
      }
      return true;
    }
    double kappa = std::exp(t[2]);
    if (!std::isfinite(kappa)) return true;
    if (model == 5) {
      if (egpd_tn_subnormal(kappa)) {
        // NLL = lpsi + kappa/2 + log M(kappa): lkappa derivatives kappa/2 - kappa/6
        out(j, 0) = 1.0;
        out(j, 2) = kappa / 3.0;
        out(j, 8) = kappa / 3.0;
        return true;
      }
      double sqk = std::sqrt(kappa);
      double denom = egpd_tn_norm(kappa);
      double phisqk = R::dnorm(sqk, 0.0, 1.0, 0);
      double ratio = sqk * phisqk / (2.0 * denom);
      out(j, 0) = 1.0;
      out(j, 2) = -0.5 + 0.5 * kappa + ratio;
      out(j, 8) = 0.5 * kappa + 0.5 * ratio * (1.0 - kappa - 2.0 * ratio);
      return true;
    }
    if (model == 6) {
      // s = 1 at y = 0: NLL = lpsi - (kappa - 1) log(1 - C) + log Z(kappa).
      TBnorm nz = tb_normaliser(t[2]);
      double dk = kappa * (nz.mean - std::log1p(-225.0 / 256.0));
      out(j, 0) = 1.0;
      out(j, 2) = dk;
      out(j, 8) = dk + kappa * kappa * nz.var;
      return true;
    }
  }
  double z = y * std::exp(-t[0]);
  bool small_y = z < (higher ? 1e-2 : 1e-4);
  // Far in the upper tail, or next to the xi < 0 endpoint, the closed-form d12
  // forms log H from H = 1 - s and loses the O(s) terms (e.g. the log-kappa
  // curvature -kappa log H); A = -log s > 16 means s < 1.2e-7.
  double A = egpd_tail_exponent(z, t[1]);
  bool far_tail = !higher && A > 16.0;
  // Carriers 3/4: the closed-form d12 forms 1 - exp(-delta A) (carrier 4 also
  // G3) by subtraction and loses ~eps/(delta A) (carrier 4: ~eps/(delta
  // min(A, 1)^2)) as these vanish; below delta ~ 1e-9 the log-delta score has
  // the wrong sign. The Jet path below uses series there.
  bool small_delta = false;
  if (!higher && (model == 3 || model == 4)) {
    double a1 = (model == 3) ? A : std::min(A, 1.0) * std::min(A, 1.0);
    small_delta = std::exp(t[2]) * a1 < 1e-4;
  }
  // Carrier 5 with subnormal kappa: the closed form's normaliser ratio
  // sqrt(kappa) phi(sqrt(kappa)) / (2 (Phi(sqrt(kappa)) - 1/2)) is formed from
  // differently rounded subnormals (see egpd_tn_subnormal); the Jet path is exact.
  bool small_kappa5 = model == 5 && egpd_tn_subnormal(std::exp(t[2]));
  if (egpd_cont_near_zero(y, t) || small_y || far_tail || small_delta || small_kappa5) {
    double f0 = egpd_cont_nll(y, t, model);
    if (!std::isfinite(f0) || f0 >= 1e19) {
      out.row(j).zeros();
      return true;
    }
    // d12: exact Jet derivatives; d34: predictor-scale stencils.
    out.row(j) = higher ? egpd_cont_derivs(y, t, model, true)
                        : egpd_cont_jet_row(y, t, model);
    if (!out.row(j).is_finite()) out.row(j).zeros();
    return true;
  }
  return false;
}

static inline void egpd_cont_finish_row(double y, const std::vector<double>& t,
                                        int model, bool higher,
                                        arma::mat& out, int j) {
  if (!out.row(j).is_finite()) {
    double f0 = egpd_cont_nll(y, t, model);
    if (std::isfinite(f0) && f0 < 1e19) {
      out.row(j) = higher ? egpd_cont_derivs(y, t, model, true)
                          : egpd_cont_jet_row(y, t, model);
      if (!out.row(j).is_finite()) out.row(j).zeros();
    } else {
      out.row(j).zeros();
    }
  }
}

#endif
