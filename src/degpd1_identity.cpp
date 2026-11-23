// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "gfunc_derivs.h"
#include "egpd_jet.h"
#include "gpd_endpoint.h"

namespace {
using namespace egpd_ad;
// Second-order forward differentiation in the predictors (never in log(xi)).
// Keeping the removable singularity inside log1p(z)/z also differentiates the
// exponential limit correctly: replacing a small xi by zero would not.
// Cells are templated on the number of active derivative directions: NPAR for
// d12 and 0 (value only) for d0, so objective-only calls skip the Hessians.
const int NPAR = 5;
// log(1 - exp(x)) for x <= 0 (log1mexp(-Inf) = 0, log1mexp(0) = -Inf); NaN for
// x > 0, so callers must keep the argument non-positive.
inline double log1mexp_d(double x) {
  return x > -0.6931471805599453 ? std::log(-std::expm1(x)) : std::log1p(-std::exp(x));
}
template<int N>
struct GpdPoint {
  Jet<N> e;   // log(1 - H(t)) in [-Inf, 0]
  Jet<N> lh;  // log(H(t))     in [-Inf, 0]
  Jet<N> lu;  // log(-e)       when -Inf < e < 0
};
template<int N>
GpdPoint<N> gpd_point(double t, const Jet<N>& ls, const Jet<N>& xi) {
  if (t <= 0.0) return {Jet<N>(0.0), Jet<N>(R_NegInf), Jet<N>(R_NegInf)};
  Jet<N> a = t * jexp(-ls), z = xi * a;
  if (z.v < -0.5) {
    // Keep the chain derivatives of 1+xi*t/scale, but evaluate its value
    // from the response-scale parameters, as in the standalone count law.
    Jet<N> gap = 1.0 + z;
    const double scale = std::exp(ls.v);
    if (scale > 0 && std::isfinite(scale)) gap.v = egpd_endpoint_gap(t, scale, xi.v);
    if (gap.v <= 0) return {Jet<N>(R_NegInf), Jet<N>(0.0), Jet<N>(R_PosInf)};
    Jet<N> e = -jlog(gap)/xi;
    return {e, log1mexp(e), jlog(-e)};
  }
  Jet<N> r;
  if (std::abs(z.v) < 0.1) {
    // At |z| <= .1 the omitted terms in r'' are < 4e-17.
    r = Jet<N>(1.0/21.0);
    for (int k = 19; k >= 0; --k)
      r = Jet<N>((k % 2 ? -1.0 : 1.0)/(k+1.0)) + z*r;
  } else {
    r = jlog1p(z)/z;
  }
  Jet<N> u = a * r;
  Jet<N> e = -u;
  Jet<N> lu = Jet<N>(std::log(t)) - ls + jlog(r);
  Jet<N> lh = (u.v < 0.1) ? (lu + jlog(expm1_over_z(-u))) : log1mexp(e);
  return {e, lh, lu};
}
template<int N>
Jet<N> log_power_sf_from_pair(const Jet<N>& e, const Jet<N>& lh, const Jet<N>& k,
                              const Jet<N>& lk) {
  if (e.v == R_NegInf) return Jet<N>(R_NegInf);
  if (lh.v == R_NegInf) return Jet<N>(0.0);
  // A subnormal k (log k below -708.4) keeps only a few significant bits while
  // lk = log(k) is exact, and a subnormal k*lh has lost its relative precision;
  // log(k) is then taken from lk (the cell objective varied by O(1) there).
  bool sub = k.v < std::numeric_limits<double>::min();
  // 1-(1-s)^k = k*s*(1+O(k*s)). Below machine precision, retaining
  // log(k)+log(s) avoids underflow; this is not a count/shape cutoff.
  if (e.v < -36.0 && e.v + std::log(std::max(1.0, k.v)) < -36.0)
    return (sub ? lk : jlog(k)) + e;
  Jet<N> klh = k * lh;
  // 1-(1-s)^k = -k*lh*(1+O(k*lh)), |k*lh| < 4e-305 here (|lh| < 1500)
  if (sub || std::abs(klh.v) < std::numeric_limits<double>::min())
    return lk + jlog(-lh);
  return log1mexp(klh);
}
template<int N>
Jet<N> carrier3_logsf(const GpdPoint<N>& pt, const Jet<N>& ldelta) {
  if (pt.e.v == R_NegInf) return Jet<N>(R_NegInf);
  if (pt.lh.v == R_NegInf) return Jet<N>(0.0);
  double d_clamp = std::exp(std::min(ldelta.v, 700.0));
  if (d_clamp * pt.e.v < -50.0) {
    return pt.e + jlog1p(jexp(-ldelta));
  }
  Jet<N> delta = jexp(ldelta);
  Jet<N> z = delta * pt.e;
  // 1-G3(1-s) = s * (1 - (s^delta - 1)/delta); factor e out near delta*e = 0.
  Jet<N> ratio = std::abs(z.v) < 0.1 ? pt.e * expm1_over_z(z) : jexpm1(z) / delta;
  return pt.e + jlog1p(-ratio);
}
template<int N>
Jet<N> carrier3_logcdf(const GpdPoint<N>& pt, const Jet<N>& ldelta) {
  if (pt.lh.v == R_NegInf) return Jet<N>(R_NegInf);
  if (pt.e.v == R_NegInf) return Jet<N>(0.0);
  Jet<N> u = -pt.e;
  if (u.v >= 1.0) return log1mexp(carrier3_logsf(pt, ldelta));
  Jet<N> ldu = ldelta + pt.lu;
  if (ldelta.v <= 40.0) {
    Jet<N> dp1 = Jet<N>(1.0) + jexp(ldelta);
    if (dp1.v * u.v < 0.2) {
      // G3(1-exp(-u)) = ((delta+1)/2)*u^2 * sum_{k>=0} (-1)^k P_{k+1}(delta) 2/(k+2)! u^k
      Jet<N> P[17];
      P[0] = Jet<N>(1.0);
      for (int k = 1; k < 17; ++k) P[k] = Jet<N>(1.0) + dp1 * P[k - 1];
      Jet<N> s = P[16];
      for (int k = 15; k >= 0; --k)
        s = P[k] - (1.0 / static_cast<double>(k + 3)) * (u * s);
      return jlog(dp1) - std::log(2.0) + 2.0 * pt.lu + jlog(s);
    }
  } else if (ldu.v < -1.6) {
    Jet<N> du = jexp(ldu);
    Jet<N> s(1.0);
    for (int k = 16; k >= 1; --k)
      s = Jet<N>(1.0) - (1.0 / static_cast<double>(k + 2)) * (du * s);
    return pt.lu + ldu - std::log(2.0) + jlog(s);
  }
  Jet<N> term2 = (ldu.v > 50.0) ? jexp(-ldu) : expm1_over_z(-jexp(ldu));
  return pt.e + pt.lu + jlog(expm1_over_z(u) - term2);
}
template<int N>
Jet<N> carrier_logsf(int m, const GpdPoint<N>& pt, const Jet<N>* p) {
  if (pt.e.v == R_NegInf) return Jet<N>(R_NegInf);
  if (pt.lh.v == R_NegInf) return Jet<N>(0.0);
  if (m == 1) return log_power_sf_from_pair(pt.e, pt.lh, jexp(p[2]), p[2]);
  if (m == 2) {
    Jet<N> k1 = jexp(p[2]), k2 = k1 + jexp(p[3]);
    Jet<N> lp, lq;
    if (p[4].v >= 0.0) {
      lp = -jlog1p(jexp(-p[4])); lq = lp-p[4];
    } else {
      lq = -jlog1p(jexp(p[4])); lp = lq+p[4];
    }
    return logadd(lp + log_power_sf_from_pair(pt.e, pt.lh, k1, p[2]),
                  lq + log_power_sf_from_pair(pt.e, pt.lh, k2, logadd(p[2], p[3])));
  }
  Jet<N> b = carrier3_logsf(pt, p[2]);
  if (m == 3) return b;
  Jet<N> lcdf3 = carrier3_logcdf(pt, p[2]);
  return log_power_sf_from_pair(b, lcdf3, 0.5 * jexp(p[3]), p[3] - std::log(2.0));
}
template<int N>
Jet<N> carrier_logcdf(int m, const GpdPoint<N>& pt, const Jet<N>* p) {
  if (pt.lh.v == R_NegInf) return Jet<N>(R_NegInf);
  if (pt.e.v == R_NegInf) return Jet<N>(0.0);
  if (m == 1) return jexp(p[2]) * pt.lh;
  if (m == 2) {
    Jet<N> k1 = jexp(p[2]), k2 = k1 + jexp(p[3]);
    Jet<N> lp, lq;
    if (p[4].v >= 0.0) {
      lp = -jlog1p(jexp(-p[4])); lq = lp-p[4];
    } else {
      lq = -jlog1p(jexp(p[4])); lp = lq+p[4];
    }
    return logadd(lp + k1 * pt.lh, lq + k2 * pt.lh);
  }
  Jet<N> f3 = carrier3_logcdf(pt, p[2]);
  if (m == 3) return f3;
  return 0.5 * jexp(p[3]) * f3;
}
template<int N, class F>
Jet<N> log_carrier_integral(const F& f, double a = 0, double b = 1, int depth = 0);

template<int N>
struct BetaDensity {
  Jet<N> k;
  Jet<N> operator()(double s) const {
    return (k-1.0)*std::log1p(-4.0*(15.0/32.0)*(15.0/32.0)*s*s);
  }
};

template<int N>
struct CarrierDensity {
  int model;
  Jet<N> k, head, s0, w, inv_omc2s0;
  CarrierDensity(int m, const Jet<N>& lo, const Jet<N>& hi, const Jet<N>& lk) :
    model(m), k(jexp(lk)) {
    Jet<N> slo = jexp(lo);
    s0 = hi.v == R_NegInf ? Jet<N>(0.0) : jexp(hi);
    w = hi.v == R_NegInf ? slo : slo * (-jexpm1(hi - lo));
    Jet<N> lognorm;
    if (m == 5) {
      Jet<N> a = jexp(0.5 * lk)/std::sqrt(2.0);
      if (a.v < 0.05) {
        Jet<N> a2 = a * a;
        const double denom_c[8] = {1.0, 3.0, 10.0, 42.0, 216.0, 1320.0, 9360.0, 75600.0};
        Jet<N> s(1.0 / 685440.0);
        for (int n = 7; n >= 0; --n)
          s = Jet<N>(1.0 / denom_c[n]) - a2 * s;
        lognorm = jlog(s);
      } else {
        double d = 2.0/std::sqrt(arma::datum::pi)*std::exp(-a.v*a.v);
        Jet<N> erf = unary(a, std::erf(a.v), d, -2*a.v*d);
        lognorm = jlog(0.5 * erf) + 0.5*std::log(2*arma::datum::pi) - 0.5 * lk;
      }
      head = -0.5 * k * (s0 * s0) - lognorm;
    } else {
      // The beta-function constant and interval Jacobian cancel. Integrating
      // the unnormalized density avoids subtracting incomplete-beta moments.
      // Truncate negligible Gaussian tail at large k so adaptive GK finishes in O(1) splits.
      double b_upper = (R_finite(k.v) && k.v > 120.0)
        ? std::min(1.0, 10.0 / ((15.0 / 16.0) * std::sqrt(k.v - 1.0)))
        : 1.0;
      lognorm = log_carrier_integral<N>(BetaDensity<N>{k}, 0.0, b_upper);
      Jet<N> c2s0 = (225.0 / 256.0) * (s0 * s0);
      inv_omc2s0 = reciprocal(Jet<N>(1.0) - c2s0);
      head = (k - 1.0) * jlog1p(-c2s0) - lognorm;
    }
  }
  Jet<N> operator()(double t) const {
    Jet<N> ds2 = w * (2.0 * s0 + t * w);
    if (model == 5) return (-0.5 * t) * (k * ds2);
    return (k - 1.0) * jlog1p((-(225.0 / 256.0) * t) * (ds2 * inv_omc2s0));
  }
};
template<int N>
bool quad_close(const Jet<N>& a, const Jet<N>& b) {
  const double tol = 2e-11;
  if (!(a.v > 0) || !R_finite(a.v) || !R_finite(b.v) ||
      std::abs(a.v-b.v) > tol*a.v) return false;
  for (int i = 0; i < N; ++i) {
    if (!R_finite(a.g[i]) || !R_finite(b.g[i]) ||
        std::abs(a.g[i]-b.g[i]) > tol*(a.v+std::abs(a.g[i]))) return false;
    for (int j = 0; j < N; ++j)
      if (!R_finite(a.h[i][j]) || !R_finite(b.h[i][j]) ||
          std::abs(a.h[i][j]-b.h[i][j]) > tol*(a.v+std::abs(a.h[i][j]))) return false;
  }
  return true;
}
// Integrate G'(1-s) over the actual survival interval, factoring its width out
// in log space. Unlike subtracting pbeta/pnorm CDFs this retains arbitrarily
// small cell probabilities. The embedded 7/15 rule checks derivatives as well
// (value only for Jet<0>, so d0 may stop at a coarser, equally accurate split).
template<int N, class F>
Jet<N> log_carrier_integral(const F& f, double a, double b, int depth) {
  double mid = (a+b)/2, half = (b-a)/2;
  Jet<N> v[15];
  v[14] = f(mid);
  int imax = 14;
  for (int i = 0; i < 7; ++i) {
    v[2*i] = f(mid-half*EGPD_XGK[i]);
    v[2*i+1] = f(mid+half*EGPD_XGK[i]);
    if (v[2*i].v > v[imax].v) imax = 2*i;
    if (v[2*i+1].v > v[imax].v) imax = 2*i+1;
  }
  const Jet<N>& base = v[imax];
  if (!R_finite(base.v)) return Jet<N>(R_NaN);
  Jet<N> c14 = jexp(v[14]-base);
  Jet<N> high = EGPD_WGK[7]*c14, low = EGPD_WG[3]*c14;
  for (int i = 0; i < 7; ++i) {
    Jet<N> pair = jexp(v[2*i]-base) + jexp(v[2*i+1]-base);
    high = high + EGPD_WGK[i]*pair;
    if (i % 2) low = low + EGPD_WG[(i-1)/2]*pair;
  }
  if (quad_close(high, low) || depth >= 20) return jlog(half*high)+base;
  return logadd(log_carrier_integral<N>(f,a,mid,depth+1),
                log_carrier_integral<N>(f,mid,b,depth+1));
}
// Very large counts: S(y) and S(y+1) then agree to ~1/(sigma + xi*y), and
// differencing them (below) loses ~eps*y*log(y) of the cell probability (the
// NLL was 5e-3 off at y = 1e12, and 1e20 once y + 1 == y). Across such a cell
// the density is flat, so integrate the continuous EGPD density (same carriers
// and predictors, egpd_cont_nll_jet) by 3-point Gauss-Legendre, accepted only
// where it agrees with the 2-point rule. Otherwise (e.g. next to the xi < 0
// endpoint, or a light tail with small sigma, where the difference is accurate)
// the cell is differenced as before.
const double BIG_COUNT = 65536.0;
// The 2- and 3-point log-integrals must agree in value and, when derivatives are
// carried (N > 0), in gradient and Hessian. Next to the xi < 0 endpoint at
// xi = -1/k the density is a polynomial, so the values agree exactly while the
// xi-derivatives carry log(1 + xi t/sigma) singularities that both rules miss
// (Hessians were 25-79% off within ~1e-8 of xi = -1). Away from the endpoint the
// derivatives agree to ~1e-12 (relative), so those cells are unaffected.
template<int N>
bool big_count_agree(const Jet<N>& a, const Jet<N>& b) {
  if (!(std::abs(a.v - b.v) <= 1e-10)) return false;
  const double tol = 1e-9;
  for (int i = 0; i < N; ++i) {
    if (!(std::abs(a.g[i] - b.g[i]) <= tol * (1.0 + std::abs(a.g[i])))) return false;
    for (int j = i; j < N; ++j)
      if (!(std::abs(a.h[i][j] - b.h[i][j]) <= tol * (1.0 + std::abs(a.h[i][j])))) return false;
  }
  return true;
}
template<int N>
bool big_count_cell(double y, const Jet<N>* p, int m, Jet<N>& out) {
  // the top attainable bin (y + 1 beyond the xi < 0 endpoint) is left to the difference
  if (p[1].v < 0) {
    const double scale = std::exp(p[0].v);
    if (!(scale > 0 && std::isfinite(scale)) ||
        egpd_endpoint_gap(y + 1.0, scale, p[1].v) <= 0) return false;
  }
  const double x2[2] = {0.5 - 0.28867513459481287, 0.5 + 0.28867513459481287};
  const double x3[3] = {0.5 - 0.38729833462074169, 0.5, 0.5 + 0.38729833462074169};
  const double w3[3] = {5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0};
  Jet<N> l2(R_NegInf), l3(R_NegInf);
  for (int i = 0; i < 2; ++i) {
    Jet<N> f = std::log(0.5) - egpd_cont_nll_jet<N>(y + x2[i], p, m);
    if (!R_finite(f.v)) return false;
    l2 = logadd(l2, f);
  }
  for (int i = 0; i < 3; ++i) {
    Jet<N> f = std::log(w3[i]) - egpd_cont_nll_jet<N>(y + x3[i], p, m);
    if (!R_finite(f.v)) return false;
    l3 = logadd(l3, f);
  }
  if (!big_count_agree(l3, l2)) return false;
  out = -l3;
  return true;
}
template<int N>
Jet<N> cell(double y, const arma::rowvec& eta, int m) {
  if (y < 0.0 || !std::isfinite(y) || !eta.is_finite()) return Jet<N>(1e20);
  Jet<N> p[NPAR];
  for (unsigned i = 0; i < eta.n_elem; ++i) p[i] = variable<N>(eta[i], i);
  GpdPoint<N> lo = gpd_point(y, p[0], p[1]);
  if (lo.e.v == R_NegInf) return Jet<N>(1e20); // constant out-of-support penalty
  if (y >= BIG_COUNT) {
    Jet<N> big;
    if (big_count_cell<N>(y, p, m, big)) return big;
  }
  GpdPoint<N> hi = gpd_point(y+1.0, p[0], p[1]);
  if (m <= 4) {
    Jet<N> sl = carrier_logsf(m, lo, p);
    // Upper tail (F(y) > 1/2): difference survival probabilities.
    if (sl.v < -0.6931471805599453) {
      Jet<N> sh = carrier_logsf(m, hi, p);
      return -logdiff(sl, sh);
    }
    // Lower tail / bulk (F(y) <= 1/2): evaluate log-CDF directly so small F(y+1)
    // never underflows exp(logF) in log(1 - F(y+1)).
    Jet<N> fl = carrier_logcdf(m, lo, p), fh = carrier_logcdf(m, hi, p);
    if (fl.v == R_NegInf) return fh.v == R_NegInf ? Jet<N>(1e20) : -fh;
    return -logdiff(fh, fl);
  }
  if (!R_finite(std::exp(p[2].v))) return Jet<N>(1e20);
  if (lo.e.v == 0.0 && hi.e.v == R_NegInf) return Jet<N>(0.0);
  if (!(hi.e.v < lo.e.v)) return Jet<N>(1e20);
  CarrierDensity<N> density(m, lo.e, hi.e, p[2]);
  if (!R_finite(density.head.v)) return Jet<N>(1e20);
  double s0 = density.s0.v;
  double w = density.w.v;
  double rate = (m == 5)
    ? (0.5 * density.k.v)
    : ((density.k.v - 1.0) * (15.0 / 16.0) * (15.0 / 16.0) * density.inv_omc2s0.v);
  double t_upper = 1.0;
  if (rate > 100.0 && w > 0.0) {
    t_upper = std::min(t_upper, 10.0 / (w * std::sqrt(rate)));
    if (s0 > 0.0) t_upper = std::min(t_upper, 50.0 / (rate * s0 * w));
  }
  if (!(t_upper > 0.0)) return Jet<N>(1e20);
  Jet<N> log_width = hi.e.v == R_NegInf ? Jet<N>(0.0) : log1mexp(hi.e - lo.e);
  return -(lo.e + log_width + density.head + log_carrier_integral<N>(density, 0.0, t_upper));
}
arma::mat predictors(const Rcpp::List& pars, const std::vector<arma::mat>& X,
                     const arma::uvec& dupid, int dcate, const Rcpp::List& offsets) {
  arma::mat eta;
  for (unsigned i = 0; i < X.size(); ++i) {
    arma::vec v = X[i] * Rcpp::as<arma::vec>(pars[i]);
    if (dcate == 1) v = v.elem(dupid);
    arma::vec off = Rcpp::as<arma::vec>(offsets[i]);
    if (off.n_elem) v += off;
    if (!i) eta.set_size(v.n_elem, X.size());
    eta.col(i) = v;
  }
  return eta;
}
}

// Shared, non-exported entry points keep the Rcpp interface unchanged.
// With value_only = true only the NLL column is computed (value-only Jets), so
// d0 does not pay for gradients and Hessians that it would discard.
arma::mat degpd_identity_eval(const Rcpp::List& pars, const std::vector<arma::mat>& X,
                             const arma::vec& y, const arma::uvec& dupid, int dcate,
                             const Rcpp::List& offsets, int m, bool value_only) {
  arma::mat eta = predictors(pars, X, dupid, dcate, offsets);
  int n = X.size();
  arma::mat out(y.n_elem, value_only ? 1 : 1+n+n*(n+1)/2, arma::fill::zeros);
  for (unsigned r = 0; r < y.n_elem; ++r) {
    if (value_only) {
      double v = cell<0>(y[r], eta.row(r), m).v;
      out(r, 0) = (v >= 1e19 || !R_finite(v)) ? 1e20 : v;
      continue;
    }
    Jet<NPAR> z = cell<NPAR>(y[r], eta.row(r), m);
    if (z.v >= 1e19 || !R_finite(z.v)) {
      out(r, 0) = 1e20;
      continue;
    }
    out(r, 0) = z.v;
    for (int i = 0; i < n; ++i) out(r, i+1) = z.g[i];
    int c = n+1;
    for (int i = 0; i < n; ++i) for (int j = i; j < n; ++j) out(r, c++) = z.h[i][j];
  }
  return out;
}

// Log-link (and bounded-shape) DEGPD evaluator for models 1..6.
// Maps lxi -> xi = exp(lxi_eff) into the stabilized identity Jet evaluator and
// applies the exact log-link chain rule, preserving full precision as xi -> 0+.
arma::mat degpd_loglink_eval(const Rcpp::List& pars, const std::vector<arma::mat>& X,
                             const arma::vec& y, const arma::uvec& dupid, int dcate,
                             const Rcpp::List& offsets, int m, double xi_max,
                             bool apply_bounded_chain, bool value_only) {
  arma::mat eta = predictors(pars, X, dupid, dcate, offsets);
  int n = X.size();
  arma::mat out(y.n_elem, value_only ? 1 : 1 + n + n * (n + 1) / 2, arma::fill::zeros);
  for (unsigned r = 0; r < y.n_elem; ++r) {
    arma::rowvec row = eta.row(r);
    if (!row.is_finite()) {
      out(r, 0) = 1e20;
      continue;
    }
    double lxi_eff = row[1], mp_xi = 1.0, mpp_xi = 0.0;
    if (R_finite(xi_max)) {
      double log_xm = std::log(xi_max);
      double t = std::exp(row[1] - log_xm);
      if (t > 700.0) {
        lxi_eff = log_xm;
        mp_xi = 0.0;
        mpp_xi = 0.0;
      } else if (t < 1e-4) {
        lxi_eff = row[1] + std::log1p(t * (-0.5 + t * (1.0 / 6.0 - t / 24.0)));
        mp_xi = 1.0 - 0.5 * t + t * t / 12.0;
        mpp_xi = mp_xi * (1.0 - t - mp_xi);
      } else {
        double em = std::exp(-t);
        double omem = -std::expm1(-t);
        lxi_eff = log_xm + std::log(omem);
        mp_xi = (t * em) / omem;
        mpp_xi = mp_xi * (1.0 - t - mp_xi);
      }
    }
    double xi = std::exp(lxi_eff);
    if (!R_finite(xi)) {
      out(r, 0) = 1e20;
      continue;
    }
    row[1] = xi;
    if (value_only) {
      double v = cell<0>(y[r], row, m).v;
      out(r, 0) = (v >= 1e19 || !R_finite(v)) ? 1e20 : v;
      continue;
    }
    Jet<NPAR> z = cell<NPAR>(y[r], row, m);
    if (z.v >= 1e19 || !R_finite(z.v)) {
      out(r, 0) = 1e20;
      continue;
    }
    double g[NPAR], H[NPAR][NPAR];
    for (int i = 0; i < n; ++i)
      g[i] = (i == 1) ? xi * z.g[1] : z.g[i];
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (i == 1 && j == 1) H[1][1] = xi * xi * z.h[1][1] + xi * z.g[1];
        else if (i == 1) H[1][j] = xi * z.h[1][j];
        else if (j == 1) H[i][1] = xi * z.h[i][1];
        else H[i][j] = z.h[i][j];
      }
    }
    if (apply_bounded_chain && (mp_xi != 1.0 || mpp_xi != 0.0)) {
      double g1 = g[1], h11 = H[1][1];
      g[1] = g1 * mp_xi;
      for (int k = 0; k < n; ++k) {
        if (k != 1) {
          H[1][k] *= mp_xi;
          H[k][1] *= mp_xi;
        }
      }
      H[1][1] = h11 * mp_xi * mp_xi + g1 * mpp_xi;
    }
    out(r, 0) = z.v;
    for (int i = 0; i < n; ++i) out(r, i + 1) = g[i];
    int c = n + 1;
    for (int i = 0; i < n; ++i)
      for (int j = i; j < n; ++j)
        out(r, c++) = H[i][j];
  }
  return out;
}

// Log-link ZIDEGPD evaluator for models 1..6.
// Combines the stabilized DEGPD log-PMF Jet with logit(pi) in log-sum-exp space.
// Preserves zidegpd4's parameter order (lsigma, lxi, lkappa, ldelta, logitpi).
arma::mat zidegpd_loglink_eval(const Rcpp::List& pars, const std::vector<arma::mat>& X,
                               const arma::vec& y, const arma::uvec& dupid, int dcate,
                               const Rcpp::List& offsets, int m, bool value_only) {
  arma::mat eta = predictors(pars, X, dupid, dcate, offsets);
  int n = X.size(), d = n - 1;
  arma::mat out(y.n_elem, value_only ? 1 : 1 + n + n * (n + 1) / 2, arma::fill::zeros);
  for (unsigned r = 0; r < y.n_elem; ++r) {
    arma::rowvec row = eta.row(r);
    if (y[r] < 0.0 || !std::isfinite(y[r]) || !row.is_finite()) {
      out(r, 0) = 1e20;
      continue;
    }
    double xi = std::exp(row[1]);
    if (!R_finite(xi)) {
      out(r, 0) = 1e20;
      continue;
    }
    arma::rowvec eta_id(d);
    eta_id[0] = row[0];
    eta_id[1] = xi;
    if (m == 4) {
      // zidegpd4 has X3 = lkappa, X4 = ldelta; carrier 4 expects (ldelta, lkappa).
      eta_id[2] = row[3];
      eta_id[3] = row[2];
    } else {
      for (int k = 2; k < d; ++k) eta_id[k] = row[k];
    }
    Jet<NPAR> z;
    double L;
    if (value_only) {
      L = cell<0>(y[r], eta_id, m).v;
    } else {
      z = cell<NPAR>(y[r], eta_id, m);
      L = z.v;
    }
    if (L >= 1e19 || !R_finite(L)) {
      out(r, 0) = 1e20;
      continue;
    }
    double lp_eta = row[d], lpi, lompi;
    if (lp_eta >= 0.0) {
      lpi = -std::log1p(std::exp(-lp_eta));
      lompi = -lp_eta + lpi;
    } else {
      lompi = -std::log1p(std::exp(lp_eta));
      lpi = lp_eta + lompi;
    }
    // y = 0: log(pi + (1 - pi) P(0)). A log-sum-exp (with log1p) is accurate
    // while the probability is below 1/2; above it, evaluate the complement
    // log(1 - (1 - pi)(1 - P(0))) so that a tiny NLL keeps full relative
    // accuracy. Either way the NLL is never negative. When P(0) is within
    // rounding of 1 the cell NLL L = -log P(0) can come out as -1e-16 (mixture
    // log-sum-exp, quadrature); P(0) <= 1, so treat that as 1 - P(0) = 0.
    double lse = std::max(lpi, lompi - L) +
      std::log1p(std::exp(-std::abs(lpi - (lompi - L))));
    if (lse > -0.6931471805599453) {
      double log_1mP0 = L > 0.0 ? log1mexp_d(-L) : R_NegInf;
      lse = log1mexp_d(lompi + log_1mP0);
    }
    double nll0 = -lse;
    if (value_only) {
      out(r, 0) = (y[r] > 0.0) ? (-lompi + L) : nll0;
      continue;
    }
    int perm[NPAR] = {0, 1, 2, 3, 4};
    if (m == 4) { perm[2] = 3; perm[3] = 2; }
    double g[NPAR], H[NPAR][NPAR];
    for (int i = 0; i < d; ++i) {
      int pi = perm[i];
      g[i] = (i == 1) ? xi * z.g[1] : z.g[pi];
    }
    for (int i = 0; i < d; ++i) {
      int pi = perm[i];
      for (int j = 0; j < d; ++j) {
        int pj = perm[j];
        if (i == 1 && j == 1) H[1][1] = xi * xi * z.h[1][1] + xi * z.g[1];
        else if (i == 1) H[1][j] = xi * z.h[1][pj];
        else if (j == 1) H[i][1] = xi * z.h[pi][1];
        else H[i][j] = z.h[pi][pj];
      }
    }
    double pi_val = std::exp(lpi), ompi_val = std::exp(lompi);
    double pi_ompi = std::exp(lpi + lompi);
    double nll_r, gf[6], Hf[6][6];
    if (y[r] > 0.0) {
      nll_r = -lompi + L;
      for (int i = 0; i < d; ++i) {
        gf[i] = g[i];
        for (int j = i; j < d; ++j) Hf[i][j] = H[i][j];
        Hf[i][d] = 0.0;
      }
      gf[d] = pi_val;
      Hf[d][d] = pi_ompi;
    } else {
      nll_r = nll0;
      double diff = lpi - (lompi - L);
      double w1, w2;
      if (diff >= 0.0) {
        double e = std::exp(-diff);
        w1 = 1.0 / (1.0 + e);
        w2 = e * w1;
      } else {
        double e = std::exp(diff);
        w2 = 1.0 / (1.0 + e);
        w1 = e * w2;
      }
      double w12 = w1 * w2;
      for (int i = 0; i < d; ++i) {
        gf[i] = w2 * g[i];
        for (int j = i; j < d; ++j)
          Hf[i][j] = w2 * H[i][j] - w12 * g[i] * g[j];
        Hf[i][d] = -w12 * g[i];
      }
      double g_eta = pi_val * w2 - ompi_val * w1;
      gf[d] = g_eta;
      Hf[d][d] = g_eta * (w2 - pi_val);
    }
    out(r, 0) = nll_r;
    for (int i = 0; i < n; ++i) out(r, i + 1) = gf[i];
    int c = n + 1;
    for (int i = 0; i < n; ++i)
      for (int j = i; j < n; ++j)
        out(r, c++) = Hf[i][j];
  }
  return out;
}

// [[Rcpp::export]]
double degpd1id0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2,
                 const arma::mat& X3, arma::vec yvec, const arma::uvec& dupid,
                 int dcate, const Rcpp::List& offsets, double xi_max) {
  return arma::accu(degpd_identity_eval(pars, {X1, X2, X3}, yvec, dupid, dcate, offsets, 1, true).col(0));
}

// [[Rcpp::export]]
arma::mat degpd1id12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3,
                    arma::vec yvec, const arma::uvec dupid, int dcate,
                    const Rcpp::List& offsets, double xi_max) {
  return degpd_identity_eval(pars, {X1, X2, X3}, yvec, dupid, dcate, offsets, 1, false).cols(1, 9);
}
