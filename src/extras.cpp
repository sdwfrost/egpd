// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "gpd_endpoint.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace {
// relative = true follows MASS::ginv: singular values d <= tol * max(d) are
// dropped; otherwise tol is an absolute cut-off (or the LAPACK default).
arma::mat svd_std_pinv(const arma::mat& x, double tol, bool default_tol,
                       bool relative = false) {
  if (!x.is_finite()) Rcpp::stop("pinv(): matrix contains non-finite values");
  if (!default_tol && !(tol >= 0.0)) Rcpp::stop("pinv(): tolerance must be >= 0");
  if (x.is_empty()) return arma::mat(x.n_cols, x.n_rows, arma::fill::zeros);
  arma::mat U, V;
  arma::vec s;
  bool ok = arma::svd_econ(U, s, V, x, 'b', "std");
  if (!ok) ok = arma::svd(U, s, V, x, "std");
  if (!ok) Rcpp::stop("pinv(): svd failed");
  if (s.n_elem == 0) return arma::mat(x.n_cols, x.n_rows, arma::fill::zeros);
  double thresh = relative ? tol * s[0]
    : (default_tol
       ? (static_cast<double>(std::max(x.n_rows, x.n_cols)) * s[0] *
          std::numeric_limits<double>::epsilon())
       : tol);
  arma::uword count = 0;
  for (arma::uword i = 0; i < s.n_elem; ++i) {
    if (relative ? (s[i] > thresh) : (s[i] >= thresh && s[i] > 0.0)) ++count;
  }
  if (count == 0) return arma::mat(x.n_cols, x.n_rows, arma::fill::zeros);
  arma::vec inv_s = 1.0 / s.head(count);
  return (V.cols(0, count - 1) * arma::diagmat(inv_s)) * U.cols(0, count - 1).t();
}
}

// [[Rcpp::export]]
arma::mat armapinv(arma::mat x, double tol)
{
  // Armadillo interprets zero as its automatic tolerance, not an absolute zero.
  if (tol == 0.0) return svd_std_pinv(x, tol, false);
  try {
    if (tol < 0) {
      return arma::pinv(x);
    } else {
      return arma::pinv(x, tol);
    }
  } catch (...) {
    return svd_std_pinv(x, tol, tol < 0);
  }
}

// [[Rcpp::export]]
arma::mat armaginv(arma::mat x, double tol)
{
  return svd_std_pinv(x, tol, false, true);
}

// Refine eta + logit(p) only near the hurdle jump, via p*exp(eta) - (1-p).
// A direct subtraction of the rounded atom or its log can give the wrong sign.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector hurdle_logodds_gap_cpp(Rcpp::NumericVector eta, double p) {
  const double ln2 = 0.69314718055994530942, ln2lo = 2.31904681384629956e-17;
  const double c = 1 - p, clo = (1 - c) - p;
  Rcpp::NumericVector ans(eta.size());
  for (R_xlen_t i = 0; i < eta.size(); ++i) {
    if (!(p > 0 && p < 1 && std::abs(eta[i]) < 746))
      Rcpp::stop("Expected an interior probability and a finite link near the zero atom.");
    const int m = static_cast<int>(std::round(eta[i] / ln2));
    const double r0 = std::fma(-m, ln2, eta[i]), t = -m * ln2lo;
    const double r = r0 + t, rl = (r0 - r) + t;
    // Compensated Taylor series for exp(r + rl), |r| <= log(2)/2.
    double hi = 1, lo = 0;
    for (int j = 24; j > 0; --j) {
      const double v = r * hi, vl = std::fma(r, hi, -v) + r * lo + rl * hi;
      const double q = v / j, ql = (std::fma(-q, j, v) + vl) / j;
      hi = 1 + q;
      lo = ((1 - hi) + q) + ql;
    }
    const double a = std::ldexp(p, m);
    ans[i] = std::log1p(((std::fma(a, hi, -c) - clo) + a * lo) / c);
  }
  return ans;
}

// 1 + xi*x/sigma near a finite upper endpoint. Multiplying mantissas in a
// fused add retains the sign even when -sigma/xi rounds to x. Normalization
// also protects the residual when sigma itself is subnormal.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector gpd_endpoint_gap_cpp(Rcpp::NumericVector x, Rcpp::NumericVector sigma,
                                         Rcpp::NumericVector xi) {
  if (sigma.size() != x.size() || xi.size() != x.size())
    Rcpp::stop("Incompatible endpoint vector lengths.");
  Rcpp::NumericVector gap(x.size());
  for (R_xlen_t i = 0; i < x.size(); ++i)
    gap[i] = egpd_endpoint_gap(x[i], sigma[i], xi[i]);
  return gap;
}
