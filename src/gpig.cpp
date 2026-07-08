// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

// Generalised Poisson-inverse Gaussian (GPIG) of Zhu & Joe (2009,
// Stat. Prob. Lett. 79:1695-1703) and its zero-inflated variant (ZIGPIG) for
// the native GAM fitter (egpd()) and the standalone d/p/q/r functions.
//
// GPIG(a, b, c) is defined by its pgf
//   G(s) = exp{ b[(1-c)^a - (1-cs)^a] },   0 < a <= 1, b > 0, 0 <= c <= 1,
// with Poisson (a=1), PIG (a=1/2) and discrete-stable (c=1) as special cases.
// It has no closed-form pmf but a stable forward recursion (Zhu & Joe eqs
// 3.2-3.3):
//   r_1 = (1-a)c,  r_{j+1} = ((j-1+a)/(j+1)) c r_j,
//   p_0 = exp(b[(1-c)^a - 1]),  p_1 = abc p_0,
//   p_{k+1} = ( abc p_k + sum_{j=1}^k j r_{k+1-j} p_j ) / (k+1).
// All terms are non-negative (r_j >= 0 for a < 1), so the recursion has no
// cancellation. The analytic first derivatives dp_k/d(a,b,c) follow eqs
// 3.4-3.5, EXCEPT the two base cases dp0/da and dp0/db, which are misprinted in
// the paper; the correct forms (verified against finite differences to ~1e-11)
// are used here:
//   dp0/da = b (1-c)^a log(1-c) p0,   dp0/db = ((1-c)^a - 1) p0,
//   dp0/dc = -a b (1-c)^{a-1} p0.
//
// Two parameterisations are supported for the fitting surfaces (selected by the
// `param` flag); a Jacobian maps the (a,b,c) score to the linear predictors eta.
//   param 0 (native): eta = (logit a, log b, logit c)   -- paper-faithful.
//   param 1 (mean)  : eta = (log mu, logit a, logit c), b = mu(1-c)^(1-a)/(ac)
//                     -- models the mean directly, consistent with the EGPD
//                     families.
//
// The score is exact; the Hessian is the empirical-Fisher / BHHH approximation
// (outer product of the per-observation NLL gradient), matching the estimator
// used by Zhu & Joe for standard errors. It is positive semi-definite; the
// fitter perturbs it when required.
//
// d12 returns per-observation derivatives of the NEGATIVE log-likelihood wrt
// the linear predictors, packed to match .gH2 / .gH3:
//   npar=3: cols [g1,g2,g3, H11,H12,H13,H22,H23,H33]                 (nobs x 9)
//   npar=4: cols [g1,g2,g3,g4, H11,H12,H13,H14,H22,H23,H24,H33,H34,H44] (x 14)

// ------------------------------------------------------------------ core ------

// Forward recursion to p_y with analytic derivatives. Returns p_y and its
// partials wrt (a,b,c) via out-parameters. Handles y = 0 through the base case.
static inline void gpig_derivs(int y, double a, double b, double c,
                               double& py, double& dpa, double& dpb, double& dpc)
{
  double omc = 1.0 - c;
  double omca = std::pow(omc, a);            // (1-c)^a
  double p0 = std::exp(b * (omca - 1.0));
  // corrected base cases (paper eq 3.5 first line is misprinted for a, b)
  double dp0a = b * omca * std::log(omc) * p0;
  double dp0b = (omca - 1.0) * p0;
  double dp0c = -a * b * std::pow(omc, a - 1.0) * p0;

  if (y == 0) {
    py = p0; dpa = dp0a; dpb = dp0b; dpc = dp0c;
    return;
  }

  std::vector<double> p(y + 1), r(y + 1);
  std::vector<double> dra(y + 1, 0.0), drc(y + 1, 0.0);
  std::vector<double> pa(y + 1), pb(y + 1), pc(y + 1);

  // r_1..r_{y-1} and their derivatives (eqs 3.2, 3.4). r[j] holds r_j.
  r[1] = (1.0 - a) * c;
  dra[1] = -c;
  drc[1] = 1.0 - a;
  for (int j = 1; j <= y - 1; j++) {
    double f = (j - 1.0 + a) / (j + 1.0);
    r[j + 1]   = f * c * r[j];
    dra[j + 1] = f * c * dra[j] + c * r[j] / (j + 1.0);
    drc[j + 1] = f * c * drc[j] + f * r[j];
  }

  // p_0, p_1 and their derivatives (eqs 3.3, 3.5).
  p[0] = p0; pa[0] = dp0a; pb[0] = dp0b; pc[0] = dp0c;
  p[1]  = a * b * c * p0;
  pa[1] = a * b * c * dp0a + b * c * p0;
  pb[1] = a * b * c * dp0b + a * c * p0;
  pc[1] = a * b * c * dp0c + a * b * p0;

  for (int k = 1; k <= y - 1; k++) {
    double s = 0.0, sa = 0.0, sb = 0.0, sc = 0.0;
    for (int j = 1; j <= k; j++) {
      int rj = k + 1 - j;
      s  += j * r[rj] * p[j];
      sa += j * (p[j] * dra[rj] + r[rj] * pa[j]);
      sb += j * (r[rj] * pb[j]);
      sc += j * (p[j] * drc[rj] + r[rj] * pc[j]);
    }
    double km = k + 1.0;
    p[k + 1]  = (a * b * c * p[k]  + s)  / km;
    pa[k + 1] = (a * b * c * pa[k] + b * c * p[k] + sa) / km;
    pb[k + 1] = (a * b * c * pb[k] + a * c * p[k] + sb) / km;
    pc[k + 1] = (a * b * c * pc[k] + a * b * p[k] + sc) / km;
  }

  py = p[y]; dpa = pa[y]; dpb = pb[y]; dpc = pc[y];
}

// Light recursion: p_y only (no derivatives), for d0 / pmf.
static inline double gpig_p(int y, double a, double b, double c)
{
  double omc = 1.0 - c;
  double p0 = std::exp(b * (std::pow(omc, a) - 1.0));
  if (y == 0) return p0;
  std::vector<double> p(y + 1), r(y + 1);
  r[1] = (1.0 - a) * c;
  for (int j = 1; j <= y - 1; j++)
    r[j + 1] = ((j - 1.0 + a) / (j + 1.0)) * c * r[j];
  p[0] = p0;
  p[1] = a * b * c * p0;
  for (int k = 1; k <= y - 1; k++) {
    double s = 0.0;
    for (int j = 1; j <= k; j++) s += j * r[k + 1 - j] * p[j];
    p[k + 1] = (a * b * c * p[k] + s) / (k + 1.0);
  }
  return p[y];
}

static inline void eta_to_abc(int param, double e1, double e2, double e3,
                              double& a, double& b, double& c)
{
  if (param == 0) {                          // native (logit a, log b, logit c)
    a = 1.0 / (1.0 + std::exp(-e1));
    b = std::exp(e2);
    c = 1.0 / (1.0 + std::exp(-e3));
  } else {                                   // mean (log mu, logit a, logit c)
    double mu = std::exp(e1);
    a = 1.0 / (1.0 + std::exp(-e2));
    c = 1.0 / (1.0 + std::exp(-e3));
    b = mu * std::pow(1.0 - c, 1.0 - a) / (a * c);
  }
}

// Map an (a,b,c)-space gradient (Fa,Fb,Fc) of a scalar to the eta-gradient.
static inline void map_grad(int param, double a, double b, double c,
                            double Fa, double Fb, double Fc, double* ge)
{
  if (param == 0) {
    ge[0] = Fa * a * (1.0 - a);
    ge[1] = Fb * b;
    ge[2] = Fc * c * (1.0 - c);
  } else {
    double dbde1 = b;                                        // d b / d eta_mu
    double dbde2 = b * (-std::log(1.0 - c) - 1.0 / a) * a * (1.0 - a);
    double dbde3 = b * (-(1.0 - a) / (1.0 - c) - 1.0 / c) * c * (1.0 - c);
    ge[0] = Fb * dbde1;
    ge[1] = Fa * a * (1.0 - a) + Fb * dbde2;
    ge[2] = Fc * c * (1.0 - c) + Fb * dbde3;
  }
}

// Pack per-obs NLL gradient g (length npar) + BHHH Hessian (= g g^T) into a row.
static inline void pack_bhhh(arma::mat& out, int j, const double* g, int npar)
{
  for (int i = 0; i < npar; i++) out(j, i) = g[i];
  int col = npar;
  for (int i = 0; i < npar; i++)
    for (int k = i; k < npar; k++)
      out(j, col++) = g[i] * g[k];
}

// ----------------------------------------------------------- GPIG (3 par) -----

// Shared d0 for both parameterisations (param flag).
static double gpig_d0_impl(int param, const Rcpp::List& pars,
                           const arma::mat& X1, const arma::mat& X2,
                           const arma::mat& X3, arma::vec yvec,
                           const arma::uvec& dupid, int dcate,
                           const Rcpp::List& offsets)
{
  arma::vec e1 = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec e2 = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec e3 = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = yvec.size();
  if (dcate == 1) { e1 = e1.elem(dupid); e2 = e2.elem(dupid); e3 = e3.elem(dupid); }
  {
    arma::vec o0 = Rcpp::as<arma::vec>(offsets[0]); if (o0.n_elem > 0) e1 += o0;
    arma::vec o1 = Rcpp::as<arma::vec>(offsets[1]); if (o1.n_elem > 0) e2 += o1;
    arma::vec o2 = Rcpp::as<arma::vec>(offsets[2]); if (o2.n_elem > 0) e3 += o2;
  }
  double nllh = 0.0;
  for (int j = 0; j < nobs; j++) {
    double a, b, c; eta_to_abc(param, e1[j], e2[j], e3[j], a, b, c);
    int y = (int)(yvec[j] + 0.5);
    double py = gpig_p(y, a, b, c);
    nllh += -std::log(py);
  }
  return nllh;
}

static arma::mat gpig_d12_impl(int param, const Rcpp::List& pars,
                               const arma::mat& X1, const arma::mat& X2,
                               const arma::mat& X3, arma::vec yvec,
                               const arma::uvec& dupid, int dcate,
                               const Rcpp::List& offsets)
{
  arma::vec e1 = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec e2 = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec e3 = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = yvec.size();
  arma::mat out(nobs, 9);
  if (dcate == 1) { e1 = e1.elem(dupid); e2 = e2.elem(dupid); e3 = e3.elem(dupid); }
  {
    arma::vec o0 = Rcpp::as<arma::vec>(offsets[0]); if (o0.n_elem > 0) e1 += o0;
    arma::vec o1 = Rcpp::as<arma::vec>(offsets[1]); if (o1.n_elem > 0) e2 += o1;
    arma::vec o2 = Rcpp::as<arma::vec>(offsets[2]); if (o2.n_elem > 0) e3 += o2;
  }
  for (int j = 0; j < nobs; j++) {
    double a, b, c; eta_to_abc(param, e1[j], e2[j], e3[j], a, b, c);
    int y = (int)(yvec[j] + 0.5);
    double py, dpa, dpb, dpc;
    gpig_derivs(y, a, b, c, py, dpa, dpb, dpc);
    double ge[3];
    map_grad(param, a, b, c, dpa / py, dpb / py, dpc / py, ge);
    double g[3] = { -ge[0], -ge[1], -ge[2] };   // NLL gradient
    pack_bhhh(out, j, g, 3);
  }
  return out;
}

// [[Rcpp::export]]
double gpig1d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2,
               const arma::mat& X3, arma::vec yvec, const arma::uvec& dupid,
               int dcate, const Rcpp::List& offsets) {
  return gpig_d0_impl(1, pars, X1, X2, X3, yvec, dupid, dcate, offsets);
}
// [[Rcpp::export]]
arma::mat gpig1d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2,
                   arma::mat X3, arma::vec yvec, const arma::uvec dupid,
                   int dcate, const Rcpp::List& offsets) {
  return gpig_d12_impl(1, pars, X1, X2, X3, yvec, dupid, dcate, offsets);
}
// [[Rcpp::export]]
double gpignat1d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2,
                  const arma::mat& X3, arma::vec yvec, const arma::uvec& dupid,
                  int dcate, const Rcpp::List& offsets) {
  return gpig_d0_impl(0, pars, X1, X2, X3, yvec, dupid, dcate, offsets);
}
// [[Rcpp::export]]
arma::mat gpignat1d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2,
                      arma::mat X3, arma::vec yvec, const arma::uvec dupid,
                      int dcate, const Rcpp::List& offsets) {
  return gpig_d12_impl(0, pars, X1, X2, X3, yvec, dupid, dcate, offsets);
}

// --------------------------------------------------------- ZIGPIG (4 par) -----

static double zigpig_d0_impl(int param, const Rcpp::List& pars,
                             const arma::mat& X1, const arma::mat& X2,
                             const arma::mat& X3, const arma::mat& X4,
                             arma::vec yvec, const arma::uvec& dupid, int dcate,
                             const Rcpp::List& offsets)
{
  arma::vec e1 = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec e2 = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec e3 = X3 * Rcpp::as<arma::vec>(pars[2]);
  arma::vec e4 = X4 * Rcpp::as<arma::vec>(pars[3]);
  int nobs = yvec.size();
  if (dcate == 1) {
    e1 = e1.elem(dupid); e2 = e2.elem(dupid);
    e3 = e3.elem(dupid); e4 = e4.elem(dupid);
  }
  {
    arma::vec o0 = Rcpp::as<arma::vec>(offsets[0]); if (o0.n_elem > 0) e1 += o0;
    arma::vec o1 = Rcpp::as<arma::vec>(offsets[1]); if (o1.n_elem > 0) e2 += o1;
    arma::vec o2 = Rcpp::as<arma::vec>(offsets[2]); if (o2.n_elem > 0) e3 += o2;
    arma::vec o3 = Rcpp::as<arma::vec>(offsets[3]); if (o3.n_elem > 0) e4 += o3;
  }
  double nllh = 0.0;
  for (int j = 0; j < nobs; j++) {
    double a, b, c; eta_to_abc(param, e1[j], e2[j], e3[j], a, b, c);
    double pi = 1.0 / (1.0 + std::exp(-e4[j]));
    int y = (int)(yvec[j] + 0.5);
    double logf;
    if (y == 0) {
      double p0 = gpig_p(0, a, b, c);
      logf = std::log(pi + (1.0 - pi) * p0);
    } else {
      logf = std::log1p(-pi) + std::log(gpig_p(y, a, b, c));
    }
    nllh += -logf;
  }
  return nllh;
}

static arma::mat zigpig_d12_impl(int param, const Rcpp::List& pars,
                                 const arma::mat& X1, const arma::mat& X2,
                                 const arma::mat& X3, const arma::mat& X4,
                                 arma::vec yvec, const arma::uvec& dupid,
                                 int dcate, const Rcpp::List& offsets)
{
  arma::vec e1 = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec e2 = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec e3 = X3 * Rcpp::as<arma::vec>(pars[2]);
  arma::vec e4 = X4 * Rcpp::as<arma::vec>(pars[3]);
  int nobs = yvec.size();
  arma::mat out(nobs, 14);
  if (dcate == 1) {
    e1 = e1.elem(dupid); e2 = e2.elem(dupid);
    e3 = e3.elem(dupid); e4 = e4.elem(dupid);
  }
  {
    arma::vec o0 = Rcpp::as<arma::vec>(offsets[0]); if (o0.n_elem > 0) e1 += o0;
    arma::vec o1 = Rcpp::as<arma::vec>(offsets[1]); if (o1.n_elem > 0) e2 += o1;
    arma::vec o2 = Rcpp::as<arma::vec>(offsets[2]); if (o2.n_elem > 0) e3 += o2;
    arma::vec o3 = Rcpp::as<arma::vec>(offsets[3]); if (o3.n_elem > 0) e4 += o3;
  }
  for (int j = 0; j < nobs; j++) {
    double a, b, c; eta_to_abc(param, e1[j], e2[j], e3[j], a, b, c);
    double pi = 1.0 / (1.0 + std::exp(-e4[j]));
    double pip = pi * (1.0 - pi);                      // d pi / d eta_pi
    int y = (int)(yvec[j] + 0.5);

    double py, dpa, dpb, dpc;
    gpig_derivs(y, a, b, c, py, dpa, dpb, dpc);
    double g[4];

    if (y > 0) {
      double ge[3];
      map_grad(param, a, b, c, dpa / py, dpb / py, dpc / py, ge);
      g[0] = -ge[0]; g[1] = -ge[1]; g[2] = -ge[2];
      g[3] = pi;                                        // -(score_pi = -pi)
    } else {
      double p0 = py;                                   // gpig pmf at 0
      double D = pi + (1.0 - pi) * p0;
      // d log D / d(a,b,c) = (1-pi) dp0/d(a,b,c) / D, then map to eta.
      double f = (1.0 - pi) / D;
      double ge[3];
      map_grad(param, a, b, c, f * dpa, f * dpb, f * dpc, ge);
      double score_pi = pip * (1.0 - p0) / D;           // d log D / d eta_pi
      g[0] = -ge[0]; g[1] = -ge[1]; g[2] = -ge[2];
      g[3] = -score_pi;
    }
    pack_bhhh(out, j, g, 4);
  }
  return out;
}

// [[Rcpp::export]]
double zigpig1d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2,
                 const arma::mat& X3, const arma::mat& X4, arma::vec yvec,
                 const arma::uvec& dupid, int dcate, const Rcpp::List& offsets) {
  return zigpig_d0_impl(1, pars, X1, X2, X3, X4, yvec, dupid, dcate, offsets);
}
// [[Rcpp::export]]
arma::mat zigpig1d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2,
                     arma::mat X3, arma::mat X4, arma::vec yvec,
                     const arma::uvec dupid, int dcate, const Rcpp::List& offsets) {
  return zigpig_d12_impl(1, pars, X1, X2, X3, X4, yvec, dupid, dcate, offsets);
}
// [[Rcpp::export]]
double zigpignat1d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2,
                    const arma::mat& X3, const arma::mat& X4, arma::vec yvec,
                    const arma::uvec& dupid, int dcate, const Rcpp::List& offsets) {
  return zigpig_d0_impl(0, pars, X1, X2, X3, X4, yvec, dupid, dcate, offsets);
}
// [[Rcpp::export]]
arma::mat zigpignat1d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2,
                        arma::mat X3, arma::mat X4, arma::vec yvec,
                        const arma::uvec dupid, int dcate, const Rcpp::List& offsets) {
  return zigpig_d12_impl(0, pars, X1, X2, X3, X4, yvec, dupid, dcate, offsets);
}

// ------------------------------------------------------ pmf for d/p/q/r --------

// pmf vector p_0..p_n in native (a,b,c) parameterisation.
// [[Rcpp::export]]
Rcpp::NumericVector gpig_pmf_cpp(int n, double a, double b, double c)
{
  Rcpp::NumericVector out(n + 1);
  double omc = 1.0 - c;
  double p0 = std::exp(b * (std::pow(omc, a) - 1.0));
  out[0] = p0;
  if (n == 0) return out;
  std::vector<double> p(n + 1), r(n + 1);
  r[1] = (1.0 - a) * c;
  for (int j = 1; j <= n - 1; j++)
    r[j + 1] = ((j - 1.0 + a) / (j + 1.0)) * c * r[j];
  p[0] = p0;
  p[1] = a * b * c * p0;
  out[1] = p[1];
  for (int k = 1; k <= n - 1; k++) {
    double s = 0.0;
    for (int j = 1; j <= k; j++) s += j * r[k + 1 - j] * p[j];
    p[k + 1] = (a * b * c * p[k] + s) / (k + 1.0);
    out[k + 1] = p[k + 1];
  }
  return out;
}
