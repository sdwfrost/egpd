// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <limits>

// Bell distribution of Castellares, Ferrari & Lemonte (2018, Applied
// Mathematical Modelling 56:172-185) and its zero-inflated variant (ZIBell) for
// the native GAM fitter (egpd()) and the standalone d/p/q/r functions.
//
// Bell(theta), theta > 0, is a one-parameter overdispersed count family with
//   Pr(Y = y) = theta^y exp(1 - e^theta) B_y / y!,   y = 0, 1, 2, ...,
// where B_y are the Bell numbers (B_y/y! = [x^y] exp(e^x - 1)). It belongs to
// the one-parameter exponential family, has mean theta*e^theta and variance
// theta(1+theta)e^theta (dispersion index 1+theta > 1, always overdispersed),
// and the ML estimator is theta.hat = W0(ybar) (Lambert W). Everything is
// closed form (no recursion is needed for the likelihood, unlike PIG/GPIG).
//
// Two parameterisations are supported for the fitting surfaces (selected by the
// `param` flag); the chain rule maps the theta-score/Hessian to eta.
//   param 0 (native): eta = log theta                       -- paper Def. 1.
//   param 1 (mean)  : eta = log mu, theta = W0(mu)           -- paper section 4,
//                     models the mean directly (E[Y] = mu), consistent with the
//                     EGPD families.
//
// Because the log-likelihood is closed form, the score AND the observed-
// information Hessian are exact (as in PIG via PigJet, and unlike the BHHH
// approximation used for GPIG).
//
// d12 returns per-observation derivatives of the NEGATIVE log-likelihood wrt
// the linear predictors, packed to match .gH1 / .gH2 (upper triangle,
// row-major):
//   npar=1: cols [g1, H11]                           (nobs x 2)
//   npar=2: cols [g1, g2, H11, H12, H22]             (nobs x 5)

// ------------------------------------------------------ Bell ratios / W0 ------

// log(B_n / n!) for n = 0..nmax, by the stable log-space recursion
//   a_n = (1/n) sum_{j=0}^{n-1} a_{n-1-j} / j!,   a_0 = 1,
// i.e. la_n = -log(n) + logsumexp_j ( la_{n-1-j} - log j! ).
//
// The table does not depend on theta, so it is built once per session and grown on
// demand (it used to be rebuilt, at O(n^2) cost, on every likelihood evaluation and
// every d/p/q/r call). Building is also O(n) rather than O(n^2): the summand
// t_j = la_{n-1-j} - log j! is concave in j -- its increments
// log(a_{m-1}/a_m) - log j decrease, because m B_{m-1}/B_m is non-decreasing in m
// (Bell numbers are log-convex) -- so it is unimodal with its peak near j = W(n),
// and the sum can stop once the terms are past the peak, more than 60 nats below it
// and falling by more than a factor e^2 per step. Everything dropped is below
// e^-60 (~1e-26) relative to the sum, far beneath double precision, so the values
// are those of the full sum (checked against exact Bell numbers).
static std::vector<double> bell_la_cache(1, 0.0);      // la_0 = log(B_0/0!) = 0
static std::vector<double> bell_lfact_cache(1, 0.0);   // log j!, grown as needed

static const std::vector<double>& bell_log_ratio(int nmax)
{
  if (nmax < 0) Rcpp::stop("Bell support must be non-negative");
  if (nmax > 100000000) Rcpp::stop("Bell support beyond 1e8 is not supported");
  std::vector<double>& la = bell_la_cache;
  std::vector<double>& lf = bell_lfact_cache;
  int have = (int)la.size() - 1;
  if (nmax <= have) return la;
  la.reserve(nmax + 1);
  for (int n = have + 1; n <= nmax; n++) {
    if ((n & 4095) == 0) Rcpp::checkUserInterrupt();
    double tmax = -std::numeric_limits<double>::infinity(), prev = tmax;
    int J = n - 1;
    for (int j = 0; j <= n - 1; j++) {
      while ((int)lf.size() <= j) lf.push_back(std::lgamma((double)lf.size() + 1.0));
      double t = la[n - 1 - j] - lf[j];
      if (t > tmax) tmax = t;
      if (j > 0 && t < tmax - 60.0 && t - prev < -2.0) { J = j; break; }
      prev = t;
    }
    double s = 0.0;
    for (int j = 0; j <= J; j++) s += std::exp(la[n - 1 - j] - lf[j] - tmax);
    la.push_back(-std::log((double)n) + tmax + std::log(s));
  }
  return la;
}

static inline double bell_logsumexp2(double a, double b)
{
  if (a == b) return a + M_LN2;
  double m = (a > b) ? a : b;
  if (!std::isfinite(m)) return m;
  return m + R::log1pexp(-std::fabs(a - b));
}

// Principal branch of the Lambert W function for x >= 0 (Halley iteration).
static inline double lambertW0(double x)
{
  if (std::isnan(x)) return x;
  if (x == 0.0) return 0.0;                       // W0(0) = 0
  if (x < 0.0) return R_NaN;
  if (!std::isfinite(x)) return R_PosInf;
  double w = (x < 1.0) ? x : std::log(x);
  if (x > 3.0) w = std::log(x) - std::log(std::log(x));
  for (int i = 0; i < 60; i++) {
    double ew = std::exp(w);
    double f  = w * ew - x;
    double wp1 = w + 1.0;
    double dw = f / (ew * wp1 - (w + 2.0) * f / (2.0 * wp1));
    w -= dw;
    if (std::fabs(dw) <= 1e-15 * (1.0 + std::fabs(w))) break;
  }
  return w;
}

static inline double bell_theta(int param, double eta)
{
  return param == 0 ? std::exp(eta) : lambertW0(std::exp(eta));
}

static inline double bell_logtheta(int param, double eta, double theta)
{
  // log(mu) = log(theta) + theta, even if exp(eta) underflows.
  return param == 0 ? eta : (theta < 1.0 ? eta - theta : std::log(theta));
}

// Differentiate in log(theta), not theta: y/theta^2 can overflow even
// when both link-scale derivatives are representable.
static inline void bell_link_d12(int param, double eta, double theta, int y,
                                 double& g, double& h)
{
  const double mu = param == 0 ? theta * std::exp(theta) : std::exp(eta);
  const double t = 1.0 + theta;
  if (param == 0) {
    g = mu - y;
    h = mu * t;
  } else {
    g = (mu - y) / t;
    h = mu / t - g * (theta / t) / t;
  }
}

// ---------------------------------------------------------- Bell (1 par) ------

static double bell_d0_impl(int param, const Rcpp::List& pars, const arma::mat& X1,
                           arma::vec yvec, const arma::uvec& dupid, int dcate,
                           const Rcpp::List& offsets)
{
  arma::vec e1 = X1 * Rcpp::as<arma::vec>(pars[0]);
  int nobs = yvec.size();
  if (dcate == 1) e1 = e1.elem(dupid);
  { arma::vec o0 = Rcpp::as<arma::vec>(offsets[0]); if (o0.n_elem > 0) e1 += o0; }

  int ymax = 0;
  for (int j = 0; j < nobs; j++) { int y = (int)(yvec[j] + 0.5); if (y > ymax) ymax = y; }
  const std::vector<double>& la = bell_log_ratio(ymax);

  double nllh = 0.0;
  for (int j = 0; j < nobs; j++) {
    double theta = bell_theta(param, e1[j]);
    int y = (int)(yvec[j] + 0.5);
    double logf = y * bell_logtheta(param, e1[j], theta) - std::expm1(theta) + la[y];
    nllh += -logf;
  }
  return nllh;
}

static arma::mat bell_d12_impl(int param, const Rcpp::List& pars, const arma::mat& X1,
                               arma::vec yvec, const arma::uvec& dupid, int dcate,
                               const Rcpp::List& offsets)
{
  arma::vec e1 = X1 * Rcpp::as<arma::vec>(pars[0]);
  int nobs = yvec.size();
  arma::mat out(nobs, 2);
  if (dcate == 1) e1 = e1.elem(dupid);
  { arma::vec o0 = Rcpp::as<arma::vec>(offsets[0]); if (o0.n_elem > 0) e1 += o0; }

  for (int j = 0; j < nobs; j++) {
    int y = (int)(yvec[j] + 0.5);
    bell_link_d12(param, e1[j], bell_theta(param, e1[j]), y, out(j, 0), out(j, 1));
  }
  return out;
}

// [[Rcpp::export]]
double bell1d0(const Rcpp::List& pars, const arma::mat& X1, arma::vec yvec,
               const arma::uvec& dupid, int dcate, const Rcpp::List& offsets) {
  return bell_d0_impl(1, pars, X1, yvec, dupid, dcate, offsets);
}
// [[Rcpp::export]]
arma::mat bell1d12(const Rcpp::List& pars, arma::mat X1, arma::vec yvec,
                   const arma::uvec dupid, int dcate, const Rcpp::List& offsets) {
  return bell_d12_impl(1, pars, X1, yvec, dupid, dcate, offsets);
}
// [[Rcpp::export]]
double bellnat1d0(const Rcpp::List& pars, const arma::mat& X1, arma::vec yvec,
                  const arma::uvec& dupid, int dcate, const Rcpp::List& offsets) {
  return bell_d0_impl(0, pars, X1, yvec, dupid, dcate, offsets);
}
// [[Rcpp::export]]
arma::mat bellnat1d12(const Rcpp::List& pars, arma::mat X1, arma::vec yvec,
                      const arma::uvec dupid, int dcate, const Rcpp::List& offsets) {
  return bell_d12_impl(0, pars, X1, yvec, dupid, dcate, offsets);
}

// -------------------------------------------------------- ZIBell (2 par) ------

static double zibell_d0_impl(int param, const Rcpp::List& pars, const arma::mat& X1,
                             const arma::mat& X2, arma::vec yvec,
                             const arma::uvec& dupid, int dcate,
                             const Rcpp::List& offsets)
{
  arma::vec e1 = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec e2 = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = yvec.size();
  if (dcate == 1) { e1 = e1.elem(dupid); e2 = e2.elem(dupid); }
  {
    arma::vec o0 = Rcpp::as<arma::vec>(offsets[0]); if (o0.n_elem > 0) e1 += o0;
    arma::vec o1 = Rcpp::as<arma::vec>(offsets[1]); if (o1.n_elem > 0) e2 += o1;
  }
  int ymax = 0;
  for (int j = 0; j < nobs; j++) { int y = (int)(yvec[j] + 0.5); if (y > ymax) ymax = y; }
  const std::vector<double>& la = bell_log_ratio(ymax);

  double nllh = 0.0;
  for (int j = 0; j < nobs; j++) {
    double theta = bell_theta(param, e1[j]);
    double eta_pi = e2[j];
    double log_pi = -R::log1pexp(-eta_pi);
    double log1m_pi = -R::log1pexp(eta_pi);
    int y = (int)(yvec[j] + 0.5);
    double lp0 = -std::expm1(theta);
    double logf;
    if (y == 0) {
      logf = lp0 > -M_LN2 ? std::log1p(std::exp(log1m_pi) * std::expm1(lp0))
                           : bell_logsumexp2(log_pi, log1m_pi + lp0);
    } else {
      logf = log1m_pi + y * bell_logtheta(param, e1[j], theta) + lp0 + la[y];
    }
    nllh += -logf;
  }
  return nllh;
}

static arma::mat zibell_d12_impl(int param, const Rcpp::List& pars, const arma::mat& X1,
                                 const arma::mat& X2, arma::vec yvec,
                                 const arma::uvec& dupid, int dcate,
                                 const Rcpp::List& offsets)
{
  arma::vec e1 = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec e2 = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = yvec.size();
  arma::mat out(nobs, 5);
  if (dcate == 1) { e1 = e1.elem(dupid); e2 = e2.elem(dupid); }
  {
    arma::vec o0 = Rcpp::as<arma::vec>(offsets[0]); if (o0.n_elem > 0) e1 += o0;
    arma::vec o1 = Rcpp::as<arma::vec>(offsets[1]); if (o1.n_elem > 0) e2 += o1;
  }
  for (int j = 0; j < nobs; j++) {
    double theta = bell_theta(param, e1[j]);
    double eta_pi = e2[j];
    double log_pi = -R::log1pexp(-eta_pi);
    double log1m_pi = -R::log1pexp(eta_pi);
    double pi = std::exp(log_pi);
    double one_m_pi = std::exp(log1m_pi);
    double pip = std::exp(log_pi + log1m_pi);
    int y = (int)(yvec[j] + 0.5);
    if (y > 0) {
      bell_link_d12(param, e1[j], theta, y, out(j, 0), out(j, 2));
      out(j, 1) = pi;
      out(j, 3) = 0.0;
      out(j, 4) = pip;
    } else {
      double lp0 = -std::expm1(theta);
      double delta = eta_pi - lp0;
      double log_w0 = -R::log1pexp(delta);
      double log_w1 = -R::log1pexp(-delta);
      double w0 = (delta >= 0.0) ? std::exp(log_w0) : (1.0 - std::exp(log_w1));
      double w1 = (delta >= 0.0) ? (1.0 - w0) : std::exp(log_w1);
      double w01 = std::exp(log_w0 + log_w1);
      double w0_eth   = std::isfinite(delta) ? std::exp(log_w0 + theta) : 0.0;
      double w01_eth  = std::isfinite(delta) ? std::exp(log_w0 + log_w1 + theta) : 0.0;
      double w01_eth2 = std::isfinite(delta) ? std::exp(log_w0 + log_w1 + 2.0 * theta) : 0.0;
      double t = param == 0 ? 1.0 : 1.0 + theta;
      double dth = theta / t, d2th = dth / (t * t);
      out(j, 0) = w0_eth * dth;
      out(j, 2) = (w0_eth - w01_eth2) * dth * dth + w0_eth * d2th;
      out(j, 3) = -w01_eth * dth;
      if (lp0 > -M_LN2) {
        // pi - w1 cancels when a Bell zero is nearly certain.
        double d = std::expm1(lp0);
        out(j, 1) = pip * d / (1.0 + one_m_pi * d);
        out(j, 4) = out(j, 1) * (one_m_pi - pi + out(j, 1));
      } else {
        out(j, 1) = (eta_pi >= 0.0) ? (w0 - one_m_pi) : (pi - w1);
        out(j, 4) = pip - w01;
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
double zibell1d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2,
                 arma::vec yvec, const arma::uvec& dupid, int dcate,
                 const Rcpp::List& offsets) {
  return zibell_d0_impl(1, pars, X1, X2, yvec, dupid, dcate, offsets);
}
// [[Rcpp::export]]
arma::mat zibell1d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2,
                     arma::vec yvec, const arma::uvec dupid, int dcate,
                     const Rcpp::List& offsets) {
  return zibell_d12_impl(1, pars, X1, X2, yvec, dupid, dcate, offsets);
}
// [[Rcpp::export]]
double zibellnat1d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2,
                    arma::vec yvec, const arma::uvec& dupid, int dcate,
                    const Rcpp::List& offsets) {
  return zibell_d0_impl(0, pars, X1, X2, yvec, dupid, dcate, offsets);
}
// [[Rcpp::export]]
arma::mat zibellnat1d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2,
                        arma::vec yvec, const arma::uvec dupid, int dcate,
                        const Rcpp::List& offsets) {
  return zibell_d12_impl(0, pars, X1, X2, yvec, dupid, dcate, offsets);
}

// ------------------------------------------------------ pmf for d/p/q/r --------

// pmf vector p_0..p_n in the native theta parameterisation.
// [[Rcpp::export]]
Rcpp::NumericVector bell_pmf_cpp(int n, double theta)
{
  if (n < 0) Rcpp::stop("Bell support must be non-negative");
  Rcpp::NumericVector out(n + 1);
  const std::vector<double>& la = bell_log_ratio(n);
  double base = -std::expm1(theta);               // log Pr(Y=0) = 1 - e^theta
  double lth = std::log(theta);
  Rcpp::NumericVector logs(n + 1);
  for (int k = 0; k <= n; k++) {
    logs[k] = k * lth + base + la[k];
    out[k] = std::exp(logs[k]);
  }
  out.attr("logpmf") = logs;
  return out;
}

// Lambert W0 exposed for the mean-parameterisation R helpers (theta = W0(mu)).
// [[Rcpp::export]]
Rcpp::NumericVector bell_W0_cpp(Rcpp::NumericVector x)
{
  int n = x.size();
  Rcpp::NumericVector out(n);
  for (int i = 0; i < n; i++) out[i] = lambertW0(x[i]);
  return out;
}
