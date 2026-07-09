// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

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
// information Hessian are exact (not the BHHH approximation used for PIG/GPIG).
//
// d12 returns per-observation derivatives of the NEGATIVE log-likelihood wrt
// the linear predictors, packed to match .gH1 / .gH2 (upper triangle,
// row-major):
//   npar=1: cols [g1, H11]                           (nobs x 2)
//   npar=2: cols [g1, g2, H11, H12, H22]             (nobs x 5)

// ------------------------------------------------------ Bell ratios / W0 ------

// Fill la[0..nmax] with log(B_k / k!) via the stable log-space recursion
//   a_n = (1/n) sum_{m=0}^{n-1} a_m / (n-1-m)!,   a_0 = 1,
// i.e. la_n = -log(n) + logsumexp_{m<n} ( la_m - lgamma(n-m) ).
static inline void log_bell_ratio(int nmax, std::vector<double>& la)
{
  la.assign(nmax + 1, 0.0);
  la[0] = 0.0;                                   // a_0 = B_0/0! = 1
  for (int n = 1; n <= nmax; n++) {
    double mx = -std::numeric_limits<double>::infinity();
    for (int m = 0; m < n; m++) {
      double t = la[m] - std::lgamma((double)(n - m));   // (n-1-m)! = Gamma(n-m)
      if (t > mx) mx = t;
    }
    double s = 0.0;
    for (int m = 0; m < n; m++)
      s += std::exp(la[m] - std::lgamma((double)(n - m)) - mx);
    la[n] = -std::log((double) n) + mx + std::log(s);
  }
}

// Principal branch of the Lambert W function for x >= 0 (Halley iteration).
static inline double lambertW0(double x)
{
  if (x <= 0.0) return 0.0;                       // W0(0) = 0
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

// eta -> theta with dtheta/deta and d2theta/deta2.
//   native (param 0): theta = exp(eta),          dtheta = theta, d2theta = theta
//   mean   (param 1): mu = exp(eta), theta = W0(mu),
//                     dtheta = theta/(1+theta),  d2theta = theta/(1+theta)^3
static inline void eta_to_theta(int param, double eta,
                                double& theta, double& dth, double& d2th)
{
  if (param == 0) {
    theta = std::exp(eta);
    dth = theta; d2th = theta;
  } else {
    double mu = std::exp(eta);
    theta = lambertW0(mu);
    double opt = 1.0 + theta;
    dth  = theta / opt;
    d2th = theta / (opt * opt * opt);
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
  std::vector<double> la; log_bell_ratio(ymax, la);

  double nllh = 0.0;
  for (int j = 0; j < nobs; j++) {
    double theta, dth, d2th; eta_to_theta(param, e1[j], theta, dth, d2th);
    int y = (int)(yvec[j] + 0.5);
    double logf = y * std::log(theta) - (std::exp(theta) - 1.0) + la[y];
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
    double theta, dth, d2th; eta_to_theta(param, e1[j], theta, dth, d2th);
    int y = (int)(yvec[j] + 0.5);
    double eth = std::exp(theta);
    double g_th = eth - y / theta;                 // dNLL/dtheta
    double h_th = eth + y / (theta * theta);       // d2NLL/dtheta2
    out(j, 0) = g_th * dth;                                 // g1
    out(j, 1) = h_th * dth * dth + g_th * d2th;             // H11
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
  std::vector<double> la; log_bell_ratio(ymax, la);

  double nllh = 0.0;
  for (int j = 0; j < nobs; j++) {
    double theta, dth, d2th; eta_to_theta(param, e1[j], theta, dth, d2th);
    double pi = 1.0 / (1.0 + std::exp(-e2[j]));
    int y = (int)(yvec[j] + 0.5);
    double logf;
    if (y == 0) {
      double p0 = std::exp(1.0 - std::exp(theta));
      logf = std::log(pi + (1.0 - pi) * p0);
    } else {
      logf = std::log1p(-pi) + y * std::log(theta) - (std::exp(theta) - 1.0) + la[y];
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
    double theta, dth, d2th; eta_to_theta(param, e1[j], theta, dth, d2th);
    double pi = 1.0 / (1.0 + std::exp(-e2[j]));
    double dpi = pi * (1.0 - pi);                  // dpi/deta_pi
    double d2pi = dpi * (1.0 - 2.0 * pi);          // d2pi/deta_pi2
    int y = (int)(yvec[j] + 0.5);
    double eth = std::exp(theta);

    double NLL_th, NLL_pi, NLL_thth, NLL_thpi, NLL_pipi;
    if (y > 0) {
      NLL_th   = eth - y / theta;
      NLL_thth = eth + y / (theta * theta);
      NLL_pi   = 1.0 / (1.0 - pi);
      NLL_pipi = 1.0 / ((1.0 - pi) * (1.0 - pi));
      NLL_thpi = 0.0;
    } else {
      double p0 = std::exp(1.0 - eth);
      double p0_th   = -eth * p0;                  // dp0/dtheta
      double p0_thth = p0 * eth * (eth - 1.0);     // d2p0/dtheta2
      double D    = pi + (1.0 - pi) * p0;
      double D_th   = (1.0 - pi) * p0_th;
      double D_pi   = 1.0 - p0;
      double D_thth = (1.0 - pi) * p0_thth;
      double D_thpi = -p0_th;                      // d2D/dtheta dpi
      NLL_th   = -D_th / D;
      NLL_pi   = -D_pi / D;
      NLL_thth = -D_thth / D + (D_th * D_th) / (D * D);
      NLL_thpi = -D_thpi / D + (D_th * D_pi) / (D * D);
      NLL_pipi = (D_pi * D_pi) / (D * D);          // D_pipi = 0
    }
    out(j, 0) = NLL_th * dth;                                   // g1
    out(j, 1) = NLL_pi * dpi;                                   // g2
    out(j, 2) = NLL_thth * dth * dth + NLL_th * d2th;           // H11
    out(j, 3) = NLL_thpi * dth * dpi;                           // H12
    out(j, 4) = NLL_pipi * dpi * dpi + NLL_pi * d2pi;           // H22
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
  Rcpp::NumericVector out(n + 1);
  std::vector<double> la; log_bell_ratio(n, la);
  double base = 1.0 - std::exp(theta);            // log Pr(Y=0) = 1 - e^theta
  double lth = std::log(theta);
  for (int k = 0; k <= n; k++)
    out[k] = std::exp(k * lth + base + la[k]);
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
