// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

// Poisson Ramos-Louzada (PRL), Alkhairy (2023), Math. Biosci. Eng. 20(8):14061-14080.
// A one-parameter mixed Poisson: X | theta ~ Poisson(theta) with theta ~ RL(tau).
//
// pmf (their Eq. 2), for x = 0,1,2,... and tau >= 2:
//
//     P(X = x; tau) = (1 + 1/tau)^-x * (x - 1 + tau(tau-1)) / ((tau-1)(1+tau)^2)
//
// Verified to sum to 1 (to 1e-6, tau in [2, 200]). Its moments, derived and then checked
// against the exact sums, are
//
//     E[X]   = tau^2 / (tau - 1)
//     Var[X] = tau^2 (tau^2 + tau - 3) / (tau - 1)^2
//     DI     = Var/E = (tau^2 + tau - 3) / (tau - 1)  ~  E[X]
//
// TWO STRUCTURAL CONSEQUENCES, both of which are the reason to fit it:
//
//  1. THE MEAN CANNOT GO BELOW 4. d/dtau [tau^2/(tau-1)] = tau(tau-2)/(tau-1)^2 vanishes at
//     tau = 2, so the mean is minimised there at mu = 4 and increases monotonically after.
//     A GAM whose fitted mean dips below 4 in low-count weeks -- routine for the rarer
//     diseases -- cannot be represented at all. The clamp records this rather than hiding it.
//
//  2. DISPERSION IS PINNED TO THE MEAN, DI ~ mu. Unlike Bell (DI ~ 1) this is geometric-like
//     rather than near-equidispersed, so PRL is a far stronger one-parameter model than Bell.
//     But it still has no free dispersion: it must be overdispersed by exactly the amount its
//     mean dictates, at every observation.
//
// PARAMETERISATION. eta = log mu, matching every other count family here so that the AIC
// comparison holds the mean structure fixed. Inverting mu = tau^2/(tau-1) gives
// tau^2 - mu tau + mu = 0, hence tau = (mu + sqrt(mu^2 - 4 mu)) / 2, the root with tau >= 2.
// dtau/dmu diverges as mu -> 4 (because dmu/dtau = 0 there), so mu is clamped a little above
// 4 to keep the score finite.
//
// The score is analytic; the Hessian is the BHHH approximation, as for the other families.
// d12 packs, per observation, the NEGATIVE log-likelihood derivatives: cols [g1, H11].

static const double PRL_MU_LO = 4.0 + 1e-4;      // dtau/dmu -> Inf as mu -> 4
static const double PRL_MU_HI = 1e12;

static inline double prl_clamp(double x, double lo, double hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

// tau >= 2 from the mean, and dtau/dmu.
static inline void prl_tau(double mu, double& tau, double& dtau_dmu) {
  double disc = std::sqrt(mu * mu - 4.0 * mu);
  tau = 0.5 * (mu + disc);
  dtau_dmu = 0.5 * (1.0 + (mu - 2.0) / disc);
}

// log pmf and d log pmf / d eta (eta = log mu).
static void prl_lpg(double y, double e1, double& lp, double* g)
{
  double mu = prl_clamp(std::exp(e1), PRL_MU_LO, PRL_MU_HI);
  double tau, dtau_dmu;
  prl_tau(mu, tau, dtau_dmu);

  double A = y - 1.0 + tau * (tau - 1.0);      // > 0 for tau >= 2 (worst case y=0: tau^2-tau-1)
  lp = -y * std::log1p(1.0 / tau) + std::log(A)
     - std::log(tau - 1.0) - 2.0 * std::log1p(tau);
  if (!g) return;

  double dl_dtau = y / (tau * (tau + 1.0))            // -y d/dtau log((tau+1)/tau)
                 + (2.0 * tau - 1.0) / A
                 - 1.0 / (tau - 1.0)
                 - 2.0 / (tau + 1.0);
  g[0] = dl_dtau * dtau_dmu * mu;                     // chain to eta = log mu
}

static inline void prl_pack(arma::mat& out, int j, double g) {
  out(j, 0) = g;
  out(j, 1) = g * g;                                  // BHHH
}

static void prl_eta(const Rcpp::List& pars, const arma::mat& X1, const arma::uvec& dupid,
                    int dcate, const Rcpp::List& offsets, arma::vec& e1)
{
  e1 = X1 * Rcpp::as<arma::vec>(pars[0]);
  if (dcate == 1) e1 = e1.elem(dupid);
  arma::vec o0 = Rcpp::as<arma::vec>(offsets[0]);
  if (o0.n_elem > 0) e1 += o0;
}

// [[Rcpp::export]]
double prl1d0(const Rcpp::List& pars, const arma::mat& X1, arma::vec yvec,
              const arma::uvec& dupid, int dcate, const Rcpp::List& offsets) {
  arma::vec e1;
  prl_eta(pars, X1, dupid, dcate, offsets, e1);
  double nllh = 0.0, lp;
  for (arma::uword j = 0; j < yvec.n_elem; j++) {
    prl_lpg(yvec[j], e1[j], lp, NULL);
    if (!std::isfinite(lp)) return std::numeric_limits<double>::infinity();
    nllh -= lp;
  }
  return nllh;
}

// [[Rcpp::export]]
arma::mat prl1d12(const Rcpp::List& pars, arma::mat X1, arma::vec yvec,
                  const arma::uvec dupid, int dcate, const Rcpp::List& offsets) {
  arma::vec e1;
  prl_eta(pars, X1, dupid, dcate, offsets, e1);
  arma::mat out(yvec.n_elem, 2);
  double lp, g[1];
  for (arma::uword j = 0; j < yvec.n_elem; j++) {
    prl_lpg(yvec[j], e1[j], lp, g);
    prl_pack(out, j, -g[0]);                          // NLL gradient
  }
  return out;
}

// ------------------------------------------------------------------ test hooks ----

// [[Rcpp::export]]
Rcpp::NumericVector prl_logpmf_cpp(Rcpp::NumericVector y, double tau) {
  Rcpp::NumericVector out(y.size());
  for (int i = 0; i < y.size(); i++) {
    double A = y[i] - 1.0 + tau * (tau - 1.0);
    out[i] = (A > 0.0)
      ? -y[i] * std::log1p(1.0 / tau) + std::log(A) - std::log(tau - 1.0) - 2.0 * std::log1p(tau)
      : -std::numeric_limits<double>::infinity();
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector prl_grad_cpp(double y, double eta) {
  double lp, g[1];
  prl_lpg(y, eta, lp, g);
  return Rcpp::NumericVector::create(Rcpp::Named("logp") = lp, Rcpp::Named("g") = g[0]);
}

// The clamp, exported so that predict.egpd() and the tests read the same number the
// likelihood enforces rather than duplicating it (cf. cf_bounds_cpp in countfams.cpp).
// [[Rcpp::export]]
Rcpp::NumericVector prl_bounds_cpp() {
  return Rcpp::NumericVector::create(Rcpp::Named("mu_lo") = PRL_MU_LO,
                                     Rcpp::Named("mu_hi") = PRL_MU_HI);
}
