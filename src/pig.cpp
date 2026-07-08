// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

// Poisson-inverse Gaussian (PIG) and zero-inflated PIG (ZIPIG) likelihoods for
// the native GAM fitter (egpd()).
//
// Parameterisation follows gamlss.dist::PIG:
//   mu > 0    is the mean,
//   sigma > 0 is the dispersion (Var = mu + sigma * mu^2),
//   pi in (0,1) (ZIPIG only) is the extra zero-inflation probability.
// Links: log(mu) (X1), log(sigma) (X2), logit(pi) (X3, ZIPIG only).
//
// log f(y) = -lgamma(y+1) + (1 - b)/sigma + sum_{k=0}^{y-1} log t_k,
//   b = sqrt(1 + 2*sigma*mu),
// where t_y = mu * E[Z | Y = y] is computed by the forward recursion of
// gamlss.dist's tofyPIG1 (M. Enea, 2014):
//   t_0 = mu / b,   t_j = (sigma*(2j-1)/mu + 1/t_{j-1}) * t_0^2.
//
// Both the score AND the exact observed Hessian are analytic. The Hessian uses
// the second posterior moment mu^2 * E[Z^2|y] = t_y * t_{y+1} (one extra
// recursion step), giving, on the (log mu, log sigma) scale,
//   d2l/deta_mu^2      = t_y*t_{y+1} - t_y^2 - t_y,
//   d2l/deta_mu deta_s = -sigma * dty_dsigma,
//   d2l/deta_sigma^2   = dty_dsigma*(1+sigma*mu)/mu + t_y - y - s2,
// with dty_dsigma = t_y * (Ssig_{y+1} - Ssig_y). Verified against numerical
// differentiation to ~1e-10. The observed Hessian is not guaranteed
// positive-definite; the fitter perturbs it when required (as for the EGPD
// families).
//
// d12 returns per-observation derivatives of the NEGATIVE log-likelihood with
// respect to the linear predictors, packed to match .gH2 / .gH3:
//   PIG (npar=2):  cols [g1, g2, H11, H12, H22]                       (nobs x 5)
//   ZIPIG (npar=3): cols [g1, g2, g3, H11, H12, H13, H22, H23, H33]   (nobs x 9)

// Forward recursion for t_y; also accumulates sumlt = sum_{k=0}^{y-1} log t_k
// and returns b via the out-parameter. (Used by d0.)
static inline double pig_ty(double y, double mu, double sigma,
                            double& sumlt, double& b) {
  b = std::sqrt(1.0 + 2.0 * sigma * mu);
  double t0 = mu / b;
  double c0sq = t0 * t0;
  double t = t0;
  sumlt = 0.0;
  int iy = (int) (y + 0.5);
  for (int j = 1; j <= iy; j++) {
    sumlt += std::log(t);                                // log t_{j-1}
    t = (sigma * (2.0 * j - 1.0) / mu + 1.0 / t) * c0sq; // t_j
  }
  return t;                                              // t_y
}

// Exact PIG log-likelihood score (s1, s2) and observed Hessian (H11, H12, H22)
// on the (log mu, log sigma) scale at count y. Also returns b and p0 = P(Y=0)
// (needed by the ZIPIG zero mixture). Runs the recursion to t_{y+1}.
static inline void pig_derivs(double y, double mu, double sigma,
                              double& s1, double& s2,
                              double& H11, double& H12, double& H22,
                              double& b, double& p0) {
  b = std::sqrt(1.0 + 2.0 * sigma * mu);
  double t0 = mu / b;
  double c0sq = t0 * t0;
  double t = t0;
  int iy = (int) (y + 0.5);
  double ty = t0, typ1 = t0;
  for (int j = 1; j <= iy + 1; j++) {
    t = (sigma * (2.0 * j - 1.0) / mu + 1.0 / t) * c0sq;  // t_j
    if (j == iy)     ty = t;
    if (j == iy + 1) typ1 = t;
  }
  if (iy == 0) ty = t0;

  double s2i = sigma * sigma;
  double ty2 = ty * typ1;                                 // mu^2 E[Z^2|y]
  double Ssig_y   = (ty   * (1.0 + sigma * mu) / mu - (1.0 + sigma * y))         / s2i;
  double Ssig_yp1 = (typ1 * (1.0 + sigma * mu) / mu - (1.0 + sigma * (y + 1.0))) / s2i;
  double dty_ds = ty * (Ssig_yp1 - Ssig_y);

  s1 = y - ty;
  s2 = sigma * Ssig_y;
  H11 = ty2 - ty * ty - ty;
  H12 = -sigma * dty_ds;
  H22 = dty_ds * (1.0 + sigma * mu) / mu + ty - y - s2;
  p0 = std::exp((1.0 - b) / sigma);
}

// ---------------------------------------------------------------- PIG ---------

// [[Rcpp::export]]
double pig1d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2,
              arma::vec yvec, const arma::uvec& dupid, int dcate,
              const Rcpp::List& offsets)
{
  arma::vec lmuvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lsigmavec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = yvec.size();

  if (dcate == 1) {
    lmuvec = lmuvec.elem(dupid);
    lsigmavec = lsigmavec.elem(dupid);
  }
  {
    arma::vec off0 = Rcpp::as<arma::vec>(offsets[0]);
    if (off0.n_elem > 0) lmuvec += off0;
    arma::vec off1 = Rcpp::as<arma::vec>(offsets[1]);
    if (off1.n_elem > 0) lsigmavec += off1;
  }

  double nllh = 0.0;
  for (int j = 0; j < nobs; j++) {
    double y = yvec[j];
    double mu = std::exp(lmuvec[j]);
    double sigma = std::exp(lsigmavec[j]);
    double sumlt, b;
    pig_ty(y, mu, sigma, sumlt, b);
    double logf = -R::lgammafn(y + 1.0) + (1.0 - b) / sigma + sumlt;
    nllh += -logf;
  }
  return nllh;
}

// [[Rcpp::export]]
arma::mat pig1d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2,
                  arma::vec yvec, const arma::uvec dupid, int dcate,
                  const Rcpp::List& offsets)
{
  arma::vec lmuvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lsigmavec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = yvec.size();
  arma::mat out = arma::mat(nobs, 5);

  if (dcate == 1) {
    lmuvec = lmuvec.elem(dupid);
    lsigmavec = lsigmavec.elem(dupid);
  }
  {
    arma::vec off0 = Rcpp::as<arma::vec>(offsets[0]);
    if (off0.n_elem > 0) lmuvec += off0;
    arma::vec off1 = Rcpp::as<arma::vec>(offsets[1]);
    if (off1.n_elem > 0) lsigmavec += off1;
  }

  for (int j = 0; j < nobs; j++) {
    double y = yvec[j];
    double mu = std::exp(lmuvec[j]);
    double sigma = std::exp(lsigmavec[j]);
    double s1, s2, H11, H12, H22, b, p0;
    pig_derivs(y, mu, sigma, s1, s2, H11, H12, H22, b, p0);
    // NLL = negative log-likelihood: negate score and Hessian.
    out(j, 0) = -s1;
    out(j, 1) = -s2;
    out(j, 2) = -H11;
    out(j, 3) = -H12;
    out(j, 4) = -H22;
  }
  return out;
}

// -------------------------------------------------------------- ZIPIG ---------

// [[Rcpp::export]]
double zipig1d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2,
                const arma::mat& X3, arma::vec yvec, const arma::uvec& dupid,
                int dcate, const Rcpp::List& offsets)
{
  arma::vec lmuvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lsigmavec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec logitpivec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = yvec.size();

  if (dcate == 1) {
    lmuvec = lmuvec.elem(dupid);
    lsigmavec = lsigmavec.elem(dupid);
    logitpivec = logitpivec.elem(dupid);
  }
  {
    arma::vec off0 = Rcpp::as<arma::vec>(offsets[0]);
    if (off0.n_elem > 0) lmuvec += off0;
    arma::vec off1 = Rcpp::as<arma::vec>(offsets[1]);
    if (off1.n_elem > 0) lsigmavec += off1;
    arma::vec off2 = Rcpp::as<arma::vec>(offsets[2]);
    if (off2.n_elem > 0) logitpivec += off2;
  }

  double nllh = 0.0;
  for (int j = 0; j < nobs; j++) {
    double y = yvec[j];
    double mu = std::exp(lmuvec[j]);
    double sigma = std::exp(lsigmavec[j]);
    double pi = 1.0 / (1.0 + std::exp(-logitpivec[j]));
    double sumlt, b;
    pig_ty(y, mu, sigma, sumlt, b);
    double logf;
    if (y < 0.5) {
      double p0 = std::exp((1.0 - b) / sigma);
      logf = std::log(pi + (1.0 - pi) * p0);
    } else {
      double logf_pig = -R::lgammafn(y + 1.0) + (1.0 - b) / sigma + sumlt;
      logf = std::log1p(-pi) + logf_pig;
    }
    nllh += -logf;
  }
  return nllh;
}

// [[Rcpp::export]]
arma::mat zipig1d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2,
                    arma::mat X3, arma::vec yvec, const arma::uvec dupid,
                    int dcate, const Rcpp::List& offsets)
{
  arma::vec lmuvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lsigmavec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec logitpivec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = yvec.size();
  arma::mat out = arma::mat(nobs, 9);

  if (dcate == 1) {
    lmuvec = lmuvec.elem(dupid);
    lsigmavec = lsigmavec.elem(dupid);
    logitpivec = logitpivec.elem(dupid);
  }
  {
    arma::vec off0 = Rcpp::as<arma::vec>(offsets[0]);
    if (off0.n_elem > 0) lmuvec += off0;
    arma::vec off1 = Rcpp::as<arma::vec>(offsets[1]);
    if (off1.n_elem > 0) lsigmavec += off1;
    arma::vec off2 = Rcpp::as<arma::vec>(offsets[2]);
    if (off2.n_elem > 0) logitpivec += off2;
  }

  for (int j = 0; j < nobs; j++) {
    double y = yvec[j];
    double mu = std::exp(lmuvec[j]);
    double sigma = std::exp(lsigmavec[j]);
    double pi = 1.0 / (1.0 + std::exp(-logitpivec[j]));
    double pip = pi * (1.0 - pi);            // dpi/deta_pi
    double pipp = pip * (1.0 - 2.0 * pi);    // d2pi/deta_pi^2

    double u1, u2, Hp11, Hp12, Hp22, b, p0;
    pig_derivs(y, mu, sigma, u1, u2, Hp11, Hp12, Hp22, b, p0);

    // log-likelihood score (se) and Hessian (He, upper triangle) on
    // (eta_mu, eta_sigma, eta_pi).
    double se1, se2, se3;
    double He11, He12, He13, He22, He23, He33;

    if (y > 0.5) {
      se1 = u1; se2 = u2; se3 = -pi;
      He11 = Hp11; He12 = Hp12; He22 = Hp22;
      He13 = 0.0; He23 = 0.0; He33 = -pip;
    } else {
      double q = (1.0 - pi) * p0;
      double D = pi + q;
      double D1 = q * u1, D2 = q * u2, D3 = pip * (1.0 - p0);
      double D11 = q * (u1 * u1 + Hp11);
      double D12 = q * (u1 * u2 + Hp12);
      double D22 = q * (u2 * u2 + Hp22);
      double D13 = -pip * p0 * u1;
      double D23 = -pip * p0 * u2;
      double D33 = pipp * (1.0 - p0);
      se1 = D1 / D; se2 = D2 / D; se3 = D3 / D;
      He11 = D11 / D - se1 * se1;
      He12 = D12 / D - se1 * se2;
      He13 = D13 / D - se1 * se3;
      He22 = D22 / D - se2 * se2;
      He23 = D23 / D - se2 * se3;
      He33 = D33 / D - se3 * se3;
    }

    // NLL = negative log-likelihood: negate score and Hessian.
    out(j, 0) = -se1;
    out(j, 1) = -se2;
    out(j, 2) = -se3;
    out(j, 3) = -He11;
    out(j, 4) = -He12;
    out(j, 5) = -He13;
    out(j, 6) = -He22;
    out(j, 7) = -He23;
    out(j, 8) = -He33;
  }
  return out;
}
