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
// Both the score AND exact observed Hessian are analytic, obtained by forward
// differentiation of the same stable recursion on the log-link scale. This
// avoids subtracting posterior moments and dividing by sigma^2 near Poisson.
// The observed Hessian is not guaranteed
// positive-definite; the fitter perturbs it when required (as for the EGPD
// families).
//
// d12 returns per-observation derivatives of the NEGATIVE log-likelihood with
// respect to the linear predictors, packed to match .gH2 / .gH3:
//   PIG (npar=2):  cols [g1, g2, H11, H12, H22]                       (nobs x 5)
//   ZIPIG (npar=3): cols [g1, g2, g3, H11, H12, H13, H22, H23, H33]   (nobs x 9)

// Forward recursion for t_y; also accumulates sumlt = sum_{k=0}^{y-1} log t_k
// and returns b via the out-parameter. (Used by d0.)
struct PigJet {
  double v, g[2], h[2][2];
  PigJet(double x) : v(x), g{0, 0}, h{{0, 0}, {0, 0}} {}
};
static PigJet operator+(const PigJet& a, const PigJet& b) {
  PigJet z(a.v + b.v);
  for (int i=0;i<2;i++) { z.g[i]=a.g[i]+b.g[i];
    for (int j=0;j<2;j++) z.h[i][j]=a.h[i][j]+b.h[i][j]; }
  return z;
}
static PigJet operator*(const PigJet& a, const PigJet& b) {
  PigJet z(a.v*b.v);
  for (int i=0;i<2;i++) { z.g[i]=a.g[i]*b.v+a.v*b.g[i];
    for (int j=0;j<2;j++) z.h[i][j]=a.h[i][j]*b.v+a.v*b.h[i][j]+a.g[i]*b.g[j]+a.g[j]*b.g[i]; }
  return z;
}
static PigJet pig_transform(const PigJet& a, double v, double d, double d2) {
  PigJet z(v);
  for (int i=0;i<2;i++) { z.g[i]=d*a.g[i];
    for (int j=0;j<2;j++) z.h[i][j]=d*a.h[i][j]+d2*a.g[i]*a.g[j]; }
  return z;
}
static PigJet pig_inv(const PigJet& a) {
  return pig_transform(a,1/a.v,-1/(a.v*a.v),2/(a.v*a.v*a.v));
}
static PigJet pig_log(const PigJet& a) {
  return pig_transform(a,std::log(a.v),1/a.v,-1/(a.v*a.v));
}

// The recursion as written squares t_0 = mu/b and differentiates through 1/t^3, and forms
// 2 sigma mu directly: t_0^2 is subnormal for mu below ~1e-154 and zero below ~2e-162 (NaN
// pmf), 2/t_0^3 overflows below ~2e-103 (NaN Hessian), and b overflows once sigma mu
// exceeds ~9e307 (p_0 = 1). It is kept where t_0 lies in [1e-100, 1e100] and 2 sigma mu
// <= 1e300; elsewhere the same recursion is carried scaled by t_0 and in logs:
//   u_j = t_j / t_0 = (2j-1) sigma/b + 1/u_{j-1},  u_0 = 1,
// so sum_{k<y} log t_k = y log t_0 + sum_{k<y} log u_k, with log b = log1pexp(log 2 sigma mu)/2
// and log t_0 = log mu - log b. Every u_j >= 1 (it is E[Z | Y = j] / E[Z | Y = 0]), so
// nothing under- or overflows.
static inline bool pig_direct(double mu, double sigma) {
  double q = 2.0 * sigma * mu;
  if (!(q <= 1e300)) return false;
  double t0 = mu / std::sqrt(1.0 + q);
  return t0 >= 1e-100 && t0 <= 1e100;
}

// log b, b, log t_0, log s with s = sigma/b, and log p_0 = -2 mu/(1+b), all formed in logs.
static inline void pig_scaled_init(double mu, double sigma, double& lb, double& b,
                                   double& lt0, double& lsb, double& lp0) {
  double lm = std::log(mu), lsg = std::log(sigma);
  lb = 0.5 * R::log1pexp(M_LN2 + lsg + lm);
  b = std::exp(lb);
  lt0 = lm - lb;
  lsb = lsg - lb;
  lp0 = -2.0 * std::exp(lm - R::log1pexp(lb));
}

// log u_j from log u_{j-1}: log((2j-1) s + exp(-log u_{j-1})).
static inline double pig_lu_next(int j, double lsb, double lu) {
  double A = std::log(2.0 * j - 1.0) + lsb, B = -lu;
  return (A > B) ? A + std::log1p(std::exp(B - A)) : B + std::log1p(std::exp(A - B));
}

static inline double pig_ty(double y, double mu, double sigma,
                            double& sumlt, double& b, double& lp0) {
  int iy = (int) (y + 0.5);
  if (!pig_direct(mu, sigma)) {
    double lb, lt0, lsb, lu = 0.0, su = 0.0;
    pig_scaled_init(mu, sigma, lb, b, lt0, lsb, lp0);
    for (int j = 1; j < iy; j++) {
      lu = pig_lu_next(j, lsb, lu);                      // log u_j
      su += lu;
    }
    sumlt = iy * lt0 + su;
    return std::exp(lt0 + (iy > 0 ? pig_lu_next(iy, lsb, lu) : 0.0));   // t_y
  }
  b = std::sqrt(1.0 + 2.0 * sigma * mu);
  lp0 = -2.0 * mu / (1.0 + b);
  double t0 = mu / b;
  double c0sq = t0 * t0;
  double t = t0;
  sumlt = 0.0;
  for (int j = 1; j <= iy; j++) {
    sumlt += std::log(t);                                // log t_{j-1}
    t = (sigma * (2.0 * j - 1.0) / mu + 1.0 / t) * c0sq; // t_j
  }
  return t;                                              // t_y
}

static inline double pig_logsumexp2(double a, double b) {
  if (a == b) return a + M_LN2;
  double m = (a > b) ? a : b;
  if (!std::isfinite(m)) return m;
  return m + R::log1pexp(-std::fabs(a - b));
}

// The scaled recursion of pig_ty, differentiated on the (log mu, log sigma) scale. Each
// log u_j is carried with the gradient and Hessian of the LOG, which the recursion maps by
// the log-sum-exp rule; with weights al + be = 1 every term stays O(1), so no derivative of
// a tiny or huge t_j (or of 1/mu) is ever formed. Used where pig_direct() is false.
static void pig_derivs_scaled(double y, double mu, double sigma,
                              double& s1, double& s2,
                              double& H11, double& H12, double& H22,
                              double& b, double& lp0) {
  double lb, lt0, lsb;
  pig_scaled_init(mu, sigma, lb, b, lt0, lsb, lp0);
  double lq = M_LN2 + std::log(sigma) + std::log(mu);
  // log b = log1pexp(lq)/2: gradient gb (1, 1), Hessian hb J (J = all ones)
  double p = R::plogis(lq, 0.0, 1.0, 1, 0), pc = R::plogis(lq, 0.0, 1.0, 0, 0);
  double gb = 0.5 * p, hb = 0.5 * p * pc;
  // log(-lp0) = log 2 + log mu - log(1 + b), with r = b/(1+b)
  double r = R::plogis(lb, 0.0, 1.0, 1, 0), rc = R::plogis(lb, 0.0, 1.0, 0, 0);
  double w[2] = { 1.0 - r * gb, -r * gb };
  double hw = -(r * hb + r * rc * gb * gb);
  double G[2], H[2][2];
  for (int i = 0; i < 2; i++) {
    G[i] = lp0 * w[i];
    for (int k = 0; k < 2; k++) H[i][k] = lp0 * (hw + w[i] * w[k]);
  }
  // y log t_0, log t_0 = log mu - log b
  int iy = (int)(y + 0.5);
  G[0] += iy * (1.0 - gb); G[1] -= iy * gb;
  for (int i = 0; i < 2; i++) for (int k = 0; k < 2; k++) H[i][k] -= iy * hb;
  // log s = log sigma - log b; log u_j = log(e^A + e^B), A = log(2j-1) + log s, B = -log u_{j-1}
  double gs[2] = { -gb, 1.0 - gb }, hs = -hb;
  double lu = 0.0, gu[2] = { 0.0, 0.0 }, hu[2][2] = { { 0.0, 0.0 }, { 0.0, 0.0 } };
  for (int j = 1; j < iy; j++) {
    double A = std::log(2.0 * j - 1.0) + lsb, B = -lu;
    double L = (A > B) ? A + std::log1p(std::exp(B - A)) : B + std::log1p(std::exp(A - B));
    double al = std::exp(A - L), be = std::exp(B - L);
    double gn[2] = { al * gs[0] - be * gu[0], al * gs[1] - be * gu[1] };
    double hn[2][2];
    for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
        hn[i][k] = al * (hs + gs[i] * gs[k]) + be * (gu[i] * gu[k] - hu[i][k]) - gn[i] * gn[k];
    lu = L;
    for (int i = 0; i < 2; i++) {
      gu[i] = gn[i]; G[i] += gn[i];
      for (int k = 0; k < 2; k++) { hu[i][k] = hn[i][k]; H[i][k] += hn[i][k]; }
    }
  }
  s1 = G[0]; s2 = G[1];
  H11 = H[0][0]; H12 = H[0][1]; H22 = H[1][1];
}

// Exact PIG log-likelihood score (s1, s2) and observed Hessian (H11, H12, H22)
// on the (log mu, log sigma) scale at count y. Also returns b and lp0 = log P(Y=0)
// (needed by the ZIPIG zero mixture).
static inline void pig_derivs(double y, double mu, double sigma,
                              double& s1, double& s2,
                              double& H11, double& H12, double& H22,
                              double& b, double& lp0) {
  if (!pig_direct(mu, sigma)) {
    pig_derivs_scaled(y, mu, sigma, s1, s2, H11, H12, H22, b, lp0);
    return;
  }
  // Differentiate the cancellation-free recursion directly on the log-link
  // scale; posterior-moment differences divided by sigma^2 fail near Poisson.
  PigJet m(mu), s(sigma);
  m.g[0]=m.h[0][0]=mu; s.g[1]=s.h[1][1]=sigma;
  PigJet b2=PigJet(1)+PigJet(2)*s*m;
  b=std::sqrt(b2.v);
  PigJet bj=pig_transform(b2,b,0.5/b,-0.25/(b*b*b));
  PigJet lp=PigJet(-2)*m*pig_inv(PigJet(1)+bj);
  lp0=lp.v;
  PigJet t=m*pig_inv(bj), t02=t*t;
  for (int j=1;j<=(int)(y+0.5);j++) {
    lp=lp+pig_log(t);
    t=(s*PigJet(2.0*j-1)*pig_inv(m)+pig_inv(t))*t02;
  }
  s1=lp.g[0]; s2=lp.g[1];
  H11=lp.h[0][0]; H12=lp.h[0][1]; H22=lp.h[1][1];
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
    double sumlt, b, lp0;
    pig_ty(y, mu, sigma, sumlt, b, lp0);
    double logf = -R::lgammafn(y + 1.0) + lp0 + sumlt;
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
    double eta_pi = logitpivec[j];
    double log_pi = -R::log1pexp(-eta_pi);
    double log1m_pi = -R::log1pexp(eta_pi);
    double sumlt, b, lp0;
    pig_ty(y, mu, sigma, sumlt, b, lp0);
    double logf;
    if (y < 0.5) {
      logf = pig_logsumexp2(log_pi, log1m_pi + lp0);
    } else {
      double logf_pig = -R::lgammafn(y + 1.0) + lp0 + sumlt;
      logf = log1m_pi + logf_pig;
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
    double eta_pi = logitpivec[j];
    double log_pi = -R::log1pexp(-eta_pi);
    double log1m_pi = -R::log1pexp(eta_pi);
    double pi = std::exp(log_pi);
    double one_m_pi = std::exp(log1m_pi);
    double pip = std::exp(log_pi + log1m_pi);    // dpi/deta_pi

    double u1, u2, Hp11, Hp12, Hp22, b, lp0;
    pig_derivs(y, mu, sigma, u1, u2, Hp11, Hp12, Hp22, b, lp0);

    // log-likelihood score (se) and Hessian (He, upper triangle) on
    // (eta_mu, eta_sigma, eta_pi).
    double se1, se2, se3;
    double He11, He12, He13, He22, He23, He33;

    if (y > 0.5) {
      se1 = u1; se2 = u2; se3 = -pi;
      He11 = Hp11; He12 = Hp12; He22 = Hp22;
      He13 = 0.0; He23 = 0.0; He33 = -pip;
    } else {
      double delta = eta_pi - lp0;
      double log_w0 = -R::log1pexp(delta);
      double log_w1 = -R::log1pexp(-delta);
      double w0 = (delta >= 0.0) ? std::exp(log_w0) : (1.0 - std::exp(log_w1));
      double w1 = (delta >= 0.0) ? (1.0 - w0) : std::exp(log_w1);
      double w01 = std::exp(log_w0 + log_w1);
      double w01_u1 = (w01 == 0.0) ? 0.0 : w01 * u1;
      double w01_u2 = (w01 == 0.0) ? 0.0 : w01 * u2;
      se1 = (w0 == 0.0) ? 0.0 : w0 * u1;
      se2 = (w0 == 0.0) ? 0.0 : w0 * u2;
      se3 = (eta_pi >= 0.0) ? (one_m_pi - w0) : (w1 - pi);
      He11 = ((w0 == 0.0) ? 0.0 : w0 * Hp11) + w01_u1 * u1;
      He12 = ((w0 == 0.0) ? 0.0 : w0 * Hp12) + w01_u1 * u2;
      He13 = -w01_u1;
      He22 = ((w0 == 0.0) ? 0.0 : w0 * Hp22) + w01_u2 * u2;
      He23 = -w01_u2;
      He33 = w01 - pip;
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

// ------------------------------------------------------ pmf for d/p/q/r --------

// Log pmf ladder log p_0, ..., log p_N of PIG(mu, sigma) for the standalone d/p/q/r
// functions, by the same cancellation-free recursion as the likelihood:
//   log p_0 = -2 mu / (1 + b),   log p_y = log p_{y-1} + log t_{y-1} - log y.
// Every term is positive, so the ladder is exact to rounding at any sigma, including the
// Poisson limit sigma -> 0 (t_j -> mu) where differencing posterior moments fails.
// [[Rcpp::export]]
Rcpp::NumericVector pig_logpmf_ladder_cpp(int N, double mu, double sigma)
{
  if (N < 0) Rcpp::stop("PIG support must be non-negative");
  Rcpp::NumericVector out(N + 1);
  if (!pig_direct(mu, sigma)) {                          // scaled by t_0, see pig_direct()
    double lb, b, lt0, lsb, lp0, lu = 0.0, sum = 0.0, comp = 0.0;
    pig_scaled_init(mu, sigma, lb, b, lt0, lsb, lp0);
    out[0] = lp0;
    for (int j = 1; j <= N; j++) {
      if ((j & 65535) == 0) Rcpp::checkUserInterrupt();
      if (j > 1) {
        lu = pig_lu_next(j - 1, lsb, lu);                // adds log u_{j-1}
        double s2 = sum + lu;
        comp += (std::fabs(sum) >= std::fabs(lu)) ? (sum - s2) + lu : (lu - s2) + sum;
        sum = s2;
      }
      out[j] = lp0 + (j * lt0 + (sum + comp)) - std::lgamma(j + 1.0);
    }
    return out;
  }
  double b = std::sqrt(1.0 + 2.0 * sigma * mu);
  double t0 = mu / b, c0sq = t0 * t0, t = t0;
  double lp0 = -2.0 * mu / (1.0 + b);
  double sum = 0.0, comp = 0.0;                          // Neumaier sum of log t_k
  out[0] = lp0;
  for (int j = 1; j <= N; j++) {
    if ((j & 65535) == 0) Rcpp::checkUserInterrupt();
    double v = std::log(t), s2 = sum + v;               // adds log t_{j-1}
    comp += (std::fabs(sum) >= std::fabs(v)) ? (sum - s2) + v : (v - s2) + sum;
    sum = s2;
    out[j] = lp0 + (sum + comp) - std::lgamma(j + 1.0);
    t = (sigma * (2.0 * j - 1.0) / mu + 1.0 / t) * c0sq; // t_j
  }
  return out;
}
