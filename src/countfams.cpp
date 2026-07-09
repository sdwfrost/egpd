// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <vector>

// Three further count families for the native GAM fitter (egpd()), added so that the
// discrete-EVT families can be compared against the standard mechanisms for
// overdispersion and for heavy tails.
//
//   gpois    Generalized Poisson (Consul & Jain), mean parameterisation
//   gwaring  Generalized Waring (a, k, rho), mean parameterisation
//   plnorm   Poisson-lognormal, mean parameterisation
//
// TAILS. Only one of the three is heavy-tailed in the extreme-value sense:
//
//   gpois    log p(y) ~ y (log lambda + 1 - lambda)   geometric decay,   xi = 0
//   plnorm   lognormal mixing                          Gumbel domain,     xi = 0
//   gwaring  log p(y) ~ -(rho+1) log y                 REGULARLY VARYING, xi = 1/rho > 0
//
// So gwaring is the direct rival to DEGPD / GPIG (a genuine power-law tail), plnorm is
// the heavy-body / light-tail contrast, and gpois is the overdispersion-only control
// alongside NB and PIG.
//
// PARAMETERISATIONS AND LINKS. predict() applies exp() to a "log" prefix and plogis() to
// a "logit" prefix, so only those two links are used; rho > 1 is imposed by a clamp
// rather than by a shifted link.
//
//   gpois    eta = (log mu, logit lambda),   mu > 0, lambda in (0,1)
//            theta = mu (1 - lambda);  Var = mu / (1-lambda)^2
//            log p = log theta + (y-1) log(theta + lambda y) - (theta + lambda y) - log y!
//
//   gwaring  eta = (log mu, log k, log rho), mu > 0, k > 0, rho > 1
//            a = mu (rho - 1) / k, so E[Y] = a k / (rho - 1) = mu
//            log p = lgamma(a+rho) + lgamma(k+rho) - lgamma(rho) - lgamma(a) - lgamma(k)
//                    + lgamma(a+y) + lgamma(k+y) - lgamma(a+k+rho+y) - lgamma(y+1)
//            Mean needs rho > 1, variance needs rho > 2. Reduces to Waring at k = 1.
//
//   plnorm   eta = (log mu, log sigma),      mu = E[Y] > 0, sigma > 0
//            Y | Z ~ Poisson(e^Z),  Z ~ N(muz, sigma^2),  muz = log(mu) - sigma^2/2
//            p(y) = int Pois(y | e^z) N(z; muz, sigma^2) dz    -- no closed form
//
// The Poisson-lognormal integral uses ADAPTIVE Gauss-Hermite: nodes are centred on the
// mode of the integrand and scaled by its curvature. A FIXED Gauss-Hermite grid fails
// outright here -- the Poisson kernel is a sharp spike when mu ~ 5e3, exactly the regime
// of the largest Tycho series, and fixed nodes then miss the mass and return zero. (The
// same failure mode bit the GPIG normalising constant; see gpig.cpp.) Its score is exact,
// not a finite difference: d log p / d(muz, sigma) are posterior moments of z, which the
// same nodes already supply.
//
// The score for every family is analytic; the Hessian is the empirical-Fisher / BHHH
// approximation (outer product of the per-observation NLL gradient), as for gpig.
//
// d12 returns per-observation derivatives of the NEGATIVE log-likelihood wrt the linear
// predictors, packed to match .gH2 / .gH3:
//   npar=2: cols [g1,g2, H11,H12,H22]                    (nobs x 5)
//   npar=3: cols [g1,g2,g3, H11,H12,H13,H22,H23,H33]     (nobs x 9)

static inline double cf_clamp(double x, double lo, double hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

static inline void cf_pack(arma::mat& out, int j, const double* g, int npar) {
  for (int i = 0; i < npar; i++) out(j, i) = g[i];
  int col = npar;
  for (int i = 0; i < npar; i++)
    for (int k = i; k < npar; k++)
      out(j, col++) = g[i] * g[k];
}

// ------------------------------------------------- Generalized Poisson -------------
// lambda is confined to [0,1). Negative lambda truncates the support, which makes the pmf
// improper and the likelihood discontinuous. Every Tycho series is overdispersed
// (variance/mean from 1.8 to 5.9e4), so the underdispersed branch is not needed.
static const double GP_LAM_LO = 1e-10, GP_LAM_HI = 1.0 - 1e-8;
static const double GP_MU_LO  = 1e-10, GP_MU_HI  = 1e12;

static void gpois_lpg(double y, double e1, double e2, double& lp, double* g)
{
  double mu  = cf_clamp(std::exp(e1), GP_MU_LO, GP_MU_HI);
  double lam = cf_clamp(1.0 / (1.0 + std::exp(-e2)), GP_LAM_LO, GP_LAM_HI);
  double th  = mu * (1.0 - lam);
  double A   = th + lam * y;
  lp = std::log(th) + (y - 1.0) * std::log(A) - A - std::lgamma(y + 1.0);
  if (!g) return;
  double P = 1.0 / th + (y - 1.0) / A - 1.0;           // dl/dtheta at fixed lambda
  double dl_dmu  = P * (1.0 - lam);
  double dl_dlam = -mu * P + y * ((y - 1.0) / A - 1.0);
  g[0] = dl_dmu * mu;                                   // d/d eta1 (log mu)
  g[1] = dl_dlam * lam * (1.0 - lam);                   // d/d eta2 (logit lambda)
}

// -------------------------------------------------- Generalized Waring -------------
// rho > 1 keeps the mean finite; rho > 2 is needed for a finite variance, but that is a
// property of the fit, not a constraint to impose -- a fitted rho in (1,2] is a meaningful
// statement (infinite variance) and is reported rather than forbidden.
static const double GW_RHO_LO = 1.0 + 1e-6, GW_RHO_HI = 1e6;
static const double GW_K_LO   = 1e-8,  GW_K_HI  = 1e8;
static const double GW_MU_LO  = 1e-10, GW_MU_HI = 1e12;

static void gwaring_lpg(double y, double e1, double e2, double e3, double& lp, double* g)
{
  double mu  = cf_clamp(std::exp(e1), GW_MU_LO,  GW_MU_HI);
  double k   = cf_clamp(std::exp(e2), GW_K_LO,   GW_K_HI);
  double rho = cf_clamp(std::exp(e3), GW_RHO_LO, GW_RHO_HI);
  double a   = mu * (rho - 1.0) / k;
  if (!(a > 1e-12)) a = 1e-12;
  double akr = a + k + rho;

  lp = std::lgamma(a + rho) + std::lgamma(k + rho) - std::lgamma(rho)
     - std::lgamma(a) - std::lgamma(k)
     + std::lgamma(a + y) + std::lgamma(k + y) - std::lgamma(akr + y)
     - std::lgamma(y + 1.0);
  if (!g) return;

  double la = R::digamma(a + rho) - R::digamma(a) + R::digamma(a + y) - R::digamma(akr + y);
  double lk = R::digamma(k + rho) - R::digamma(k) + R::digamma(k + y) - R::digamma(akr + y);
  double lr = R::digamma(a + rho) + R::digamma(k + rho) - R::digamma(rho) - R::digamma(akr + y);

  // a = mu (rho-1)/k  =>  da/dmu = (rho-1)/k,  da/dk = -a/k,  da/drho = mu/k
  double dl_dmu  = la * (rho - 1.0) / k;
  double dl_dk   = lk - la * a / k;
  double dl_drho = lr + la * mu / k;

  g[0] = dl_dmu  * mu;                                  // d/d eta1 (log mu)
  g[1] = dl_dk   * k;                                   // d/d eta2 (log k)
  g[2] = dl_drho * rho;                                 // d/d eta3 (log rho)
}

// -------------------------------------------------- Poisson-lognormal -------------
// Gauss-Hermite nodes/weights by Golub-Welsch on the Hermite Jacobi matrix, cached.
static const int PLN_NQ = 21;

static void gh_nodes(int n, std::vector<double>& x, std::vector<double>& lw)
{
  static int cached_n = -1;
  static std::vector<double> cx, clw;
  if (n == cached_n) { x = cx; lw = clw; return; }
  arma::mat J(n, n, arma::fill::zeros);
  for (int i = 0; i < n - 1; i++) {
    double v = std::sqrt((i + 1) / 2.0);
    J(i, i + 1) = v; J(i + 1, i) = v;
  }
  arma::vec eval; arma::mat evec;
  arma::eig_sym(eval, evec, J);
  x.assign(n, 0.0); lw.assign(n, 0.0);
  for (int i = 0; i < n; i++) {
    x[i]  = eval(i);
    lw[i] = std::log(std::sqrt(M_PI) * evec(0, i) * evec(0, i));
  }
  cx = x; clw = lw; cached_n = n;
}

static const double PLN_MU_LO = 1e-10, PLN_MU_HI = 1e12;
static const double PLN_SG_LO = 1e-6,  PLN_SG_HI = 1e3;

static void plnorm_lpg(double y, double e1, double e2, double& lp, double* g)
{
  double mu = cf_clamp(std::exp(e1), PLN_MU_LO, PLN_MU_HI);
  double sg = cf_clamp(std::exp(e2), PLN_SG_LO, PLN_SG_HI);
  double s2 = sg * sg;
  double muz = std::log(mu) - 0.5 * s2;

  // h(z) = -e^z + y z - (z-muz)^2/(2 s2) is strictly concave; Newton from log(y+1/2).
  double z = std::log(y + 0.5);
  for (int it = 0; it < 100; it++) {
    double ez = std::exp(z);
    double h1 = -ez + y - (z - muz) / s2;
    double h2 = -ez - 1.0 / s2;
    double step = h1 / h2;
    z -= step;
    if (std::fabs(step) < 1e-12) break;
  }
  double ez = std::exp(z);
  double s  = 1.0 / std::sqrt(ez + 1.0 / s2);          // curvature scale at the mode

  std::vector<double> x, lw;
  gh_nodes(PLN_NQ, x, lw);

  std::vector<double> lt(PLN_NQ), zz(PLN_NQ);
  double mx = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < PLN_NQ; i++) {
    double zi = z + M_SQRT2 * s * x[i];
    zz[i] = zi;
    double d  = zi - muz;
    double hi = -std::exp(zi) + y * zi - d * d / (2.0 * s2);
    lt[i] = lw[i] + x[i] * x[i] + hi;                  // undo the e^{-x^2} of the rule
    if (lt[i] > mx) mx = lt[i];
  }
  double S = 0.0;
  for (int i = 0; i < PLN_NQ; i++) S += std::exp(lt[i] - mx);

  lp = -std::lgamma(y + 1.0) - 0.5 * std::log(2.0 * M_PI) - std::log(sg)
     + std::log(M_SQRT2 * s) + mx + std::log(S);
  if (!g) return;

  // The score is a posterior moment of z, available from the same nodes:
  //   dl/dmuz   = E[z - muz] / s2
  //   dl/dsigma = E[(z-muz)^2]/sigma^3 - 1/sigma        (at fixed muz)
  double m1 = 0.0, m2 = 0.0;
  for (int i = 0; i < PLN_NQ; i++) {
    double p = std::exp(lt[i] - mx) / S;
    double d = zz[i] - muz;
    m1 += p * d; m2 += p * d * d;
  }
  double s_muz = m1 / s2;
  double s_sg  = m2 / (s2 * sg) - 1.0 / sg;
  // muz = log(mu) - sigma^2/2, so d muz/d eta1 = 1 and d muz/d eta2 = -sigma^2.
  g[0] = s_muz;
  g[1] = s_sg * sg - s_muz * s2;
}

// ------------------------------------------------------- fitter plumbing ----------

static void cf_eta(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2,
                   const arma::mat* X3, const arma::uvec& dupid, int dcate,
                   const Rcpp::List& offsets,
                   arma::vec& e1, arma::vec& e2, arma::vec& e3)
{
  e1 = X1 * Rcpp::as<arma::vec>(pars[0]);
  e2 = X2 * Rcpp::as<arma::vec>(pars[1]);
  if (X3) e3 = (*X3) * Rcpp::as<arma::vec>(pars[2]);
  if (dcate == 1) {
    e1 = e1.elem(dupid); e2 = e2.elem(dupid);
    if (X3) e3 = e3.elem(dupid);
  }
  arma::vec o0 = Rcpp::as<arma::vec>(offsets[0]); if (o0.n_elem > 0) e1 += o0;
  arma::vec o1 = Rcpp::as<arma::vec>(offsets[1]); if (o1.n_elem > 0) e2 += o1;
  if (X3) { arma::vec o2 = Rcpp::as<arma::vec>(offsets[2]); if (o2.n_elem > 0) e3 += o2; }
}

// ---- gpois (2 par) ----
// [[Rcpp::export]]
double gpois1d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2,
                arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets) {
  arma::vec e1, e2, e3;
  cf_eta(pars, X1, X2, NULL, dupid, dcate, offsets, e1, e2, e3);
  double nllh = 0.0, lp;
  for (arma::uword j = 0; j < yvec.n_elem; j++) {
    gpois_lpg(yvec[j], e1[j], e2[j], lp, NULL);
    if (!std::isfinite(lp)) return std::numeric_limits<double>::infinity();
    nllh -= lp;
  }
  return nllh;
}
// [[Rcpp::export]]
arma::mat gpois1d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2,
                    arma::vec yvec, const arma::uvec dupid, int dcate, const Rcpp::List& offsets) {
  arma::vec e1, e2, e3;
  cf_eta(pars, X1, X2, NULL, dupid, dcate, offsets, e1, e2, e3);
  arma::mat out(yvec.n_elem, 5);
  double lp, gg[2];
  for (arma::uword j = 0; j < yvec.n_elem; j++) {
    gpois_lpg(yvec[j], e1[j], e2[j], lp, gg);
    double g[2] = { -gg[0], -gg[1] };                   // NLL gradient
    cf_pack(out, j, g, 2);
  }
  return out;
}

// ---- gwaring (3 par) ----
// [[Rcpp::export]]
double gwaring1d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2,
                  const arma::mat& X3, arma::vec yvec, const arma::uvec& dupid, int dcate,
                  const Rcpp::List& offsets) {
  arma::vec e1, e2, e3;
  cf_eta(pars, X1, X2, &X3, dupid, dcate, offsets, e1, e2, e3);
  double nllh = 0.0, lp;
  for (arma::uword j = 0; j < yvec.n_elem; j++) {
    gwaring_lpg(yvec[j], e1[j], e2[j], e3[j], lp, NULL);
    if (!std::isfinite(lp)) return std::numeric_limits<double>::infinity();
    nllh -= lp;
  }
  return nllh;
}
// [[Rcpp::export]]
arma::mat gwaring1d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3,
                      arma::vec yvec, const arma::uvec dupid, int dcate, const Rcpp::List& offsets) {
  arma::vec e1, e2, e3;
  cf_eta(pars, X1, X2, &X3, dupid, dcate, offsets, e1, e2, e3);
  arma::mat out(yvec.n_elem, 9);
  double lp, gg[3];
  for (arma::uword j = 0; j < yvec.n_elem; j++) {
    gwaring_lpg(yvec[j], e1[j], e2[j], e3[j], lp, gg);
    double g[3] = { -gg[0], -gg[1], -gg[2] };
    cf_pack(out, j, g, 3);
  }
  return out;
}

// ---- plnorm (2 par) ----
// [[Rcpp::export]]
double plnorm1d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2,
                 arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets) {
  arma::vec e1, e2, e3;
  cf_eta(pars, X1, X2, NULL, dupid, dcate, offsets, e1, e2, e3);
  double nllh = 0.0, lp;
  for (arma::uword j = 0; j < yvec.n_elem; j++) {
    plnorm_lpg(yvec[j], e1[j], e2[j], lp, NULL);
    if (!std::isfinite(lp)) return std::numeric_limits<double>::infinity();
    nllh -= lp;
  }
  return nllh;
}
// [[Rcpp::export]]
arma::mat plnorm1d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2,
                     arma::vec yvec, const arma::uvec dupid, int dcate, const Rcpp::List& offsets) {
  arma::vec e1, e2, e3;
  cf_eta(pars, X1, X2, NULL, dupid, dcate, offsets, e1, e2, e3);
  arma::mat out(yvec.n_elem, 5);
  double lp, gg[2];
  for (arma::uword j = 0; j < yvec.n_elem; j++) {
    plnorm_lpg(yvec[j], e1[j], e2[j], lp, gg);
    double g[2] = { -gg[0], -gg[1] };
    cf_pack(out, j, g, 2);
  }
  return out;
}

// ------------------------------------------------------------ test hooks ----------
// Natural-parameter log pmf, vectorised over y. Used by the d/p/q/r wrappers and by the
// test suite (sum-to-one, moment identities, gradient vs finite differences).

// [[Rcpp::export]]
Rcpp::NumericVector gpois_logpmf_cpp(Rcpp::NumericVector y, double mu, double lambda) {
  Rcpp::NumericVector out(y.size());
  double e1 = std::log(mu), e2 = std::log(lambda / (1.0 - lambda)), lp;
  for (int i = 0; i < y.size(); i++) { gpois_lpg(y[i], e1, e2, lp, NULL); out[i] = lp; }
  return out;
}
// [[Rcpp::export]]
Rcpp::NumericVector gwaring_logpmf_cpp(Rcpp::NumericVector y, double mu, double k, double rho) {
  Rcpp::NumericVector out(y.size());
  double e1 = std::log(mu), e2 = std::log(k), e3 = std::log(rho), lp;
  for (int i = 0; i < y.size(); i++) { gwaring_lpg(y[i], e1, e2, e3, lp, NULL); out[i] = lp; }
  return out;
}
// [[Rcpp::export]]
Rcpp::NumericVector plnorm_logpmf_cpp(Rcpp::NumericVector y, double mu, double sigma) {
  Rcpp::NumericVector out(y.size());
  double e1 = std::log(mu), e2 = std::log(sigma), lp;
  for (int i = 0; i < y.size(); i++) { plnorm_lpg(y[i], e1, e2, lp, NULL); out[i] = lp; }
  return out;
}

// Gradient of log pmf wrt the linear predictors at one y, for the finite-difference test.
// [[Rcpp::export]]
Rcpp::NumericVector cf_grad_cpp(std::string fam, double y, Rcpp::NumericVector eta) {
  double lp, g[3] = {0, 0, 0};
  int np = 2;
  if (fam == "gpois")        gpois_lpg(y, eta[0], eta[1], lp, g);
  else if (fam == "plnorm")  plnorm_lpg(y, eta[0], eta[1], lp, g);
  else if (fam == "gwaring") { gwaring_lpg(y, eta[0], eta[1], eta[2], lp, g); np = 3; }
  else Rcpp::stop("unknown family");
  Rcpp::NumericVector out(np + 1);
  out[0] = lp;
  for (int i = 0; i < np; i++) out[i + 1] = g[i];
  return out;
}
