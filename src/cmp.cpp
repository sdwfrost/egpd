// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <vector>

// Conway-Maxwell-Poisson (CMP), mean-parameterised, for the native GAM fitter.
//
//   P(Y = y) = lambda^y / ((y!)^nu Z(lambda, nu)),   Z = sum_{j>=0} lambda^j / (j!)^nu
//
// nu < 1 gives overdispersion, nu = 1 the Poisson, nu > 1 underdispersion. The tail decays
// like exp(-nu y log y), i.e. LIGHTER THAN GEOMETRIC for any nu > 0, so CMP has extreme-value
// index xi = 0. It is included precisely as a control: if a family whose tail is
// asymptotically lighter than geometric fits these counts as well as the power-law families,
// the heavy-tail conclusion is weaker than it looks. (At the nu ~ 1e-5 these data force, the
// factorial damping only bites far beyond the largest observed count, so whether that happens
// is an empirical question, not a foregone one.)
//
// WHY NOT THE ASYMPTOTIC NORMALISING CONSTANT. Gaunt et al.'s expansion for Z is governed by
// u*nu, not by u = lambda^(1/nu). These series need u*nu in [0.008, 0.4], where its mean is
// wrong by 2%-100% (at nu = 8.3e-5, u = 100 the exact mean is 3058 while the asymptotic says
// 6124). Nor is u the mean there. So Z is computed EXACTLY, by log-sum-exp over j = 0..J with
// J ~ mode + 12 widths; J reaches ~1e5-2e5 on the largest series.
//
// WHY THAT IS STILL AFFORDABLE. nu is intercept-only in every model fitted here, so the exact
// sums need not be recomputed per observation. Once per likelihood evaluation we build a grid
// over log lambda and interpolate. The interpolation is cubic Hermite with EXACT derivatives,
// because the same sums hand them to us:
//
//     d logZ / d loglambda = E[Y]
//     d E[Y]  / d loglambda = Var[Y]
//     d logZ / d nu         = -E[log Y!]
//     d E[Y]  / d nu        = -Cov(Y, log Y!)
//
// The mean parameterisation then costs nothing extra: E[Y] is strictly increasing in
// log lambda (its derivative is Var[Y] > 0), so lambda(mu) comes from inverting the same
// monotone interpolant. The scores follow by the chain rule,
//
//     d log p / d log mu = (y - mu) mu / Var[Y]
//     d log p / d nu     = E[log Y!] - log y! + (y - mu) Cov(Y, log Y!) / Var[Y]
//
// THE PRECONDITION IS CHECKED, NOT ASSUMED. If nu varies across observations (i.e. it has been
// given covariates) the grid is invalid, and the code falls back to a per-observation root
// solve. That is far slower, but it is never silently wrong.
//
// d12 packs per-observation NEGATIVE log-likelihood derivatives: cols [g1,g2, H11,H12,H22]
// (BHHH), matching .gH2.

static const double CMP_NU_LO = 1e-6,  CMP_NU_HI = 1e3;
static const double CMP_MU_LO = 1e-8,  CMP_MU_HI = 1e12;
static const double CMP_JMAX  = 4e6;          // refuse rather than hang
static const int    CMP_GRID  = 192;          // nodes in log lambda

static inline double cmp_clamp(double x, double lo, double hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

// Exact sums at (loglam, nu). Returns false if the required truncation exceeds CMP_JMAX.
static bool cmp_sums(double ll, double nu,
                     double& logZ, double& EY, double& ElgY, double& VarY, double& Cov)
{
  // f(j) = j*ll - nu*lgamma(j+1) is strictly concave (f'' = -nu*trigamma(j+1) < 0), so it has
  // a single mode. Truncation is by the ACTUAL log-summand, not a Gaussian-width heuristic:
  // when nu is tiny the summand decays like exp(-nu j log j), far slower than sqrt(2 j*/nu)
  // would suggest -- at nu = 8e-5 the mode sits at j ~ 2 but the mean is ~1800, and a
  // 12-width rule truncates while the summand is still ~15% of its max. We walk j upward,
  // tracking the running max, and stop once past the mode AND the summand has fallen LOG_TOL
  // below it (exp(-LOG_TOL) negligible against a max term of order 1).
  const double LOG_TOL = 45.0;                           // exp(-45) ~ 3e-20
  std::vector<double> f, lg;
  f.reserve(4096); lg.reserve(4096);
  double lgj = 0.0, fmax = -std::numeric_limits<double>::infinity();
  int jmode = 0;
  for (int j = 0; ; j++) {
    if (j > 0) lgj += std::log((double)j);               // lgamma(j+1), accumulated
    double fj = j * ll - nu * lgj;
    lg.push_back(lgj); f.push_back(fj);
    if (fj > fmax) { fmax = fj; jmode = j; }
    if (j > jmode && fj < fmax - LOG_TOL) break;         // tail is now negligible
    if (j >= (int)CMP_JMAX) return false;                // refuse rather than hang
  }
  int J = (int)f.size() - 1;

  double S = 0.0, s1 = 0.0, s2 = 0.0, sl = 0.0, sjl = 0.0;
  for (int j = 0; j <= J; j++) {
    double w = std::exp(f[j] - fmax);
    S += w; s1 += j * w; s2 += (double)j * j * w; sl += lg[j] * w; sjl += j * lg[j] * w;
  }
  if (!(S > 0.0) || !std::isfinite(S)) return false;
  logZ = fmax + std::log(S);
  EY   = s1 / S;
  ElgY = sl / S;
  VarY = s2 / S - EY * EY;
  Cov  = sjl / S - EY * ElgY;
  if (!(VarY > 0.0) || !std::isfinite(VarY)) return false;
  return true;
}

// Invert E[Y](loglam; nu) = mu. E[Y] is strictly increasing (dE/dll = Var > 0).
static bool cmp_ll_from_mu(double mu, double nu, double& ll_out)
{
  double lo = -60.0, hi = 60.0;
  double logZ, EY, ElgY, VarY, Cov;
  for (int i = 0; i < 200 && !cmp_sums(hi, nu, logZ, EY, ElgY, VarY, Cov); i++) hi *= 0.8;
  if (!cmp_sums(hi, nu, logZ, EY, ElgY, VarY, Cov)) return false;
  if (EY < mu) { ll_out = hi; return true; }             // saturated: best available
  for (int i = 0; i < 200; i++) {
    double mid = 0.5 * (lo + hi);
    if (!cmp_sums(mid, nu, logZ, EY, ElgY, VarY, Cov)) { hi = mid; continue; }
    if (EY < mu) lo = mid; else hi = mid;
    if (hi - lo < 1e-13 * (1.0 + std::fabs(hi))) break;
  }
  ll_out = 0.5 * (lo + hi);
  return true;
}

// ---------------------------------------------------------------------- the grid ----
struct CmpGrid {
  double nu;
  std::vector<double> ll, logZ, EY, ElgY, VarY, Cov;
};

static bool cmp_build(CmpGrid& G, double nu, double mu_lo, double mu_hi)
{
  double a, b;
  if (!cmp_ll_from_mu(mu_lo * 0.5, nu, a)) return false;
  if (!cmp_ll_from_mu(mu_hi * 2.0, nu, b)) return false;
  if (!(b > a)) b = a + 1e-6;
  int n = CMP_GRID;
  G.nu = nu;
  G.ll.resize(n); G.logZ.resize(n); G.EY.resize(n);
  G.ElgY.resize(n); G.VarY.resize(n); G.Cov.resize(n);
  for (int i = 0; i < n; i++) {
    double ll = a + (b - a) * i / (n - 1.0);
    if (!cmp_sums(ll, nu, G.logZ[i], G.EY[i], G.ElgY[i], G.VarY[i], G.Cov[i])) return false;
    G.ll[i] = ll;
  }
  return true;
}

// cubic Hermite on [x0,x1] with values y0,y1 and exact slopes m0,m1
static inline double hermite(double x, double x0, double x1, double y0, double y1,
                             double m0, double m1) {
  double h = x1 - x0, s = (x - x0) / h, s2 = s * s, s3 = s2 * s;
  return (2 * s3 - 3 * s2 + 1) * y0 + (s3 - 2 * s2 + s) * h * m0
       + (-2 * s3 + 3 * s2) * y1 + (s3 - s2) * h * m1;
}

static int cmp_bracket(const std::vector<double>& x, double v) {
  int lo = 0, hi = (int)x.size() - 1;
  if (v <= x[0]) return 0;
  if (v >= x[hi]) return hi - 1;
  while (hi - lo > 1) { int m = (lo + hi) / 2; if (x[m] <= v) lo = m; else hi = m; }
  return lo;
}

static void cmp_at_ll(const CmpGrid& G, double ll,
                      double& logZ, double& EY, double& ElgY, double& VarY, double& Cov)
{
  int i = cmp_bracket(G.ll, ll);
  double x0 = G.ll[i], x1 = G.ll[i + 1];
  logZ = hermite(ll, x0, x1, G.logZ[i], G.logZ[i + 1], G.EY[i], G.EY[i + 1]);     // dlogZ/dll = EY
  EY   = hermite(ll, x0, x1, G.EY[i],   G.EY[i + 1],   G.VarY[i], G.VarY[i + 1]); // dEY/dll = Var
  double s = (ll - x0) / (x1 - x0);
  ElgY = (1 - s) * G.ElgY[i] + s * G.ElgY[i + 1];
  VarY = (1 - s) * G.VarY[i] + s * G.VarY[i + 1];
  Cov  = (1 - s) * G.Cov[i]  + s * G.Cov[i + 1];
}

// Invert the monotone EY interpolant for loglam given mu.
static double cmp_ll_from_mu_grid(const CmpGrid& G, double mu)
{
  int lo = 0, hi = (int)G.EY.size() - 1;
  if (mu <= G.EY[0]) return G.ll[0];
  if (mu >= G.EY[hi]) return G.ll[hi];
  while (hi - lo > 1) { int m = (lo + hi) / 2; if (G.EY[m] <= mu) lo = m; else hi = m; }
  double ll = G.ll[lo], logZ, EY, ElgY, VarY, Cov;
  double last = G.ll[G.ll.size() - 1];
  for (int it = 0; it < 60; it++) {                      // Newton, derivative = Var
    cmp_at_ll(G, ll, logZ, EY, ElgY, VarY, Cov);
    double step = (EY - mu) / VarY;
    ll -= step;
    if (ll < G.ll[0]) ll = G.ll[0];
    if (ll > last)    ll = last;
    if (std::fabs(step) < 1e-14) break;
  }
  return ll;
}

// ---------------------------------------------------------------- per-observation ----
static bool cmp_lpg_at(const CmpGrid& G, double y, double mu, double nu,
                       double& lp, double* g)
{
  double ll = cmp_ll_from_mu_grid(G, mu);
  double logZ, EY, ElgY, VarY, Cov;
  cmp_at_ll(G, ll, logZ, EY, ElgY, VarY, Cov);
  double lgy = std::lgamma(y + 1.0);
  lp = y * ll - nu * lgy - logZ;
  if (!std::isfinite(lp)) return false;
  if (!g) return true;
  g[0] = (y - EY) * EY / VarY;                           // d/d log mu  (EY == mu at the node)
  g[1] = nu * ((ElgY - lgy) + (y - EY) * Cov / VarY);    // d/d log nu
  return true;
}

// Slow, correct fallback when nu is not constant across observations.
static bool cmp_lpg_slow(double y, double mu, double nu, double& lp, double* g)
{
  double ll;
  if (!cmp_ll_from_mu(mu, nu, ll)) return false;
  double logZ, EY, ElgY, VarY, Cov;
  if (!cmp_sums(ll, nu, logZ, EY, ElgY, VarY, Cov)) return false;
  double lgy = std::lgamma(y + 1.0);
  lp = y * ll - nu * lgy - logZ;
  if (!std::isfinite(lp)) return false;
  if (!g) return true;
  g[0] = (y - EY) * EY / VarY;
  g[1] = nu * ((ElgY - lgy) + (y - EY) * Cov / VarY);
  return true;
}

// ------------------------------------------------------------------- fitter glue ----

static void cmp_eta(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2,
                    const arma::uvec& dupid, int dcate, const Rcpp::List& offsets,
                    arma::vec& e1, arma::vec& e2)
{
  e1 = X1 * Rcpp::as<arma::vec>(pars[0]);
  e2 = X2 * Rcpp::as<arma::vec>(pars[1]);
  if (dcate == 1) { e1 = e1.elem(dupid); e2 = e2.elem(dupid); }
  arma::vec o0 = Rcpp::as<arma::vec>(offsets[0]); if (o0.n_elem > 0) e1 += o0;
  arma::vec o1 = Rcpp::as<arma::vec>(offsets[1]); if (o1.n_elem > 0) e2 += o1;
}

// True when nu is the same for every observation (the intercept-only case), which is what
// makes the grid valid. Checked, never assumed.
static bool cmp_nu_constant(const arma::vec& e2) {
  double m = e2[0];
  for (arma::uword j = 1; j < e2.n_elem; j++)
    if (std::fabs(e2[j] - m) > 1e-12 * (1.0 + std::fabs(m))) return false;
  return true;
}

// [[Rcpp::export]]
double cmp1d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2,
              arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets) {
  arma::vec e1, e2;
  cmp_eta(pars, X1, X2, dupid, dcate, offsets, e1, e2);
  const double INF = std::numeric_limits<double>::infinity();
  arma::uword n = yvec.n_elem;
  double nllh = 0.0, lp;

  if (cmp_nu_constant(e2)) {
    double nu = cmp_clamp(std::exp(e2[0]), CMP_NU_LO, CMP_NU_HI);
    double mlo = INF, mhi = 0.0;
    for (arma::uword j = 0; j < n; j++) {
      double mu = cmp_clamp(std::exp(e1[j]), CMP_MU_LO, CMP_MU_HI);
      if (mu < mlo) mlo = mu;
      if (mu > mhi) mhi = mu;
    }
    CmpGrid G;
    if (!cmp_build(G, nu, mlo, mhi)) return INF;
    for (arma::uword j = 0; j < n; j++) {
      double mu = cmp_clamp(std::exp(e1[j]), CMP_MU_LO, CMP_MU_HI);
      if (!cmp_lpg_at(G, yvec[j], mu, nu, lp, NULL)) return INF;
      nllh -= lp;
    }
    return nllh;
  }
  for (arma::uword j = 0; j < n; j++) {
    double mu = cmp_clamp(std::exp(e1[j]), CMP_MU_LO, CMP_MU_HI);
    double nu = cmp_clamp(std::exp(e2[j]), CMP_NU_LO, CMP_NU_HI);
    if (!cmp_lpg_slow(yvec[j], mu, nu, lp, NULL)) return INF;
    nllh -= lp;
  }
  return nllh;
}

// [[Rcpp::export]]
arma::mat cmp1d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2,
                  arma::vec yvec, const arma::uvec dupid, int dcate, const Rcpp::List& offsets) {
  arma::vec e1, e2;
  cmp_eta(pars, X1, X2, dupid, dcate, offsets, e1, e2);
  arma::uword n = yvec.n_elem;
  arma::mat out(n, 5, arma::fill::zeros);
  const double INF = std::numeric_limits<double>::infinity();
  double lp, gg[2];

  bool constnu = cmp_nu_constant(e2);
  CmpGrid G;
  double nu0 = cmp_clamp(std::exp(e2[0]), CMP_NU_LO, CMP_NU_HI);
  if (constnu) {
    double mlo = INF, mhi = 0.0;
    for (arma::uword j = 0; j < n; j++) {
      double mu = cmp_clamp(std::exp(e1[j]), CMP_MU_LO, CMP_MU_HI);
      if (mu < mlo) mlo = mu;
      if (mu > mhi) mhi = mu;
    }
    if (!cmp_build(G, nu0, mlo, mhi)) constnu = false;
  }
  for (arma::uword j = 0; j < n; j++) {
    double mu = cmp_clamp(std::exp(e1[j]), CMP_MU_LO, CMP_MU_HI);
    double nu = constnu ? nu0 : cmp_clamp(std::exp(e2[j]), CMP_NU_LO, CMP_NU_HI);
    bool ok = constnu ? cmp_lpg_at(G, yvec[j], mu, nu, lp, gg)
                      : cmp_lpg_slow(yvec[j], mu, nu, lp, gg);
    if (!ok) { gg[0] = 0.0; gg[1] = 0.0; }
    double g[2] = { -gg[0], -gg[1] };                    // NLL gradient
    out(j, 0) = g[0]; out(j, 1) = g[1];
    out(j, 2) = g[0] * g[0]; out(j, 3) = g[0] * g[1]; out(j, 4) = g[1] * g[1];
  }
  return out;
}

// ------------------------------------------------------------------- test hooks ----

// Exact sums at the NATURAL parameters, for validation against a brute-force sum in R.
// [[Rcpp::export]]
Rcpp::NumericVector cmp_sums_cpp(double loglambda, double nu) {
  double logZ, EY, ElgY, VarY, Cov;
  if (!cmp_sums(loglambda, nu, logZ, EY, ElgY, VarY, Cov))
    Rcpp::stop("cmp_sums failed (truncation exceeds CMP_JMAX)");
  return Rcpp::NumericVector::create(Rcpp::Named("logZ") = logZ, Rcpp::Named("EY") = EY,
                                     Rcpp::Named("ElgY") = ElgY, Rcpp::Named("VarY") = VarY,
                                     Rcpp::Named("Cov") = Cov);
}

// [[Rcpp::export]]
double cmp_ll_from_mu_cpp(double mu, double nu) {
  double ll;
  if (!cmp_ll_from_mu(mu, nu, ll)) Rcpp::stop("cmp_ll_from_mu failed");
  return ll;
}

// log pmf on the mean scale, via the SLOW exact path (no grid), for validation.
// [[Rcpp::export]]
Rcpp::NumericVector cmp_logpmf_cpp(Rcpp::NumericVector y, double mu, double nu) {
  Rcpp::NumericVector out(y.size());
  // (mu, nu) are scalar, so solve lambda and the sums ONCE and reuse them across all y --
  // otherwise dcmp(0:N) would re-run the mu->lambda bisection N times (quadratic).
  double ll, logZ, EY, ElgY, VarY, Cov;
  if (!cmp_ll_from_mu(mu, nu, ll) || !cmp_sums(ll, nu, logZ, EY, ElgY, VarY, Cov)) {
    for (int i = 0; i < y.size(); i++) out[i] = -std::numeric_limits<double>::infinity();
    return out;
  }
  for (int i = 0; i < y.size(); i++)
    out[i] = y[i] * ll - nu * std::lgamma(y[i] + 1.0) - logZ;
  return out;
}

// Analytic score at one y, via the slow exact path.
// [[Rcpp::export]]
Rcpp::NumericVector cmp_grad_cpp(double y, double mu, double nu) {
  double lp, g[2];
  if (!cmp_lpg_slow(y, mu, nu, lp, g)) Rcpp::stop("cmp_lpg_slow failed");
  return Rcpp::NumericVector::create(Rcpp::Named("logp") = lp, Rcpp::Named("g_logmu") = g[0],
                                     Rcpp::Named("g_lognu") = g[1]);
}

// [[Rcpp::export]]
Rcpp::NumericVector cmp_bounds_cpp() {
  return Rcpp::NumericVector::create(Rcpp::Named("mu_lo") = CMP_MU_LO, Rcpp::Named("mu_hi") = CMP_MU_HI,
                                     Rcpp::Named("nu_lo") = CMP_NU_LO, Rcpp::Named("nu_hi") = CMP_NU_HI);
}
