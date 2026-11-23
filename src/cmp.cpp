// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
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
// wrong by 2%-100% (at nu = 8.3e-5, u = 100 the exact mean is 3159.09 -- a brute-force sum --
// while the asymptotic says 6124). Nor is u the mean there. So Z is computed EXACTLY, by
// log-sum-exp over j = 0..J with J ~ mode + 12 widths; J reaches ~1e5-2e5 on the largest series.
//
// WHY THAT IS STILL AFFORDABLE. nu is intercept-only in every model fitted here, so the exact
// sums need not be recomputed per observation. Once per likelihood evaluation we tabulate, on a
// lattice in v = log mu, the two functions the likelihood needs,
//
//     L(v) = log lambda(mu, nu)        P(v) = log Z(lambda(mu, nu), nu) / mu,
//
// and interpolate them by cubic Hermite with EXACT node slopes, which the same sums supply:
//
//     d log Z / d loglambda = E[Y]      d E[Y] / d loglambda = Var[Y]
//     d log Z / d nu        = -E[log Y!]
//     d E[Y]  / d nu        = -Cov(Y, log Y!)
//     d Var[Y]/ d loglambda = E[(Y-EY)^3],  d Var[Y] / d nu = -E[(Y-EY)^2 (log Y! - E log Y!)]
//
// so dL/dv = mu/Var, d(mu P)/dv = mu^2/Var, and the nu-derivatives of every node value and slope
// follow at fixed v (e.g. dL/dnu = Cov/Var). The lattice is v_k = k h with h = 1/64, anchored at
// multiples of h: the nodes depend neither on nu nor on the coefficients (only which nodes exist
// changes, never the interval that brackets an observation). The objective
//     log p(y) = y L(v) - nu log y! - mu P(v)
// is therefore an explicit smooth function of (log mu, nu) whose derivatives are those of the
// interpolants, and the analytic score IS the derivative of the objective that is optimised. (A
// lattice in log lambda whose ends tracked the extreme fitted means made the objective depend on
// the coefficients through the lattice itself, and the score differed from it by ~1e-5.) At a
// node the score reduces to the exact formulas
//
//     d log p / d log mu = (y - mu) mu / Var[Y]
//     d log p / d nu     = E[log Y!] - log y! + (y - mu) Cov(Y, log Y!) / Var[Y]
//
// and P varies slowly (it tends to 1 as mu -> 0 and to nu as mu -> Inf). For laws that are
// not too concentrated that keeps the interpolation error small: a fine scan of mu in
// [0.05, 400] gives at most 3e-8 in log p wherever Var[Y] >= 2, 3e-8 for nu <= 4 and 6e-6 for
// nu <= 8 (near mu = 1, where Var[Y] dips to 0.1).
//
// THE GRID IS VALIDATED, NOT TRUSTED. A strongly underdispersed law (Var[Y] well below 1) sits
// on one or two integers, and log lambda(mu) is then nearly singular at integer means: it rises
// by about nu log((j+1)/j) across a transition far narrower than the lattice, so no cubic can
// follow it (at nu = 1000, mu = 100.2 the gridded law had total mass 1.1, and at nu >= 64 the
// node at mu = 1 could not even be solved). An interval whose end nodes both have Var[Y] >= 2 is
// accepted when its counts lie within 100 sd of the mean (the discreteness ripple, of relative
// size ~exp(-2 pi^2 Var), is then negligible; the smooth part is accurate to ~5e-7 there even
// at nu = 1000). Any other interval that holds an observation is checked at its midpoint
// against an exact solve: if the interpolant's log-pmf error there could exceed CMP_TOL for the
// bulk of the law or for the counts observed in it, its observations use the exact sums instead,
// solved once per distinct mean. d0 and d12 share one plan, so values and scores always come
// from the same objective, and wherever the grid validates the result is unchanged.
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

static inline double cmp_clamp(double x, double lo, double hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

// The CENTRED mean equation at a target mean mu, sum_j (j - mu) w_j = 0. A strongly
// underdispersed law puts almost all its mass on the integer nearest its mean, and at an
// integer mean E[Y] then rounds to mu over a wide range of loglam (at mu = 1, nu = 128 for
// every loglam in [36.7, 52.0], about the root 44.36), so comparing the rounded mean with mu
// cannot locate the root. The positive and negative parts of the centred sum are therefore
// kept apart, as logs of their totals relative to Z: the j = mu term vanishes exactly, the
// factors j - mu near an integer mean are exact (Sterbenz), and each part keeps full
// relative accuracy however small it is. E[Y] < mu iff lpos < lneg, and
// (E[Y] - mu) / E[(Y - mu)^2] estimates the distance to the root in loglam.
struct CmpCentre {
  double mu;                  // target mean (input)
  double lpos, lneg;          // log sum_{j > mu} (j - mu) p_j and log sum_{j < mu} (mu - j) p_j
  double lv2;                 // log E[(Y - mu)^2]
};

// The three logs from their linear sums, or, where a part is too small for those (a law
// more concentrated than about e^-640 on mu), by log-sum-exp over the summands.
static void cmp_centre_logs(CmpCentre& C, const std::vector<double>& f, double lZ,
                            double pos, double neg, double v2)
{
  const double SAFE = 1e-280;
  if (pos >= SAFE && neg >= SAFE && v2 >= SAFE) {
    C.lpos = std::log(pos); C.lneg = std::log(neg); C.lv2 = std::log(v2);
    return;
  }
  const double NINF = -std::numeric_limits<double>::infinity();
  double mp = NINF, mn = NINF, mv = NINF;
  int J = (int)f.size() - 1;
  for (int j = 0; j <= J; j++) {
    double e = j - C.mu;
    if (e == 0.0) continue;
    double t = std::log(std::fabs(e)) + f[j];
    if (e > 0.0) mp = std::max(mp, t); else mn = std::max(mn, t);
    mv = std::max(mv, t + std::log(std::fabs(e)));
  }
  double sp = 0.0, sn = 0.0, sv = 0.0;
  for (int j = 0; j <= J; j++) {
    double e = j - C.mu;
    if (e == 0.0) continue;
    double t = std::log(std::fabs(e)) + f[j];
    if (e > 0.0) sp += std::exp(t - mp); else sn += std::exp(t - mn);
    sv += std::exp(t + std::log(std::fabs(e)) - mv);
  }
  C.lpos = (mp > NINF) ? mp + std::log(sp) - lZ : NINF;
  C.lneg = (mn > NINF) ? mn + std::log(sn) - lZ : NINF;
  C.lv2  = (mv > NINF) ? mv + std::log(sv) - lZ : NINF;
}

// Exact sums at (loglam, nu). Returns false if the required truncation exceeds CMP_JMAX.
// When K3/M21 are supplied they receive E[(Y-EY)^3] and E[(Y-EY)^2 (log Y! - E log Y!)];
// when C is supplied it receives the centred sums about C->mu (see CmpCentre).
static bool cmp_sums(double ll, double nu,
                     double& logZ, double& EY, double& ElgY, double& VarY, double& Cov,
                     double* K3 = NULL, double* M21 = NULL, CmpCentre* C = NULL)
{
  // f(j) = j*ll - nu*lgamma(j+1) is strictly concave (f'' = -nu*trigamma(j+1) < 0), so it has
  // a single mode. Truncation is by the ACTUAL log-summand, not a Gaussian-width heuristic:
  // when nu is tiny the summand decays like exp(-nu j log j), far slower than sqrt(2 j*/nu)
  // would suggest -- at nu = 8e-5 the mode sits at j ~ 2 but the mean is ~1800, and a
  // 12-width rule truncates while the summand is still ~15% of its max. We walk j upward,
  // tracking the running max, and stop once past the mode AND the summand has fallen LOG_TOL
  // below it (exp(-LOG_TOL) negligible against a max term of order 1).
  // Fast infeasibility bail. The mode is where f'(j) = ll - nu*digamma(j+1) = 0, i.e.
  // digamma(j+1) = ll/nu, so j* ~ exp(ll/nu) for large argument. If that already exceeds
  // CMP_JMAX the entire computable range lies on the rising part of f and the sum cannot
  // converge -- return at once rather than walk 4e6 terms to discover it. Without this the
  // optimiser's probes into infeasible (lambda, nu) regions each cost the full JMAX walk,
  // which on overdispersed data (small nu) turns a fit into hours. lambda <= 1 (ll <= 0) has
  // its mode at 0 and is always fine.
  if (ll > 0.0 && ll / nu > std::log(CMP_JMAX)) return false;

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

  double S = 0.0, s1 = 0.0, sl = 0.0;
  for (int j = 0; j <= J; j++) {
    double w = std::exp(f[j] - fmax);
    S += w; s1 += j * w; sl += lg[j] * w;
  }
  if (!(S > 0.0) || !std::isfinite(S)) return false;
  logZ = fmax + std::log(S);
  EY   = s1 / S;
  ElgY = sl / S;
  VarY = Cov = 0.0;
  double k3 = 0.0, m21 = 0.0;
  for (int j = 0; j <= J; j++) {
    double w = std::exp(f[j] - fmax) / S;
    double d = j - EY, d2 = d * d;
    VarY += w * d2;
    Cov += w * d * (lg[j] - ElgY);
    if (K3) { k3 += w * d2 * d; m21 += w * d2 * (lg[j] - ElgY); }
  }
  if (K3) { *K3 = k3; *M21 = m21; }
  if (C) {
    double cpos = 0.0, cneg = 0.0, cv2 = 0.0;
    for (int j = 0; j <= J; j++) {
      double w = std::exp(f[j] - fmax) / S, e = j - C->mu;
      if (e > 0.0) cpos += e * w; else cneg -= e * w;
      cv2 += e * e * w;
    }
    cmp_centre_logs(*C, f, logZ, cpos, cneg, cv2);
  }
  if (!(VarY > 0.0) || !std::isfinite(VarY)) return false;
  return true;
}

static const double CMP_ROOT_TOL = 1e-12;     // accepted distance to the centred root, x (1 + |loglam|)

// Does loglam lie within CMP_ROOT_TOL (1 + |loglam|) of the root of the centred mean
// equation? The Newton correction (E[Y] - mu) / E[(Y - mu)^2], formed from the centred
// parts, estimates the distance. The solvers below locate the root far more closely than
// this wherever the rounded mean resolves it, and their result is then kept as it is.
static inline bool cmp_root_near(const CmpCentre& C, double ll)
{
  double d = std::exp(C.lpos - C.lv2) - std::exp(C.lneg - C.lv2);
  return std::fabs(d) <= CMP_ROOT_TOL * (1.0 + std::fabs(ll));
}

// Solve the centred mean equation by bisection on its sign over the bracket of
// cmp_ll_from_mu, down to adjacent doubles. The fallback of both solvers below, so the
// exact path and the grid nodes share one root. A loglam whose truncation fits but whose
// variance underflows (an integer mean at nu beyond ~2000) still has a well-defined sign and
// is bisected on, not taken as infeasible; the root itself must have usable sums.
static bool cmp_ll_centred(double mu, double nu, double& ll_out)
{
  double lo = std::min(-60.0, std::log(mu) - 2.0), hi = nu * std::log1p(mu);
  hi = std::min(hi, nu * std::log(CMP_JMAX));
  double logZ, EY, ElgY, VarY, Cov;
  CmpCentre C;
  C.mu = mu;
  C.lpos = NAN;
  cmp_sums(hi, nu, logZ, EY, ElgY, VarY, Cov, NULL, NULL, &C);
  if (!std::isnan(C.lpos) && C.lpos < C.lneg) return false;     // mu is not attained
  for (int i = 0; i < 400; i++) {
    double mid = 0.5 * (lo + hi);
    if (!(mid > lo && mid < hi)) break;
    C.lpos = NAN;
    cmp_sums(mid, nu, logZ, EY, ElgY, VarY, Cov, NULL, NULL, &C);
    if (std::isnan(C.lpos)) { hi = mid; continue; }             // infeasible: lambda too large
    if (C.lpos < C.lneg) lo = mid; else hi = mid;
  }
  ll_out = 0.5 * (lo + hi);
  return cmp_sums(ll_out, nu, logZ, EY, ElgY, VarY, Cov, NULL, NULL, &C) &&
    std::fabs(std::exp(C.lpos) - std::exp(C.lneg)) <= 1e-8 * mu;
}

// Invert E[Y](loglam; nu) = mu. E[Y] is strictly increasing (dE/dll = Var > 0). The
// bisection on the rounded mean is kept wherever the centred sums confirm its root;
// elsewhere (means within rounding of an integer at large nu) the centred root is used.
static bool cmp_ll_from_mu(double mu, double nu, double& ll_out)
{
  if (!(mu > 0.0 && mu < CMP_JMAX && nu > 0.0) ||
      !std::isfinite(mu) || !std::isfinite(nu)) return false;
  double lo = std::min(-60.0, std::log(mu) - 2.0), hi = nu * std::log1p(mu);
  double logZ, EY, ElgY, VarY, Cov;
  // An infeasible upper endpoint still bounds every attainable mean. Never
  // shrink it to a computable but nonattaining value.
  hi = std::min(hi, nu * std::log(CMP_JMAX));
  if (cmp_sums(hi, nu, logZ, EY, ElgY, VarY, Cov) && EY < mu) return false;
  for (int i = 0; i < 200; i++) {
    double mid = 0.5 * (lo + hi);
    if (!cmp_sums(mid, nu, logZ, EY, ElgY, VarY, Cov)) { hi = mid; continue; }
    if (EY < mu) lo = mid; else hi = mid;
    if (hi - lo < 1e-13 * (1.0 + std::fabs(hi))) break;
  }
  ll_out = 0.5 * (lo + hi);
  CmpCentre C;
  C.mu = mu;
  if (!cmp_sums(ll_out, nu, logZ, EY, ElgY, VarY, Cov, NULL, NULL, &C) ||
      !(std::fabs(EY - mu) <= 1e-8 * mu)) return false;
  return cmp_root_near(C, ll_out) || cmp_ll_centred(mu, nu, ll_out);
}

// ---------------------------------------------------------------------- the grid ----
static const double CMP_H = 1.0 / 64.0;       // lattice spacing in v = log mu
static const double CMP_TOL = 1e-5;           // largest accepted grid error in log p
static const double CMP_VAR_SAFE = 2.0;       // node variance that needs no check ...
static const double CMP_SD_SAFE = 100.0;      // ... for counts within this many sd of mu

// Solve E[Y](loglam; nu) = mu by Newton's method on log E (derivative Var/E), safeguarded by
// bisection once a bracket is known, starting from `ll`. Where the law is too concentrated
// for the rounded mean to resolve the root, the root is checked against the centred sums
// and, if they reject it, replaced by the centred root, as in cmp_ll_from_mu. On success the
// exact sums at the returned loglam are left in the outputs.
static bool cmp_ll_newton(double mu, double nu, double& ll,
                          double& logZ, double& EY, double& ElgY, double& VarY, double& Cov,
                          double& K3, double& M21)
{
  const double INF = std::numeric_limits<double>::infinity();
  double lo = -INF, hi = INF, target = std::log(mu);
  bool ok = false;
  for (int it = 0; it < 200; it++) {
    if (!cmp_sums(ll, nu, logZ, EY, ElgY, VarY, Cov, &K3, &M21)) {
      hi = ll;                                       // infeasible: lambda too large
      ll = std::isfinite(lo) ? 0.5 * (lo + hi) : ll - 1.0;
      continue;
    }
    double g = std::log(EY) - target;
    if (std::fabs(g) <= 4e-16) { ok = true; break; }
    if (g < 0.0) lo = ll; else hi = ll;
    double nl = ll - g * EY / VarY;
    if (!(nl > lo && nl < hi))
      nl = (std::isfinite(lo) && std::isfinite(hi)) ? 0.5 * (lo + hi)
                                                  : (g < 0.0 ? ll + 1.0 : ll - 1.0);
    if (std::fabs(nl - ll) <= 1e-15 * (1.0 + std::fabs(ll))) {
      ll = nl;
      ok = cmp_sums(ll, nu, logZ, EY, ElgY, VarY, Cov, &K3, &M21) &&
        std::fabs(EY / mu - 1.0) <= 1e-12;
      break;
    }
    ll = nl;
  }
  if (!ok) return false;
  // The rounded mean resolves the root to CMP_ROOT_TOL (1 + |loglam|) unless the variance is
  // tiny next to its rounding (64 ulps of mu); only then is the root checked against the
  // centred sums, so ordinary nodes cost nothing extra.
  const double EPS = std::numeric_limits<double>::epsilon();
  if (VarY * CMP_ROOT_TOL * (1.0 + std::fabs(ll)) >= 64.0 * EPS * mu) return true;
  double lz, ey, elg, vy, cv;
  CmpCentre C;
  C.mu = mu;
  if (cmp_sums(ll, nu, lz, ey, elg, vy, cv, NULL, NULL, &C) && cmp_root_near(C, ll))
    return true;
  return cmp_ll_centred(mu, nu, ll) &&
    cmp_sums(ll, nu, logZ, EY, ElgY, VarY, Cov, &K3, &M21) &&
    std::fabs(EY / mu - 1.0) <= 1e-12;
}

// Node k of the lattice carries L, P and their v-slopes (Ls, Ps), the nu-derivatives of
// all four at fixed v (Ln, Pn, Lsn, Psn), Var[Y], and whether it could be solved (ok).
struct CmpGrid {
  double nu;
  int k0;
  std::vector<double> L, Ls, Ln, Lsn, P, Ps, Pn, Psn, Var;
  std::vector<char> ok;
};

// Solves every node it can; an interval with an unsolved end is left to the exact path.
// Returns false only when no node could be solved.
static bool cmp_build(CmpGrid& G, double nu, double mu_lo, double mu_hi)
{
  int k0 = (int)std::floor(std::log(mu_lo) / CMP_H);
  int k1 = (int)std::floor(std::log(mu_hi) / CMP_H) + 1;
  int n = k1 - k0 + 1;
  G.nu = nu; G.k0 = k0;
  G.L.resize(n); G.Ls.resize(n); G.Ln.resize(n); G.Lsn.resize(n);
  G.P.resize(n); G.Ps.resize(n); G.Pn.resize(n); G.Psn.resize(n);
  G.Var.assign(n, 0.0); G.ok.assign(n, 0);
  double ll = 0.0;
  bool any = false;
  for (int i = 0; i < n; i++) {
    double mu = std::exp((k0 + i) * CMP_H);
    if (!(mu > 0.0 && mu < CMP_JMAX) || !std::isfinite(mu)) continue;
    if (i == 0 || !G.ok[i - 1]) {
      // Starting value from the two regimes: lambda ~ mu for small means, lambda^(1/nu)
      // ~ mu + (nu-1)/(2 nu) for large ones, and the geometric limit lambda ~ mu/(1+mu).
      double g1 = mu + (nu - 1.0) / (2.0 * nu);
      ll = (mu <= 1.0) ? std::log(mu) : (g1 > 0.5 * mu ? nu * std::log(g1)
                                                        : std::log(mu / (1.0 + mu)));
    } else {
      ll += CMP_H * G.Ls[i - 1];                     // predictor along dL/dv
    }
    double logZ, EY, ElgY, VarY, Cov, K3, M21;
    if (!cmp_ll_newton(mu, nu, ll, logZ, EY, ElgY, VarY, Cov, K3, M21)) continue;
    double Ln = Cov / VarY;                          // d loglam / d nu at fixed mean
    double dVar = K3 * Ln - M21;                     // d Var / d nu at fixed mean
    double Phi_n = EY * Ln - ElgY;                   // d log Z / d nu at fixed mean
    G.L[i] = ll;
    G.Ls[i] = EY / VarY;
    G.Ln[i] = Ln;
    G.Lsn[i] = -EY * dVar / (VarY * VarY);
    G.P[i] = logZ / mu;
    G.Ps[i] = EY * EY / (VarY * mu) - G.P[i];
    G.Pn[i] = Phi_n / mu;
    G.Psn[i] = G.Lsn[i] * EY / mu - G.Pn[i];
    G.Var[i] = VarY;
    G.ok[i] = 1;
    any = true;
  }
  return any;
}

// ---------------------------------------------------------------- per-observation ----
// Cubic Hermite value, v-derivative and nu-derivative on the bracketing interval.
static inline void cmp_herm(const std::vector<double>& F, const std::vector<double>& Fs,
                            const std::vector<double>& Fn, const std::vector<double>& Fsn,
                            int i, double s, double& f, double& fv, double& fn)
{
  double s2 = s * s, s3 = s2 * s, h = CMP_H;
  double h00 = 2 * s3 - 3 * s2 + 1, h10 = s3 - 2 * s2 + s, h01 = -2 * s3 + 3 * s2, h11 = s3 - s2;
  double d00 = 6 * s2 - 6 * s, d10 = 3 * s2 - 4 * s + 1, d01 = -6 * s2 + 6 * s, d11 = 3 * s2 - 2 * s;
  f  = h00 * F[i] + h10 * h * Fs[i] + h01 * F[i + 1] + h11 * h * Fs[i + 1];
  fv = (d00 * F[i] + d10 * h * Fs[i] + d01 * F[i + 1] + d11 * h * Fs[i + 1]) / h;
  fn = h00 * Fn[i] + h10 * h * Fsn[i] + h01 * Fn[i + 1] + h11 * h * Fsn[i + 1];
}

// The lattice interval that brackets mu (clamped to the grid's ends).
static inline int cmp_interval(const CmpGrid& G, double mu)
{
  int i = (int)std::floor(std::log(mu) / CMP_H) - G.k0;
  int last = (int)G.L.size() - 2;
  if (i < 0) i = 0;
  if (i > last) i = last;
  return i;
}

static bool cmp_lpg_at(const CmpGrid& G, double y, double mu, double nu,
                       double& lp, double* g)
{
  double v = std::log(mu);
  int i = cmp_interval(G, mu);
  double s = (v - (G.k0 + i) * CMP_H) / CMP_H;
  double L, Lv, Ln, P, Pv, Pn;
  cmp_herm(G.L, G.Ls, G.Ln, G.Lsn, i, s, L, Lv, Ln);
  cmp_herm(G.P, G.Ps, G.Pn, G.Psn, i, s, P, Pv, Pn);
  double lgy = std::lgamma(y + 1.0);
  lp = y * L - nu * lgy - mu * P;
  if (!std::isfinite(lp)) return false;
  if (!g) return true;
  g[0] = y * Lv - mu * (Pv + P);                         // d/d log mu
  g[1] = nu * (y * Ln - lgy - mu * Pn);                  // d/d log nu
  return true;
}

// Exact sums at the mean mu: the same mean inversion and sums as dcmp().
struct CmpExact {
  bool ok;
  double ll, logZ, EY, ElgY, VarY, Cov;
};

static void cmp_exact_solve(double mu, double nu, CmpExact& E)
{
  E.ok = cmp_ll_from_mu(mu, nu, E.ll) &&
    cmp_sums(E.ll, nu, E.logZ, E.EY, E.ElgY, E.VarY, E.Cov);
}

static bool cmp_lpg_exact(const CmpExact& E, double y, double nu, double& lp, double* g)
{
  if (!E.ok) return false;
  double lgy = std::lgamma(y + 1.0);
  lp = y * E.ll - nu * lgy - E.logZ;
  if (!std::isfinite(lp)) return false;
  if (!g) return true;
  g[0] = (y - E.EY) * E.EY / E.VarY;
  g[1] = nu * ((E.ElgY - lgy) + (y - E.EY) * E.Cov / E.VarY);
  return true;
}

// Slow, correct fallback when nu is not constant across observations.
static bool cmp_lpg_slow(double y, double mu, double nu, double& lp, double* g)
{
  CmpExact E;
  cmp_exact_solve(mu, nu, E);
  return cmp_lpg_exact(E, y, nu, lp, g);
}

// Is the interpolant on interval i accurate for counts in [ylo, yhi]? See the header: both
// end nodes solved, and either both variances >= CMP_VAR_SAFE with the counts within
// CMP_SD_SAFE sd of the mean, or, against an exact solve at the midpoint, the log-pmf error
// there plus half an interval times the error of its v-slope at most CMP_TOL, for the bulk
// of the law and these counts. The slope term matters where the discreteness ripple (period
// 1/mu in v) is not resolved by the lattice: the value can match at the midpoint by chance,
// but not the value and the slope together. The count limit matters because the error grows
// like |y - mu| times that of log lambda: at nu = 1000, mu = 4800 (Var 4.8) it is 5e-7 in the
// bulk but 2.5e-5 at y = 0, some 2000 sd away.
static bool cmp_interval_ok(const CmpGrid& G, int i, double ylo, double yhi)
{
  if (!G.ok[i] || !G.ok[i + 1]) return false;
  double mu = std::exp((G.k0 + i + 0.5) * CMP_H);
  if (G.Var[i] >= CMP_VAR_SAFE && G.Var[i + 1] >= CMP_VAR_SAFE) {
    double sd = std::sqrt(std::min(G.Var[i], G.Var[i + 1]));
    if (ylo >= mu - CMP_SD_SAFE * sd && yhi <= mu + CMP_SD_SAFE * sd) return true;
  }
  double L, Lv, Ln, P, Pv, Pn;
  cmp_herm(G.L, G.Ls, G.Ln, G.Lsn, i, 0.5, L, Lv, Ln);
  cmp_herm(G.P, G.Ps, G.Pn, G.Psn, i, 0.5, P, Pv, Pn);
  double ll = L, logZ, EY, ElgY, VarY, Cov, K3, M21;
  if (!std::isfinite(L) || !std::isfinite(P) ||
      !cmp_ll_newton(mu, G.nu, ll, logZ, EY, ElgY, VarY, Cov, K3, M21)) return false;
  double Pm = logZ / mu, dL = L - ll, dP = P - Pm;
  double dLv = Lv - EY / VarY, dPv = Pv - (EY * EY / (VarY * mu) - Pm);
  double sd = std::sqrt(VarY);
  double lo = std::min(ylo, std::max(0.0, std::floor(mu - 8.0 * sd) - 2.0));
  double hi = std::max(yhi, std::ceil(mu + 8.0 * sd) + 2.0);
  double err = 0.0;
  for (double y : {lo, hi})
    err = std::max(err, std::fabs(y * dL - mu * dP) +
                        0.5 * CMP_H * std::fabs(y * dLv - mu * (dPv + dP)));
  return err <= CMP_TOL;
}

// How each observation is evaluated when nu is constant: on the grid where its interval
// validates (how[j] = interval), otherwise from the exact sums at its mean, solved once per
// distinct mean (how[j] = -1 - index into ex). d0 and d12 both evaluate through this plan.
struct CmpPlan {
  CmpGrid G;
  std::vector<int> how;
  std::vector<CmpExact> ex;
};

static void cmp_plan(const std::vector<double>& mu, double nu, const arma::vec& yvec,
                     CmpPlan& P)
{
  const double INF = std::numeric_limits<double>::infinity();
  int n = (int)mu.size();
  double mlo = INF, mhi = 0.0;
  for (int j = 0; j < n; j++) {
    if (mu[j] < mlo) mlo = mu[j];
    if (mu[j] > mhi) mhi = mu[j];
  }
  bool grid = cmp_build(P.G, nu, mlo, mhi);
  P.how.assign(n, -1);
  if (grid) {
    int nint = (int)P.G.L.size() - 1;
    std::vector<double> ylo(nint, INF), yhi(nint, -INF);
    for (int j = 0; j < n; j++) {
      int i = cmp_interval(P.G, mu[j]);
      P.how[j] = i;
      if (yvec[j] < ylo[i]) ylo[i] = yvec[j];
      if (yvec[j] > yhi[i]) yhi[i] = yvec[j];
    }
    std::vector<char> good(nint, 0);
    for (int i = 0; i < nint; i++)
      if (ylo[i] <= yhi[i]) good[i] = cmp_interval_ok(P.G, i, ylo[i], yhi[i]);
    for (int j = 0; j < n; j++)
      if (!good[P.how[j]]) P.how[j] = -1;
  }
  std::vector<double> emu;
  for (int j = 0; j < n; j++) if (P.how[j] < 0) emu.push_back(mu[j]);
  if (emu.empty()) return;
  std::sort(emu.begin(), emu.end());
  emu.erase(std::unique(emu.begin(), emu.end()), emu.end());
  P.ex.resize(emu.size());
  for (size_t k = 0; k < emu.size(); k++) cmp_exact_solve(emu[k], nu, P.ex[k]);
  for (int j = 0; j < n; j++)
    if (P.how[j] < 0) {
      size_t k = std::lower_bound(emu.begin(), emu.end(), mu[j]) - emu.begin();
      P.how[j] = -1 - (int)k;
    }
}

static inline bool cmp_lpg_plan(const CmpPlan& P, int j, double y, double mu, double nu,
                                double& lp, double* g)
{
  int h = P.how[j];
  return h >= 0 ? cmp_lpg_at(P.G, y, mu, nu, lp, g)
                : cmp_lpg_exact(P.ex[-1 - h], y, nu, lp, g);
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
    std::vector<double> mu(n);
    for (arma::uword j = 0; j < n; j++) mu[j] = cmp_clamp(std::exp(e1[j]), CMP_MU_LO, CMP_MU_HI);
    CmpPlan P;
    cmp_plan(mu, nu, yvec, P);
    for (arma::uword j = 0; j < n; j++) {
      if (!cmp_lpg_plan(P, (int)j, yvec[j], mu[j], nu, lp, NULL)) return INF;
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
  double lp, gg[2];

  bool constnu = cmp_nu_constant(e2);
  CmpPlan P;
  double nu0 = cmp_clamp(std::exp(e2[0]), CMP_NU_LO, CMP_NU_HI);
  std::vector<double> muv(n);
  for (arma::uword j = 0; j < n; j++) muv[j] = cmp_clamp(std::exp(e1[j]), CMP_MU_LO, CMP_MU_HI);
  if (constnu) cmp_plan(muv, nu0, yvec, P);
  for (arma::uword j = 0; j < n; j++) {
    double mu = muv[j];
    double nu = constnu ? nu0 : cmp_clamp(std::exp(e2[j]), CMP_NU_LO, CMP_NU_HI);
    bool ok = constnu ? cmp_lpg_plan(P, (int)j, yvec[j], mu, nu, lp, gg)
                      : cmp_lpg_slow(yvec[j], mu, nu, lp, gg);
    if (!ok) { gg[0] = 0.0; gg[1] = 0.0; }
    if (std::exp(e1[j]) < CMP_MU_LO || std::exp(e1[j]) > CMP_MU_HI) gg[0] = 0.0;
    if (std::exp(e2[j]) < CMP_NU_LO || std::exp(e2[j]) > CMP_NU_HI) gg[1] = 0.0;
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
    Rcpp::stop("CMP mean inversion or normalization exceeds numerical limits");
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
