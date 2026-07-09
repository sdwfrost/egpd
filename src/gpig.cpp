// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <limits>

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

// ------------------------------------------------------------- fit methods ----
//
// The Zhu & Joe forward recursion costs O(y^2) per observation and is carried in
// raw probability space, seeded by p_0 = exp(b[(1-c)^a - 1]). Both properties bite
// at large counts: a likelihood evaluation costs sum_i y_i^2, and p_0 underflows to
// zero once b[(1-c)^a - 1] < -745 (mean of order 1e3), after which every p_k is
// identically zero and the NLL is +Inf. gpig_method selects among five evaluation
// routes, set from R by gpig_control():
//
//   0 "legacy"      raw unscaled recursion, exactly as before (reproducibility)
//   1 "recursion"   same recursion carried with dynamic rescaling and a log offset;
//                   exact, O(y^2), immune to the p_0 underflow
//   2 "trunc"       (1) plus truncation of the severity convolution. The severity
//                   r_j decays geometrically at rate c, so the convolution is
//                   accumulated from its largest terms downwards and stopped once a
//                   term falls below `eps` times the partial sum. The criterion must
//                   be relative: p_y itself may be of order e^-300, so an absolute
//                   bound on the dropped mass destroys all relative accuracy.
//   3 "saddlepoint" Daniels' saddlepoint applied to the closed-form CGF
//                   K(t) = b[(1-c)^a - (1-c e^t)^a]; O(1) per observation
//   4 "hybrid"      (default) exact recursion for y <= yswitch; above it, the
//                   saddlepoint, but only where it is trustworthy.
//
// Accuracy of the saddlepoint is governed by the large parameter b (the compound-
// Poisson rate is b[1-(1-c)^a]), NOT by y: at b = 300 the second-order log-pmf error
// is ~1e-6 at every y, while at b = 5 it grows with y. The asymptotic series is
// therefore policed directly by the size of its own first correction term,
// corr = kappa4/8 - 5 kappa3^2/24 = g/b: when |corr| > sptol the expansion has not
// converged and "hybrid" falls back to the exact recursion.
//
// No fixed-order recurrence exists for generic a: exp(g) is holonomic only when g'
// is rational, and here g' = abc(1-cs)^(a-1). At a = 1/2 the family reduces to PIG
// and the Bessel three-term recurrence applies; otherwise the pmf must be reached
// by convolution or by transform/saddlepoint methods.
enum GpigMethod { GPIG_LEGACY = 0, GPIG_RECURSION = 1, GPIG_TRUNC = 2,
                  GPIG_SADDLE = 3, GPIG_HYBRID = 4 };
static int    gpig_method  = GPIG_HYBRID;
static int    gpig_yswitch = 200;
static int    gpig_order   = 2;      // saddlepoint order: 1 or 2
static double gpig_eps     = 1e-12;  // relative severity-truncation tolerance
static double gpig_sptol   = 1e-3;   // max |corr| for "hybrid" to trust the saddlepoint
static int    gpig_ymax    = 1000;   // largest y at which "hybrid" may fall back to O(y^2)
static int    gpig_nquad   = 256;    // quadrature nodes for the normalising constant
// Normalisation is OPT-IN. It restores a proper pmf and sharpens the saddlepoint by 50-100x
// where the expansion has already converged, but it buys nothing at the c -> 1 boundary
// (the bias there is y-dependent, not the constant S), and the quadrature is costly. The
// likelihood is kept safe instead by the GPIG_CORR_MAX gate and the log p <= 0 invariant.
static int    gpig_normalize = 0;    // divide the saddlepoint pmf by sum_k phat(k)

// [[Rcpp::export]]
void gpig_set_opts(int method, int yswitch, int order, double eps, double sptol, int ymax,
                   int nquad, int normalize) {
  gpig_method  = method;
  gpig_yswitch = yswitch;
  gpig_order   = order;
  gpig_eps     = eps;
  gpig_sptol   = sptol;
  gpig_ymax    = ymax;
  gpig_nquad   = nquad;
  gpig_normalize = normalize;
}

// [[Rcpp::export]]
Rcpp::List gpig_get_opts() {
  return Rcpp::List::create(Rcpp::Named("method")  = gpig_method,
                            Rcpp::Named("yswitch") = gpig_yswitch,
                            Rcpp::Named("order")   = gpig_order,
                            Rcpp::Named("eps")     = gpig_eps,
                            Rcpp::Named("sptol")   = gpig_sptol,
                            Rcpp::Named("ymax")    = gpig_ymax,
                            Rcpp::Named("nquad")   = gpig_nquad,
                            Rcpp::Named("normalize") = gpig_normalize);
}

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

// ------------------------------------- exact recursion, carried in log space ---
//
// Identical mathematics to gpig_derivs(), but p_k and dp_k/d(a,b,c) are held only
// up to a common positive factor: the recursion is linear and homogeneous in the
// stacked vector (p_k, dp_k/da, dp_k/db, dp_k/dc), so rescaling every entry by the
// same constant leaves it invariant. We seed with p~_0 = 1 and logscale = log p_0,
// rescale whenever p~ leaves [1e-250, 1e250], and return
//     log p_y = logscale + log p~_y,   d log p_y/dtheta = (dp~_y/dtheta) / p~_y,
// both of which are invariant to the scale. This removes the p_0 underflow.
//
// use_trunc additionally truncates the severity convolution. The terms
// j * r_{k+1-j} * p_j are summed from j = k downwards, which visits the severity in
// increasing index and so in decreasing magnitude; the sum stops once every one of
// the four accumulators gains less than eps times its own partial sum. The stopping
// rule has to be relative rather than absolute: p_y may be of order e^-300, so
// bounding the dropped mass by eps * max(r) says nothing about the accuracy of
// log p_y.
//
// LIMITATION. A single common scale cannot represent the whole ladder once
// p_y / max_j p_j falls below ~1e-308: the convolution needs every p_j at once, and
// no rescaling can hold both ends in double precision. Rescaling is therefore capped
// so that no entry can overflow to Inf, and the far tail is allowed to underflow to
// zero, whereupon the routine reports failure by returning false with lpy = -Inf.
// It never returns a saturated finite value. "hybrid" then falls through to the
// saddlepoint, which is exactly the regime the saddlepoint handles well.
static bool gpig_recur(int y, double a, double b, double c, bool want_grad,
                       double& lpy, double& dla, double& dlb, double& dlc,
                       bool use_trunc, double eps)
{
  double omc = 1.0 - c, lomc = std::log(omc);
  double omca = std::pow(omc, a);
  double lp0 = b * (omca - 1.0);
  double d0a = b * omca * lomc;             // d log p0 / da
  double d0b = omca - 1.0;                  // d log p0 / db
  double d0c = -a * b * std::pow(omc, a - 1.0);

  if (y == 0) { lpy = lp0; dla = d0a; dlb = d0b; dlc = d0c; return true; }

  std::vector<double> p(y + 1), r(y + 1);
  std::vector<double> pa, pb, pc, dra, drc;
  if (want_grad) {
    pa.assign(y + 1, 0.0); pb.assign(y + 1, 0.0); pc.assign(y + 1, 0.0);
    dra.assign(y + 1, 0.0); drc.assign(y + 1, 0.0);
  }

  r[1] = (1.0 - a) * c;
  if (want_grad) { dra[1] = -c; drc[1] = 1.0 - a; }
  for (int j = 1; j <= y - 1; j++) {
    double f = (j - 1.0 + a) / (j + 1.0);
    r[j + 1] = f * c * r[j];
    if (want_grad) {
      dra[j + 1] = f * c * dra[j] + c * r[j] / (j + 1.0);
      drc[j + 1] = f * c * drc[j] + f * r[j];
    }
  }
  // The severity itself underflows: r_j ~ c^j, so once j > log(DBL_MIN)/log(c)
  // (j ~ 1380 at c = 0.6) it drops into the subnormals. It never reaches exact zero
  // -- it sticks at the smallest denormal, since 4.94e-324 * c rounds back to
  // 4.94e-324 -- so from there on r_j is stored many orders of magnitude too LARGE
  // and injects spurious mass into the convolution. That is what makes the raw
  // recursion return a saturated finite log-pmf in the far tail rather than failing.
  // zi is the first subnormal index and rlast bounds the error in every severity at
  // or beyond it; the two certify or reject the result at the end.
  const double TINY = std::numeric_limits<double>::min();   // 2.225e-308
  int zi = y;  double rlast = 0.0;
  for (int j = 1; j <= y - 1; j++) {
    if (r[j] < TINY) { zi = j; rlast = (j > 1) ? r[j - 1] : 0.0; break; }
  }

  double logscale = lp0, abc = a * b * c;
  p[0] = 1.0;                                  // p~_0 = p_0 / exp(logscale)
  if (want_grad) { pa[0] = d0a; pb[0] = d0b; pc[0] = d0c; }
  p[1] = abc;
  if (want_grad) {
    pa[1] = abc * d0a + b * c; pb[1] = abc * d0b + a * c; pc[1] = abc * d0c + a * b;
  }

  const double HI = 1e250, LO = 1e-250, RS = 1e250, CAP = 1e300;
  double mx = 1.0;                                    // running max over p[0..k]
  for (int k = 1; k <= y - 1; k++) {
    double s = 0.0, sa = 0.0, sb = 0.0, sc = 0.0;
    for (int j = k; j >= 1; j--) {          // descending j == ascending r-index
      int rj = k + 1 - j;
      double term = j * r[rj] * p[j];
      s += term;
      double ta = 0.0, tb = 0.0, tc = 0.0;
      if (want_grad) {
        ta = j * (p[j] * dra[rj] + r[rj] * pa[j]);
        tb = j * (r[rj] * pb[j]);
        tc = j * (p[j] * drc[rj] + r[rj] * pc[j]);
        sa += ta; sb += tb; sc += tc;
      }
      if (use_trunc && rj >= 8 && s > 0.0 && term < eps * s &&
          (!want_grad || (std::fabs(ta) < eps * std::fabs(sa) &&
                          std::fabs(tb) < eps * std::fabs(sb) &&
                          std::fabs(tc) < eps * std::fabs(sc)))) break;
    }
    double km = k + 1.0;
    p[k + 1] = (abc * p[k] + s) / km;
    if (want_grad) {
      pa[k + 1] = (abc * pa[k] + b * c * p[k] + sa) / km;
      pb[k + 1] = (abc * pb[k] + a * c * p[k] + sb) / km;
      pc[k + 1] = (abc * pc[k] + a * b * p[k] + sc) / km;
    }
    double pk = p[k + 1];
    if (!std::isfinite(pk)) { lpy = -INFINITY; dla = dlb = dlc = 0.0; return false; }
    if (pk > mx) mx = pk;

    if (pk > HI || (pk > 0.0 && pk < LO)) {
      double f;
      if (pk > HI) {
        f = 1.0 / RS;                                 // scaling down is always safe
      } else {
        f = RS;                                       // cap so no entry reaches Inf
        if (mx * f > CAP) f = CAP / mx;
        // No headroom left: max_j p_j and p_{k+1} are already >300 orders apart, so
        // no common scale holds both. Fail rather than limp on into the denormals,
        // which would silently return a saturated finite log-pmf.
        if (!(f > 1.0)) { lpy = -INFINITY; dla = dlb = dlc = 0.0; return false; }
      }
      for (int j = 0; j <= k + 1; j++) {
        p[j] *= f;
        if (want_grad) { pa[j] *= f; pb[j] *= f; pc[j] *= f; }
      }
      mx *= f;
      logscale -= std::log(f);
    }
  }

  if (!(p[y] > 0.0) || !std::isfinite(p[y])) {        // far tail lost to underflow
    lpy = -INFINITY; dla = dlb = dlc = 0.0; return false;
  }
  // Certify against the severity underflow. The error in each corrupted term is
  // bounded by rlast * max_j p_j and there are at most y of them, so the test below
  // reduces to the scale-invariant condition p_y / max_j p_j > 1e8 * y * rlast --
  // "p_y lies within ~297 orders of magnitude of the mode". If it does not, p_y is
  // not resolvable here; refuse rather than return a saturated value.
  if (zi <= y - 1 && !(p[y] > 1e8 * (double)y * rlast * mx)) {
    lpy = -INFINITY; dla = dlb = dlc = 0.0; return false;
  }
  lpy = logscale + std::log(p[y]);
  if (want_grad) { dla = pa[y] / p[y]; dlb = pb[y] / p[y]; dlc = pc[y] / p[y]; }
  return true;
}

// ------------------------------------------------ saddlepoint, O(1) per obs ----
//
// With z = c e^t and u = 1 - z, the CGF derivatives are K_n(t) = b * E_n(z;a) for
// n >= 1, where E_n = -d^n/dt^n u^a. Each E_n is a finite sum of terms g * z^m *
// u^(a-k); differentiating in t maps (g,m,k) -> (m g, m, k) and (-(a-k) g, m+1, k+1),
// so E_1..E_5 and their a-partials are built exactly by iterating that map. Two
// identities fall out and are used below as consistency checks:
//     dE_n/dc |_t = E_{n+1} / c        (E_n depends on c and t only through z)
//     dt^/dc      = -1/c               (z^ solves a b z (1-z)^(a-1) = y, free of c)
// which together force the c-derivatives of the -1/2 log K_2 and correction terms
// to vanish, leaving d log p/dc = -a b (1-c)^(a-1) + y/c in closed form.
// Largest |corr| at which the second-order Daniels term may be added at all. Beyond
// this the asymptotic series has diverged and the term is noise, not a correction.
static const double GPIG_CORR_MAX = 0.5;

// d^n/dt^n of u^a is a sum of terms g * z^m * u^(a-k). Differentiating maps
// (g,m,k) -> (m g, m, k) and (-(a-k) g, m+1, k+1), so k - m is invariant: starting from
// (1,0,0) every term has k = m. There are therefore at most n+1 distinct terms, not 2^n,
// and they can be indexed by m alone -- no container, no allocation. This matters: the
// normalising constant calls this routine nquad * 7 times per observation, and the naive
// heap-allocating version made a normalised evaluation 36x more expensive than a raw one.
//   g'[m]  = m g[m] - (a-m+1) g[m-1]
//   ga'[m] = m ga[m] - (g[m-1] + (a-m+1) ga[m-1])
static void gpig_En(double z, double u, double a, int nmax, double* E, double* dEa)
{
  double g[8] = {1.0, 0, 0, 0, 0, 0, 0, 0}, ga[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  double lu = std::log(u);
  for (int n = 1; n <= nmax; n++) {
    for (int m = n; m >= 0; m--) {                    // descend so g[m-1] is still the old value
      double gm  = (m <= n - 1) ? g[m]  : 0.0;
      double gam = (m <= n - 1) ? ga[m] : 0.0;
      double gp  = (m >= 1) ? g[m - 1]  : 0.0;
      double gap = (m >= 1) ? ga[m - 1] : 0.0;
      double ak  = a - (m - 1);
      g[m]  = m * gm  - (m >= 1 ? ak * gp : 0.0);
      ga[m] = m * gam - (m >= 1 ? (gp + ak * gap) : 0.0);
    }
    double v = 0.0, dv = 0.0, zm = 1.0;
    for (int m = 0; m <= n; m++) {
      double w = zm * std::exp((a - m) * lu);
      v  += g[m] * w;
      dv += (ga[m] + g[m] * lu) * w;
      zm *= z;
    }
    E[n] = -v; dEa[n] = -dv;
  }
}

// Stable z = logistic(w), u = 1 - z.
static inline void gpig_zu(double w, double& z, double& u) {
  if (w >= 0.0) { double e = std::exp(-w); z = 1.0 / (1.0 + e); u = e / (1.0 + e); }
  else          { double e = std::exp(w);  z = e / (1.0 + e);  u = 1.0 / (1.0 + e); }
}

// Saddlepoint equation K'(t) = y  <=>  a b z (1-z)^(a-1) = y, monotone in z.
static double gpig_what(double y, double a, double b) {
  double lo = -745.0, hi = 745.0;
  double target = std::log(y) - std::log(a * b);
  for (int i = 0; i < 200; i++) {
    double w = 0.5 * (lo + hi), z, u;
    gpig_zu(w, z, u);
    double f = std::log(z) + (a - 1.0) * std::log(u) - target;
    if (f < 0.0) lo = w; else hi = w;
  }
  return 0.5 * (lo + hi);
}

static void gpig_saddle(int y, double a, double b, double c, int order, bool want_grad,
                        double& lpy, double& dla, double& dlb, double& dlc,
                        double* acorr = NULL)
{
  double w = gpig_what((double)y, a, b), z, u;
  gpig_zu(w, z, u);
  double lu = std::log(u);
  double E[7] = {0}, dEa[7] = {0};
  gpig_En(z, u, a, 5, E, dEa);   // E_5 is needed for the correction and its diagnostic

  double omc = 1.0 - c, lomc = std::log(omc), omca = std::pow(omc, a);
  double ua = std::exp(a * lu);
  double K = b * (omca - ua);
  double that = std::log(z) - std::log(c);
  double E1 = E[1], E2 = E[2], E3 = E[3], E4 = E[4], E5 = E[5];
  double K2 = b * E2;

  double lp = K - that * (double)y - 0.5 * std::log(2.0 * M_PI * K2);

  // second-order (Daniels) correction: 1 + kappa4/8 - 5 kappa3^2/24, with
  // kappa3 = K3/K2^{3/2}, kappa4 = K4/K2^2, so the bracket equals g/b.
  double gg = E4 / (8.0 * E2 * E2) - 5.0 * E3 * E3 / (24.0 * E2 * E2 * E2);
  double corr = gg / b;
  if (acorr != NULL) *acorr = std::fabs(corr);   // size of the leading neglected term

  // The Daniels term is the leading element of an ASYMPTOTIC series: it may only be
  // applied where it is small. Admitting a large corr (it reaches 1e8 as c -> 1, the
  // discrete-stable boundary) adds log1p(corr) ~ +20 to log p, producing "pmf" values
  // far above 1 and an unbounded pseudo-likelihood that an optimiser will chase
  // straight to the corner. Order 1 alone remains a valid sub-probability there.
  bool use2 = (order >= 2);
  double dgt = 0.0, dga = 0.0;
  if (use2 && !(std::fabs(corr) <= GPIG_CORR_MAX)) use2 = false;   // fall back to order 1
  if (use2) {
    double dg2 = -E4 / (4.0 * E2 * E2 * E2) + 5.0 * E3 * E3 / (8.0 * std::pow(E2, 4.0));
    double dg3 = -5.0 * E3 / (12.0 * E2 * E2 * E2);
    double dg4 = 1.0 / (8.0 * E2 * E2);
    dgt = dg2 * E3 + dg3 * E4 + dg4 * E5;                  // dg/dt
    dga = dg2 * dEa[2] + dg3 * dEa[3] + dg4 * dEa[4];      // dg/da at fixed t
    lp += std::log1p(corr);
  } else {
    corr = 0.0;
  }
  // Hard invariant: a log-pmf cannot exceed zero. The saddlepoint is unnormalised, so
  // a small positive value is tolerable rounding, but anything above 1e-8 means the
  // approximation has failed here. Reject it -- returning a "probability" above 1
  // hands the optimiser a likelihood it can drive to infinity.
  if (!(lp <= 1e-8)) { lpy = -INFINITY; dla = dlb = dlc = 0.0; return; }
  lpy = lp;
  if (!want_grad) return;

  // dt^/dtheta = -(d K'/dtheta) / K''
  double dt_b = -E1 / (b * E2), dt_a = -dEa[1] / E2, dt_c = -1.0 / c;

  // envelope theorem: d[K(t^) - t^ y]/dtheta = dK/dtheta at fixed t
  double Kb = omca - ua;
  double Ka = b * (omca * lomc - ua * lu);
  double Kc = b * (-a * std::pow(omc, a - 1.0) + E1 / c);

  double t2b = -0.5 * (E2 + b * E3 * dt_b) / (b * E2);
  double t2a = -0.5 * (b * dEa[2] + b * E3 * dt_a) / (b * E2);
  double t2c = -0.5 * (b * E3 / c + b * E3 * dt_c) / (b * E2);   // == 0

  double t3a = 0.0, t3b = 0.0, t3c = 0.0;
  if (use2) {
    double f = 1.0 / (1.0 + corr);
    t3b = f * (-gg / (b * b) + dgt * dt_b / b);
    t3a = f * ((dga + dgt * dt_a) / b);
    t3c = f * ((dgt / c + dgt * dt_c) / b);                     // == 0
  }
  dla = Ka + t2a + t3a;
  dlb = Kb + t2b + t3b;
  dlc = Kc + t2c + t3c;
}

// ------------------------------------------- normalising the saddlepoint pmf ---
//
// The saddlepoint pmf is unnormalised: sum_k phat(k) = 1 + O(1/b). Near the
// discrete-stable boundary (c -> 1) b is not large and the residual is substantial --
// measured at ~0.15 log units per observation on the Tycho MEASLES fit, i.e. hundreds
// of AIC units over the sample. Worse, an unnormalised pmf is not a likelihood, so the
// optimiser can exploit it. Dividing by S = sum_k phat(k) restores a proper pmf and
// removes the leading relative bias.
//
// S is computed as the integral int phat(y) dy, reparameterised by the saddlepoint
// variable t: since y = K'(t) and dy = K''(t) dt,
//     S = int_{-inf}^{T} exp(K(t) - t K'(t)) sqrt(K''(t) / 2pi) dt,   T = -log c.
// Parameterising the nodes by t rather than y means each node is a direct evaluation
// -- no root solve -- so the whole constant costs O(nquad) cheap evaluations. The
// substitution t = T - e^u maps (-inf, T) to the whole real line and clusters nodes
// geometrically towards the boundary, where the heavy tail lives. The integrand peaks
// at t = 0 (d/dt[K - tK'] = -t K''), i.e. at u = log T, so the grid is centred there.
// The integrand MUST be the very pmf the likelihood uses, including the conditional
// second-order term -- normalising by the integral of a different function reintroduces
// exactly the bias it is meant to remove (an order-1 integrand against an order-2 pmf
// was off by 7% at c = 0.999936).
static double gpig_F(double t, double a, double b, double c, int order)
{
  double z = c * std::exp(t);
  if (!(z > 0.0) || !(z < 1.0)) return 0.0;
  double u = 1.0 - z, lu = std::log(u);
  double E[7] = {0}, dEa[7] = {0};
  gpig_En(z, u, a, 5, E, dEa);
  double E2 = E[2], E3 = E[3], E4 = E[4];
  double K2 = b * E2;
  if (!(K2 > 0.0) || !std::isfinite(K2)) return 0.0;
  double omca = std::pow(1.0 - c, a), ua = std::exp(a * lu);
  double lF = b * (omca - ua) - t * (b * E[1]) + 0.5 * std::log(K2 / (2.0 * M_PI));
  if (order >= 2) {
    double corr = (E4 / (8.0 * E2 * E2) - 5.0 * E3 * E3 / (24.0 * E2 * E2 * E2)) / b;
    if (std::fabs(corr) <= GPIG_CORR_MAX) lF += std::log1p(corr);   // same gate as gpig_saddle
  }
  if (!std::isfinite(lF) || lF < -700.0) return 0.0;
  return std::exp(lF);
}

// Trapezoid on the u-grid at a given resolution.
static double gpig_Squad(double a, double b, double c, int N, int order, double T)
{
  double u0 = std::log(T), lo = u0 - 25.0, hi = u0 + 25.0;
  double du = (hi - lo) / (N - 1);
  double S = 0.0;
  for (int i = 0; i < N; i++) {
    double u = lo + i * du, eu = std::exp(u);
    double w = (i == 0 || i == N - 1) ? 0.5 : 1.0;
    S += w * gpig_F(T - eu, a, b, c, order) * eu;   // Jacobian dt = -e^u du
  }
  return S * du;
}

// Returns log S, or exactly 0.0 to mean "do not normalise" -- either because S is
// indistinguishable from 1, or because the quadrature could not be trusted.
//
// Two traps, both of which silently corrupt the likelihood if ignored:
//
//  (1) The integrand is a peak at t = 0 of width 1/sqrt(K''(0)), which NARROWS as b
//      grows. A fixed uniform grid therefore fails for large b: at b = 4e4 a 256-node
//      rule returns S = 0, and the likelihood is divided by zero. The rule is refined
//      until it converges, and if it will not converge we decline to normalise rather
//      than return a wrong constant.
//  (2) Refinement is bounded. Beyond NMAX the peak is unresolvable at acceptable cost;
//      that happens only for large b, where S = 1 + O(1/b) is close to 1 anyway, so
//      declining to normalise is the safe action. (An earlier version tried to PREDICT
//      this from corr(t=0) and skip early -- unsound: at (a,b,c) = (0.609,187,0.999936)
//      corr(0) is below 1e-6 while S - 1 = 0.034. Convergence is tested, not predicted.)
static double gpig_logS(double a, double b, double c, int N, int order)
{
  double T = -std::log(c);                            // c < 1 by the clamp
  if (!(T > 0.0) || !std::isfinite(T)) return 0.0;

  const int NMAX = 1 << 13;
  double S = gpig_Squad(a, b, c, N, order, T);
  for (int it = 0; it < 8; it++) {
    int N2 = 2 * N;
    if (N2 > NMAX) return 0.0;                            // cannot resolve: leave as-is
    double S2 = gpig_Squad(a, b, c, N2, order, T);
    if (std::isfinite(S2) && S2 > 0.0 && std::fabs(S2 - S) <= 1e-5 * std::fabs(S2)) {
      S = S2;
      if (!(S > 0.0) || !std::isfinite(S)) return 0.0;
      return std::log(S);
    }
    S = S2; N = N2;
  }
  return 0.0;                                             // never converged: do not normalise
}

// Saddlepoint, normalised. log S is a smooth O(1) function of (a,b,c), so its gradient
// is taken by central differences on the quadrature -- far less error-prone than
// differentiating under the integral sign, and accurate to ~1e-8, well inside what the
// BHHH Hessian needs. Perturbations of a and c are taken in the interior of (0,1).
static void gpig_saddle_norm(int y, double a, double b, double c, int order, bool want_grad,
                             double& lpy, double& dla, double& dlb, double& dlc, double* acorr)
{
  gpig_saddle(y, a, b, c, order, want_grad, lpy, dla, dlb, dlc, acorr);
  if (!gpig_normalize || !std::isfinite(lpy)) return;
  double lS = gpig_logS(a, b, c, gpig_nquad, order);
  if (lS == 0.0) return;                          // quadrature failed: leave as-is
  lpy -= lS;
  if (!std::isfinite(lpy)) { lpy = -INFINITY; dla = dlb = dlc = 0.0; return; }
  if (!want_grad) return;
  const double h = 1e-5;
  double ha = h * a * (1.0 - a), hb = h * b, hc = h * c * (1.0 - c);
  dla -= (gpig_logS(a + ha, b, c, gpig_nquad, order) - gpig_logS(a - ha, b, c, gpig_nquad, order)) / (2.0 * ha);
  dlb -= (gpig_logS(a, b + hb, c, gpig_nquad, order) - gpig_logS(a, b - hb, c, gpig_nquad, order)) / (2.0 * hb);
  dlc -= (gpig_logS(a, b, c + hc, gpig_nquad, order) - gpig_logS(a, b, c - hc, gpig_nquad, order)) / (2.0 * hc);
}

// Test hook: the normalising constant under the current nquad.
// [[Rcpp::export]]
double gpig_logS_cpp(double a, double b, double c) { return gpig_logS(a, b, c, gpig_nquad, gpig_order); }

// --------------------------------------------------------------- dispatcher ----
static void gpig_logp_grad(int y, double a, double b, double c, bool want_grad,
                           double& lpy, double& dla, double& dlb, double& dlc)
{
  int m = gpig_method;
  if (y > 0 && m == GPIG_SADDLE) {
    gpig_saddle_norm(y, a, b, c, gpig_order, want_grad, lpy, dla, dlb, dlc, NULL);
    return;
  }
  if (y > 0 && m == GPIG_HYBRID && y > gpig_yswitch) {
    // Trust the saddlepoint only where its own asymptotic series has converged;
    // |corr| is the leading neglected term and is small exactly when b is large.
    double acorr = 0.0, slp, sda, sdb, sdc;
    gpig_saddle_norm(y, a, b, c, gpig_order, want_grad, slp, sda, sdb, sdc, &acorr);
    if (acorr <= gpig_sptol) {
      lpy = slp; dla = sda; dlb = sdb; dlc = sdc;
      return;
    }
    // Series not converged: prefer the exact recursion -- but only where it can
    // actually finish. The recursion costs O(y^2) PER OBSERVATION, so an unbounded
    // fallback at y ~ 1e5 never returns. An optimiser routinely transits (small b,
    // large y) regions where |corr| is large, so this bound is load-bearing, not a
    // formality. Above ymax we keep the saddlepoint: it is the only computable
    // option, and such regions are transient -- at a converged GAM b grows with the
    // fitted mean, which drives |corr| back below sptol.
    if (y <= gpig_ymax &&
        gpig_recur(y, a, b, c, want_grad, lpy, dla, dlb, dlc, false, gpig_eps)) return;
    lpy = slp; dla = sda; dlb = sdb; dlc = sdc;
    return;
  }
  if (m == GPIG_LEGACY) {
    double py, dpa, dpb, dpc;
    gpig_derivs(y, a, b, c, py, dpa, dpb, dpc);
    lpy = std::log(py);
    if (want_grad) { dla = dpa / py; dlb = dpb / py; dlc = dpc / py; }
    return;
  }
  gpig_recur(y, a, b, c, want_grad, lpy, dla, dlb, dlc,
             m == GPIG_TRUNC, gpig_eps);
}

// Test/validation hook: log pmf and d log pmf / d(a,b,c) under the current method.
// [[Rcpp::export]]
Rcpp::NumericVector gpig_logp_cpp(int y, double a, double b, double c) {
  double lp, da, db, dc, acorr = NA_REAL;
  gpig_logp_grad(y, a, b, c, true, lp, da, db, dc);
  if (y > 0) { double l2, a2, b2, c2; gpig_saddle(y, a, b, c, 2, false, l2, a2, b2, c2, &acorr); }
  (void)0;
  return Rcpp::NumericVector::create(Rcpp::Named("logp") = lp, Rcpp::Named("da") = da,
                                     Rcpp::Named("db") = db, Rcpp::Named("dc") = dc,
                                     Rcpp::Named("corr") = acorr);
}

// Boundary guard. The tail exponent a and down-weighting c enter the
// likelihood through terms like (1-c)^a and, in the mean parameterisation,
// b = mu(1-c)^(1-a)/(a c). As a -> 0 this b diverges (O(1/a)); at astronomically
// large b the pmf recursion overflows and returns a spurious (too-low) NLL,
// which the optimiser then chases to the degenerate corner a ~ 0. Clamping a
// and c to a safe interior keeps b finite and the recursion accurate. The
// bounds are wide enough not to bind on any well-identified fit: a = 1e-3 is
// already an extremely heavy tail (the a -> 0 limit is a proper distribution
// whose pmf differs from a = 1e-3 only in the ~4th decimal), and c = 1 - 1e-6
// is indistinguishable from the discrete-stable boundary.
static const double GPIG_A_LO = 1e-3, GPIG_A_HI = 1.0 - 1e-6;
static const double GPIG_C_LO = 1e-6, GPIG_C_HI = 1.0 - 1e-6;
static inline double clamp01(double x, double lo, double hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

static inline void eta_to_abc(int param, double e1, double e2, double e3,
                              double& a, double& b, double& c)
{
  if (param == 0) {                          // native (logit a, log b, logit c)
    a = clamp01(1.0 / (1.0 + std::exp(-e1)), GPIG_A_LO, GPIG_A_HI);
    b = std::exp(e2);
    c = clamp01(1.0 / (1.0 + std::exp(-e3)), GPIG_C_LO, GPIG_C_HI);
  } else {                                   // mean (log mu, logit a, logit c)
    double mu = std::exp(e1);
    a = clamp01(1.0 / (1.0 + std::exp(-e2)), GPIG_A_LO, GPIG_A_HI);
    c = clamp01(1.0 / (1.0 + std::exp(-e3)), GPIG_C_LO, GPIG_C_HI);
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
    double lpy, da, db, dc;
    gpig_logp_grad(y, a, b, c, false, lpy, da, db, dc);
    nllh += -lpy;
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
    double lpy, dla, dlb, dlc;
    gpig_logp_grad(y, a, b, c, true, lpy, dla, dlb, dlc);
    double ge[3];
    map_grad(param, a, b, c, dla, dlb, dlc, ge);
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
    double lpy, da, db, dc, logf;
    gpig_logp_grad(y, a, b, c, false, lpy, da, db, dc);
    if (y == 0) {
      logf = std::log(pi + (1.0 - pi) * std::exp(lpy));   // exp(lpy) -> 0 safely
    } else {
      logf = std::log1p(-pi) + lpy;
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

    double lpy, dla, dlb, dlc;
    gpig_logp_grad(y, a, b, c, true, lpy, dla, dlb, dlc);
    double g[4];

    if (y > 0) {
      double ge[3];
      map_grad(param, a, b, c, dla, dlb, dlc, ge);
      g[0] = -ge[0]; g[1] = -ge[1]; g[2] = -ge[2];
      g[3] = pi;                                        // -(score_pi = -pi)
    } else {
      double p0 = std::exp(lpy);                        // gpig pmf at 0
      double D = pi + (1.0 - pi) * p0;
      // d log D / d(a,b,c) = (1-pi) dp0/d(a,b,c) / D, then map to eta.
      // dp0/dtheta = p0 * dlog p0/dtheta, so f * dp0/dtheta = f * p0 * dl.
      double f = (1.0 - pi) * p0 / D;
      double ge[3];
      map_grad(param, a, b, c, f * dla, f * dlb, f * dlc, ge);
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
