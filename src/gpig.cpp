// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>

// Generalised Poisson-inverse Gaussian (GPIG) of Zhu & Joe (2009,
// Stat. Prob. Lett. 79:1695-1703) and its zero-inflated variant (ZIGPIG) for
// the native GAM fitter (egpd()) and the standalone d/p/q/r functions.
//
// GPIG(a, b, c) is defined by its pgf
//   G(s) = exp{ b[(1-c)^a - (1-cs)^a] },   0 < a <= 1, b > 0, 0 <= c <= 1,
// with Poisson (a=1), PIG (a=1/2) and discrete-stable (c=1) as special cases.
// Only c=1, a<1 is power-law; c<1 exponentially tempers the tail.
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
// The score is the exact derivative of the log pmf actually evaluated (the
// recursion, the tapered saddlepoint, or the hybrid's blend of the two); the
// Hessian is the empirical-Fisher / BHHH approximation (outer product of the
// per-observation NLL gradient), matching the estimator used by Zhu & Joe for
// standard errors. It is positive semi-definite; the
// fitter perturbs it when required.
//
// d12 returns per-observation derivatives of the NEGATIVE log-likelihood wrt
// the linear predictors, packed to match .gH2 / .gH3:
//   npar=3: cols [g1,g2,g3, H11,H12,H13,H22,H23,H33]                 (nobs x 9)
//   npar=4: cols [g1,g2,g3,g4, H11,H12,H13,H14,H22,H23,H24,H33,H34,H44] (x 14)

// ------------------------------------------------------------- fit methods ----
//
// The Zhu & Joe forward recursion costs O(y^2) and is carried in raw probability
// space, seeded by p_0 = exp(b[(1-c)^a - 1]). Both properties bite at large counts:
// the cost grows with the square of the largest count, and p_0 underflows to zero
// once b[(1-c)^a - 1] < -745 (mean of order 1e3), after which every p_k is
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
//                   saddlepoint where its own series has converged and the exact
//                   recursion everywhere else.
//
// ONE LADDER PER PARAMETER ROW. The recursion for p_y passes through p_0..p_{y-1},
// so every count of a row with the same (a,b,c) can be read off a single ladder.
// The exact routes therefore group observations by their parameter row, extend one
// ladder to the largest count of the row that needs it, and read the rest off it:
// an intercept-only fit costs O(max y^2) per evaluation rather than O(sum y^2).
// Ladders are also kept in a small LRU cache across calls (they depend on (a,b,c)
// alone), so d0 and d12 at the same parameters, or a sweep over y at fixed
// parameters, share the work.
//
// Accuracy of the saddlepoint is governed by the large parameter b (the compound-
// Poisson rate is b[1-(1-c)^a]), NOT by y: at b = 300 the second-order log-pmf error
// is ~1e-6 at every y, while at b = 5 it grows with y. The asymptotic series is
// therefore policed directly by the size of its own first correction term,
// corr = kappa4/8 - 5 kappa3^2/24 = g/b. "hybrid" uses the saddlepoint where
// |corr| <= sptol/2, the exact recursion where |corr| >= sptol, and a smooth
// (C^1) blend of the two in between, so the likelihood is continuous in the
// parameters across the switch. The exact branch is used up to ymax (default 1000);
// above it the saddlepoint, with its tapered second-order term, is kept even where
// |corr| > sptol. That bounds the exact work at O(ymax^2) per parameter row -- with
// covariates every observation is its own row, and an unbounded exact branch cost
// O(sum y^2) per evaluation (448 s against 0.24 s on a Tycho measles fit) -- at the
// price of an approximate likelihood above ymax where b is small. Small b is not
// transient: in the mean parameterisation b = mu (1-c)^(1-a) / (a c) -> 0 as c -> 1,
// the heavy-tailed regime GPIG is fitted for. ymax = INT_MAX (gpig_control(ymax = Inf))
// makes the hybrid exact wherever the saddlepoint has not converged. The likelihood is
// continuous in the parameters for every ymax: the branch depends on the count alone.
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
static double gpig_sptol   = 1e-3;   // |corr| above which "hybrid" uses the recursion
static int    gpig_ymax    = 1000;   // exact-branch bound (INT_MAX: no bound)
static int    gpig_nquad   = 256;    // minimum lattice truncation (legacy option name)
// Normalisation is opt-in and bounded; it removes total-mass, not shape, error.
static int    gpig_normalize = 0;

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


// log p_0 = b[(1-c)^a - 1] and its a- and b-derivatives. The direct form loses relative
// accuracy like eps / |a log(1-c)| as a log(1-c) -> 0 (small a, or small c), and b, of order
// 1/(a c) there, turns that into an error of order b eps in every log p_y (at a = 1e-20,
// b = 1e20 the whole pmf was doubled). Below |a log(1-c)| = 2^-8, where that error would
// exceed 6e-14, the bracket is formed as expm1(a log1p(-c)); above, the direct form is kept.
static const double GPIG_P0_SWITCH = 1.0 / 256.0;
static inline void gpig_p0(double a, double b, double c, double& lp0, double& d0a, double& d0b)
{
  double omc = 1.0 - c, l1 = std::log1p(-c), x = a * l1;
  if (std::fabs(x) < GPIG_P0_SWITCH) {
    double em = std::expm1(x);
    lp0 = b * em;
    d0a = b * (1.0 + em) * l1;
    d0b = em;
  } else {
    double omca = std::pow(omc, a);
    lp0 = b * (omca - 1.0);
    d0a = b * omca * std::log(omc);
    d0b = omca - 1.0;
  }
}


// ------------------------------------- exact recursion, carried in log space ---
//
// Identical mathematics to the raw recursion above, but p_k and dp_k/d(a,b,c) are
// held only up to a common positive factor: the recursion is linear and homogeneous in the
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
// It never returns a saturated finite value. The public pmf (gpig_pmf_cpp) passes
// `logpmf` and finishes such a ladder in log space; the fitting routes use the
// incremental ladder below, which performs the same arithmetic index by index, and
// "hybrid" falls through to the saddlepoint for an index the ladder cannot resolve.
static bool gpig_recur(int y, double a, double b, double c, bool want_grad,
                       double& lpy, double& dla, double& dlb, double& dlc,
                       bool use_trunc, double eps, std::vector<double>* logpmf = NULL)
{
  double omc = 1.0 - c;
  double lp0, d0a, d0b;                     // log p0 and d log p0 / d(a, b)
  gpig_p0(a, b, c, lp0, d0a, d0b);
  double d0c = -a * b * std::pow(omc, a - 1.0);
  if (logpmf) {
    logpmf->assign(y + 1, -INFINITY);
    (*logpmf)[0] = lp0;
  }

  if (y == 0) {
    lpy = lp0; dla = d0a; dlb = d0b; dlc = d0c;
    return true;
  }

  std::vector<double> p(y + 1), r(y + 1);
  std::vector<double> lr;
  if (logpmf) {
    lr.assign(y + 1, -INFINITY);
    lr[1] = std::log1p(-a) + std::log(c);
    for (int j = 1; j < y; j++)
      lr[j + 1] = lr[j] + std::log((j - 1.0 + a) / (j + 1.0)) + std::log(c);
  }
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
  double logabc = std::log(a) + std::log(b) + std::log(c);
  // Public probability ladders may extend hundreds of orders beyond the mode.
  // Finish that part of the same recursion in log space when a common scale
  // can no longer resolve it. Fitting's score recursion retains its refusal.
  auto finish_log = [&](int start) {
    for (int k = start; k < y; k++) {
      double v = logabc + (*logpmf)[k];
      for (int j = 1; j <= k; j++) {
        double term = std::log((double)j) + lr[k + 1 - j] + (*logpmf)[j];
        if (term > v) v = term + std::log1p(std::exp(v - term));
        else v += std::log1p(std::exp(term - v));
      }
      (*logpmf)[k + 1] = v - std::log(k + 1.0);
    }
    lpy = (*logpmf)[y];
    return std::isfinite(lpy);
  };
  p[0] = 1.0;                                  // p~_0 = p_0 / exp(logscale)
  if (want_grad) { pa[0] = d0a; pb[0] = d0b; pc[0] = d0c; }
  p[1] = abc;
  if (logpmf) (*logpmf)[1] = lp0 + logabc;
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
    if (logpmf && !(pk > 0.0)) return finish_log(k);
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
        if (!(f > 1.0)) {
          if (logpmf) return finish_log(k);
          lpy = -INFINITY; dla = dlb = dlc = 0.0; return false;
        }
      }
      for (int j = 0; j <= k + 1; j++) {
        p[j] *= f;
        if (want_grad) { pa[j] *= f; pb[j] *= f; pc[j] *= f; }
      }
      mx *= f;
      logscale -= std::log(f);
    }
    if (logpmf) {
      if (zi <= k && !(p[k + 1] > 1e8 * (k + 1.0) * rlast * mx))
        return finish_log(k);
      (*logpmf)[k + 1] = logscale + std::log(p[k + 1]);
    }
  }

  if (!(p[y] > 0.0) || !std::isfinite(p[y])) {        // far tail lost to underflow
    if (logpmf) return finish_log(y - 1);
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

// ------------------------------------------------- exact ladders, cached ----
//
// The fitting routes reach p_y through an incremental version of gpig_recur(): a
// ladder holds the recursion state (p~_k and dp~_k/d(a,b,c) at a common scale) and
// records log p_k and d log p_k/d(a,b,c) as each index is computed, so it can be
// extended on demand and read at any k it has passed. Index by index it performs the
// same arithmetic as gpig_recur() (scaled modes) or the original raw recursion
// (legacy mode), including the refusal to certify a p_k that the severity underflow
// has corrupted; a refused or unresolvable index is recorded as log p = -Inf.
enum GpigLadderMode { GPIG_LAD_LEGACY = 0, GPIG_LAD_SCALED = 1, GPIG_LAD_TRUNC = 2 };

struct GpigLadder {
  double a, b, c, eps;
  int mode;
  bool grad;
  std::vector<double> r, dra, drc;     // severity r_j (index 0 unused) and partials
  int zi;                              // first j with subnormal r_j (0 = none yet)
  double rlast;                        // r_{zi-1}
  std::vector<double> p, pa, pb, pc;   // recursion state at the current scale
  double logscale, mx, d0a, d0b, d0c;
  bool dead;                           // no common scale can continue the ladder
  int K;                               // p_0 .. p_K have been computed
  std::vector<double> lp, ga, gb, gc;  // recorded log p_k and its (a,b,c) gradient
  std::vector<unsigned char> ok;
  unsigned long stamp;
};

static void gpig_ladder_push(GpigLadder& L, bool good, double lp, double ga, double gb,
                             double gc)
{
  L.ok.push_back(good ? 1 : 0);
  L.lp.push_back(good ? lp : -INFINITY);
  if (L.grad) {
    L.ga.push_back(good ? ga : 0.0);
    L.gb.push_back(good ? gb : 0.0);
    L.gc.push_back(good ? gc : 0.0);
  }
}

static void gpig_ladder_init(GpigLadder& L, double a, double b, double c, int mode,
                             double eps, bool grad)
{
  const double TINY = std::numeric_limits<double>::min();
  L.a = a; L.b = b; L.c = c; L.mode = mode; L.eps = eps; L.grad = grad;
  L.r.assign(2, 0.0); L.dra.assign(2, 0.0); L.drc.assign(2, 0.0);
  L.r[1] = (1.0 - a) * c; L.dra[1] = -c; L.drc[1] = 1.0 - a;
  L.zi = 0; L.rlast = 0.0;
  if (L.r[1] < TINY) { L.zi = 1; L.rlast = 0.0; }
  double omc = 1.0 - c, lp0;
  if (mode == GPIG_LAD_LEGACY) {                      // raw semantics, exactly as before
    double omca = std::pow(omc, a);
    lp0 = b * (omca - 1.0);
    L.d0a = b * omca * std::log(omc);
    L.d0b = omca - 1.0;
  } else {
    gpig_p0(a, b, c, lp0, L.d0a, L.d0b);
  }
  L.d0c = -a * b * std::pow(omc, a - 1.0);
  L.p.clear(); L.pa.clear(); L.pb.clear(); L.pc.clear();
  L.lp.clear(); L.ga.clear(); L.gb.clear(); L.gc.clear(); L.ok.clear();
  L.dead = false; L.K = 0; L.mx = 1.0; L.logscale = 0.0;
  if (mode == GPIG_LAD_LEGACY) {
    double p0 = std::exp(lp0);
    L.p.push_back(p0);
    if (grad) { L.pa.push_back(L.d0a * p0); L.pb.push_back(L.d0b * p0); L.pc.push_back(L.d0c * p0); }
    // Raw semantics, including the p_0 underflow (log 0 = -Inf, 0/0 gradients).
    L.ok.push_back(1); L.lp.push_back(std::log(p0));
    if (grad) { L.ga.push_back(L.pa[0] / p0); L.gb.push_back(L.pb[0] / p0); L.gc.push_back(L.pc[0] / p0); }
  } else {
    L.p.push_back(1.0);                               // p~_0 = p_0 / exp(logscale)
    if (grad) { L.pa.push_back(L.d0a); L.pb.push_back(L.d0b); L.pc.push_back(L.d0c); }
    L.logscale = lp0;
    gpig_ladder_push(L, true, lp0, L.d0a, L.d0b, L.d0c);
  }
}

// Severity r_1..r_jmax (and partials when the ladder carries gradients).
static void gpig_ladder_severity(GpigLadder& L, int jmax)
{
  const double TINY = std::numeric_limits<double>::min();
  for (int j = (int)L.r.size() - 1; j < jmax; j++) {
    double f = (j - 1.0 + L.a) / (j + 1.0);
    double rn = f * L.c * L.r[j];
    L.r.push_back(rn);
    if (L.grad) {
      L.dra.push_back(f * L.c * L.dra[j] + L.c * L.r[j] / (j + 1.0));
      L.drc.push_back(f * L.c * L.drc[j] + f * L.r[j]);
    }
    if (L.zi == 0 && rn < TINY) { L.zi = j + 1; L.rlast = L.r[j]; }
  }
}

static void gpig_ladder_extend(GpigLadder& L, int Kt)
{
  if (Kt <= L.K) return;
  const double HI = 1e250, LO = 1e-250, RS = 1e250, CAP = 1e300;
  const double a = L.a, b = L.b, c = L.c, abc = a * b * c;
  const bool grad = L.grad, legacy = (L.mode == GPIG_LAD_LEGACY);
  const bool use_trunc = (L.mode == GPIG_LAD_TRUNC);
  const double eps = L.eps;
  gpig_ladder_severity(L, Kt);
  L.p.reserve(Kt + 1); L.lp.reserve(Kt + 1); L.ok.reserve(Kt + 1);
  if (grad) {
    L.pa.reserve(Kt + 1); L.pb.reserve(Kt + 1); L.pc.reserve(Kt + 1);
    L.ga.reserve(Kt + 1); L.gb.reserve(Kt + 1); L.gc.reserve(Kt + 1);
  }
  std::vector<double>& p = L.p; std::vector<double>& pa = L.pa;
  std::vector<double>& pb = L.pb; std::vector<double>& pc = L.pc;
  const std::vector<double>& r = L.r; const std::vector<double>& dra = L.dra;
  const std::vector<double>& drc = L.drc;
  while (L.K < Kt) {
    int k = L.K;                                        // computes index k + 1
    if ((k & 1023) == 1023) Rcpp::checkUserInterrupt();
    if (L.dead) {
      p.push_back(0.0);
      if (grad) { pa.push_back(0.0); pb.push_back(0.0); pc.push_back(0.0); }
      gpig_ladder_push(L, false, 0.0, 0.0, 0.0, 0.0);
      L.K++;
      continue;
    }
    double s = 0.0, sa = 0.0, sb = 0.0, sc = 0.0;
    if (legacy) {
      for (int j = 1; j <= k; j++) {
        int rj = k + 1 - j;
        s += j * r[rj] * p[j];
        if (grad) {
          sa += j * (p[j] * dra[rj] + r[rj] * pa[j]);
          sb += j * (r[rj] * pb[j]);
          sc += j * (p[j] * drc[rj] + r[rj] * pc[j]);
        }
      }
    } else {
      for (int j = k; j >= 1; j--) {                    // descending j == ascending r-index
        int rj = k + 1 - j;
        double term = j * r[rj] * p[j];
        s += term;
        double ta = 0.0, tb = 0.0, tc = 0.0;
        if (grad) {
          ta = j * (p[j] * dra[rj] + r[rj] * pa[j]);
          tb = j * (r[rj] * pb[j]);
          tc = j * (p[j] * drc[rj] + r[rj] * pc[j]);
          sa += ta; sb += tb; sc += tc;
        }
        if (use_trunc && rj >= 8 && s > 0.0 && term < eps * s &&
            (!grad || (std::fabs(ta) < eps * std::fabs(sa) &&
                       std::fabs(tb) < eps * std::fabs(sb) &&
                       std::fabs(tc) < eps * std::fabs(sc)))) break;
      }
    }
    double km = k + 1.0;
    p.push_back((abc * p[k] + s) / km);
    if (grad) {
      pa.push_back((abc * pa[k] + b * c * p[k] + sa) / km);
      pb.push_back((abc * pb[k] + a * c * p[k] + sb) / km);
      pc.push_back((abc * pc[k] + a * b * p[k] + sc) / km);
    }
    L.K++;
    double pk = p[k + 1];
    if (legacy) {
      L.ok.push_back(1); L.lp.push_back(std::log(pk));
      if (grad) { L.ga.push_back(pa[k + 1] / pk); L.gb.push_back(pb[k + 1] / pk);
                  L.gc.push_back(pc[k + 1] / pk); }
      continue;
    }
    if (!std::isfinite(pk)) { L.dead = true; gpig_ladder_push(L, false, 0, 0, 0, 0); continue; }
    if (pk > L.mx) L.mx = pk;
    if (pk > HI || (pk > 0.0 && pk < LO)) {
      double f;
      if (pk > HI) {
        f = 1.0 / RS;                                   // scaling down is always safe
      } else {
        f = RS;                                         // cap so no entry reaches Inf
        if (L.mx * f > CAP) f = CAP / L.mx;
        // No headroom: max_j p_j and p_{k+1} are already >300 orders apart.
        if (!(f > 1.0)) { L.dead = true; gpig_ladder_push(L, false, 0, 0, 0, 0); continue; }
      }
      for (int j = 0; j <= k + 1; j++) {
        p[j] *= f;
        if (grad) { pa[j] *= f; pb[j] *= f; pc[j] *= f; }
      }
      L.mx *= f;
      L.logscale -= std::log(f);
    }
    double pv = p[k + 1];
    bool good = pv > 0.0 && std::isfinite(pv);
    // Certify against the severity underflow exactly as gpig_recur() does at y = k+1.
    if (good && L.zi > 0 && L.zi <= k && !(pv > 1e8 * km * L.rlast * L.mx)) good = false;
    gpig_ladder_push(L, good, L.logscale + std::log(pv),
                     grad ? pa[k + 1] / pv : 0.0, grad ? pb[k + 1] / pv : 0.0,
                     grad ? pc[k + 1] / pv : 0.0);
  }
}

// LRU cache of ladders keyed by (a, b, c, mode[, eps]). A ladder that carries
// gradients also serves value-only requests. Total storage is bounded.
static std::vector<GpigLadder*> gpig_ladders;
static unsigned long gpig_tick = 0;
static const size_t GPIG_LADDER_SLOTS = 16;
static const double GPIG_LADDER_ENTRIES = 5e5;  // ~40 MB with gradients

static void gpig_ladder_trim(const GpigLadder* keep)
{
  for (;;) {
    double tot = 0.0;
    for (size_t i = 0; i < gpig_ladders.size(); i++) tot += gpig_ladders[i]->K + 1.0;
    if (gpig_ladders.size() <= 1 ||
        (tot <= GPIG_LADDER_ENTRIES && gpig_ladders.size() <= GPIG_LADDER_SLOTS)) return;
    size_t victim = gpig_ladders.size();
    for (size_t i = 0; i < gpig_ladders.size(); i++) {
      if (gpig_ladders[i] == keep) continue;
      if (victim == gpig_ladders.size() || gpig_ladders[i]->stamp < gpig_ladders[victim]->stamp)
        victim = i;
    }
    if (victim == gpig_ladders.size()) return;
    delete gpig_ladders[victim];
    gpig_ladders.erase(gpig_ladders.begin() + victim);
  }
}

static GpigLadder* gpig_ladder_get(double a, double b, double c, int mode, double eps,
                                   bool grad)
{
  gpig_tick++;
  for (size_t i = 0; i < gpig_ladders.size(); i++) {
    GpigLadder* L = gpig_ladders[i];
    if (L->a == a && L->b == b && L->c == c && L->mode == mode &&
        (mode != GPIG_LAD_TRUNC || L->eps == eps)) {
      if (grad && !L->grad) gpig_ladder_init(*L, a, b, c, mode, eps, true);
      L->stamp = gpig_tick;
      return L;
    }
  }
  GpigLadder* L = new GpigLadder();
  gpig_ladder_init(*L, a, b, c, mode, eps, grad);
  L->stamp = gpig_tick;
  gpig_ladders.push_back(L);
  gpig_ladder_trim(L);
  return L;
}

static inline int gpig_ladder_mode() {
  if (gpig_method == GPIG_LEGACY) return GPIG_LAD_LEGACY;
  if (gpig_method == GPIG_TRUNC) return GPIG_LAD_TRUNC;
  return GPIG_LAD_SCALED;
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
//
// The second-order (Daniels) term is the leading element of an ASYMPTOTIC series and
// may only be applied where it is small; beyond |corr| ~ 0.5 the series has diverged
// and the term is noise. It is TAPERED rather than switched off: corr is multiplied by
// a C^1 weight that is 1 for |corr| <= 0.25 and falls smoothly to 0 at |corr| = 0.5,
// so the saddlepoint log-pmf stays continuous (with a continuous score) as the
// parameters carry |corr| across the cut-off. A hard switch there made log p jump by
// ~0.7 nats.
static const double GPIG_TAPER_LO = 0.25, GPIG_TAPER_HI = 0.5;

// C^1 smoothstep S(t) = t^2 (3 - 2t) on [0,1], clamped outside; dS returns S'(t).
static inline double gpig_smoothstep(double t, double& dS) {
  if (!(t > 0.0)) { dS = 0.0; return 0.0; }
  if (t >= 1.0)   { dS = 0.0; return 1.0; }
  dS = 6.0 * t * (1.0 - t);
  return t * t * (3.0 - 2.0 * t);
}

// d^n/dt^n of u^a is a sum of terms g * z^m * u^(a-k). Differentiating maps
// (g,m,k) -> (m g, m, k) and (-(a-k) g, m+1, k+1), so k - m is invariant: starting from
// (1,0,0) every term has k = m. There are therefore at most n+1 distinct terms, not 2^n,
// and they can be indexed by m alone -- no container, no allocation. This matters: the
// normalising constant evaluates this routine at every lattice point.
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

// Saddlepoint equation K'(t) = y  <=>  a b z (1-z)^(a-1) = y with z = logistic(w).
// f(w) = log z + (a-1) log(1-z) - log(y/(ab)) is increasing (f' = 1 - a z > 0) and
// concave, so Newton's method, safeguarded by bisection on [-745, 745], converges in a
// handful of steps to the root the former 200-step bisection found.
static double gpig_what(double y, double a, double b) {
  double lo = -745.0, hi = 745.0;
  double target = std::log(y) - std::log(a * b);
  double w = (target <= 0.0 || a >= 1.0) ? target : target / (1.0 - a);
  if (!(w > lo && w < hi)) w = (w >= hi) ? hi - 1.0 : (w <= lo ? lo + 1.0 : 0.0);
  for (int i = 0; i < 400; i++) {
    double z, u;
    gpig_zu(w, z, u);
    double f = std::log(z) + (a - 1.0) * std::log(u) - target;
    if (f == 0.0) return w;
    if (f < 0.0) lo = w; else hi = w;
    double wn = w - f / (1.0 - a * z);
    if (!(wn > lo && wn < hi)) wn = 0.5 * (lo + hi);
    if (std::fabs(wn - w) <= 1e-15 * (1.0 + std::fabs(w)) || !(hi - lo > 0.0)) return wn;
    w = wn;
  }
  return w;
}

struct GpigSP {
  double lp, da, db, dc;     // log pmf and its (a,b,c) gradient
  double corr, ca, cb, cc;   // signed corr = kappa4/8 - 5 kappa3^2/24 and its gradient
};

// Returns false when the approximation is rejected (log p above zero); corr is set
// regardless, as a diagnostic.
static bool gpig_sp(int y, double a, double b, double c, int order, bool want_grad,
                    GpigSP& o)
{
  double w = gpig_what((double)y, a, b), z, u;
  gpig_zu(w, z, u);
  double lu = std::log(u);
  double E[7] = {0}, dEa[7] = {0};
  gpig_En(z, u, a, 5, E, dEa);   // E_5 is needed for the gradient of the correction

  double omc = 1.0 - c, lomc = std::log(omc), omca = std::pow(omc, a);
  double ua = std::exp(a * lu);
  double K = b * (omca - ua);
  double that = std::log(z) - std::log(c);
  double E1 = E[1], E2 = E[2], E3 = E[3], E4 = E[4], E5 = E[5];
  double K2 = b * E2;

  double lp = K - that * (double)y - 0.5 * std::log(2.0 * M_PI * K2);

  // Daniels correction 1 + kappa4/8 - 5 kappa3^2/24, with kappa3 = K3/K2^{3/2},
  // kappa4 = K4/K2^2, so the bracket equals g/b.
  double gg = E4 / (8.0 * E2 * E2) - 5.0 * E3 * E3 / (24.0 * E2 * E2 * E2);
  o.corr = gg / b;
  double dt_b = -E1 / (b * E2), dt_a = -dEa[1] / E2, dt_c = -1.0 / c;
  double dg2 = -E4 / (4.0 * E2 * E2 * E2) + 5.0 * E3 * E3 / (8.0 * std::pow(E2, 4.0));
  double dg3 = -5.0 * E3 / (12.0 * E2 * E2 * E2);
  double dg4 = 1.0 / (8.0 * E2 * E2);
  double dgt = dg2 * E3 + dg3 * E4 + dg4 * E5;                  // dg/dt
  double dga = dg2 * dEa[2] + dg3 * dEa[3] + dg4 * dEa[4];      // dg/da at fixed t
  o.cb = -gg / (b * b) + dgt * dt_b / b;
  o.ca = (dga + dgt * dt_a) / b;
  o.cc = (dgt / c + dgt * dt_c) / b;                            // == 0

  double fac = 0.0;
  if (order >= 2) {
    double x = std::fabs(o.corr), dS;
    double S = gpig_smoothstep((x - GPIG_TAPER_LO) / (GPIG_TAPER_HI - GPIG_TAPER_LO), dS);
    double tau = 1.0 - S, dtau = -dS / (GPIG_TAPER_HI - GPIG_TAPER_LO);
    double ce = o.corr * tau;                                   // |ce| <= 0.28 < 1
    lp += std::log1p(ce);
    fac = (tau + x * dtau) / (1.0 + ce);                        // d log1p(ce) / d corr
  }
  // Hard invariant: a log-pmf cannot exceed zero. The saddlepoint is unnormalised, so
  // a small positive value is tolerable rounding, but anything above 1e-8 means the
  // approximation has failed here. Reject it -- returning a "probability" above 1
  // hands the optimiser a likelihood it can drive to infinity.
  if (!(lp <= 1e-8)) { o.lp = -INFINITY; o.da = o.db = o.dc = 0.0; return false; }
  o.lp = lp;
  if (!want_grad) { o.da = o.db = o.dc = 0.0; return true; }

  // envelope theorem: d[K(t^) - t^ y]/dtheta = dK/dtheta at fixed t
  double Kb = omca - ua;
  double Ka = b * (omca * lomc - ua * lu);
  double Kc = b * (-a * std::pow(omc, a - 1.0) + E1 / c);

  double t2b = -0.5 * (E2 + b * E3 * dt_b) / (b * E2);
  double t2a = -0.5 * (b * dEa[2] + b * E3 * dt_a) / (b * E2);
  double t2c = -0.5 * (b * E3 / c + b * E3 * dt_c) / (b * E2);   // == 0

  o.da = Ka + t2a + fac * o.ca;
  o.db = Kb + t2b + fac * o.cb;
  o.dc = Kc + t2c + fac * o.cc;
  return true;
}


// --------------------------------------------------------------- dispatcher ----
struct GpigPt { double lp, ga, gb, gc; };

static inline void gpig_pt_fail(GpigPt& o) { o.lp = -INFINITY; o.ga = o.gb = o.gc = 0.0; }

static inline void gpig_pt_exact(GpigLadder*& L, int y, double a, double b, double c,
                                 bool want_grad, GpigPt& o, bool& ok)
{
  if (L == NULL) L = gpig_ladder_get(a, b, c, gpig_ladder_mode(), gpig_eps, want_grad);
  gpig_ladder_extend(*L, y);
  ok = L->ok[y] != 0;
  o.lp = L->lp[y];
  if (want_grad) { o.ga = L->ga[y]; o.gb = L->gb[y]; o.gc = L->gc[y]; }
  else o.ga = o.gb = o.gc = 0.0;
}

// The raw (unnormalised) log pmf at y under the current method, with its gradient.
// L is the row's ladder, acquired on first use; every y of a row shares it.
static void gpig_point(int y, double a, double b, double c, bool want_grad,
                       GpigLadder*& L, GpigPt& o)
{
  const int m = gpig_method;
  bool eok;
  if (m == GPIG_LEGACY || m == GPIG_RECURSION || m == GPIG_TRUNC || y == 0 ||
      (m == GPIG_HYBRID && y <= gpig_yswitch)) {
    gpig_pt_exact(L, y, a, b, c, want_grad, o, eok);
    return;                                // a refused index is already -Inf, 0
  }
  GpigSP s;
  bool spok = gpig_sp(y, a, b, c, gpig_order, want_grad, s);
  if (m == GPIG_SADDLE || y > gpig_ymax) {
    o.lp = s.lp; o.ga = s.da; o.gb = s.db; o.gc = s.dc;
    if (!spok) gpig_pt_fail(o);
    return;
  }
  // hybrid above yswitch: saddlepoint where |corr| <= sptol/2, exact recursion where
  // |corr| >= sptol, and a C^1 blend between, so log p is continuous in (a,b,c).
  double x = std::fabs(s.corr), h = 0.5 * gpig_sptol;
  if (spok && x <= h) { o.lp = s.lp; o.ga = s.da; o.gb = s.db; o.gc = s.dc; return; }
  GpigPt e;
  gpig_pt_exact(L, y, a, b, c, want_grad, e, eok);
  if (!eok) {
    // The exact ladder cannot resolve p_y (more than ~300 orders below the mode);
    // only the saddlepoint can, so it is kept even though |corr| is large.
    o.lp = s.lp; o.ga = s.da; o.gb = s.db; o.gc = s.dc;
    if (!spok) gpig_pt_fail(o);
    return;
  }
  double dS, w = gpig_smoothstep((x - h) / h, dS);
  if (!spok || w >= 1.0) { o = e; return; }
  o.lp = w * e.lp + (1.0 - w) * s.lp;
  if (!want_grad) { o.ga = o.gb = o.gc = 0.0; return; }
  double dwx = (s.corr >= 0.0 ? dS : -dS) / h, diff = e.lp - s.lp;
  o.ga = w * e.ga + (1.0 - w) * s.da + diff * dwx * s.ca;
  o.gb = w * e.gb + (1.0 - w) * s.db + diff * dwx * s.cb;
  o.gc = w * e.gc + (1.0 - w) * s.dc + diff * dwx * s.cc;
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
// The constant is a sum of the actual selected lattice masses, including zero and the
// hybrid's exact branch (read off the row's ladder). Its gradient is the same sum
// weighted by the per-point scores, d log S = sum_k phat_k d log phat_k / S, so the
// normalised score is exact and costs one pass rather than six extra sums. Results are
// cached per (a, b, c, options), so a row shared by many observations pays once.
struct GpigLogS { std::vector<double> key; double lS, ga, gb, gc; bool grad; };
static std::vector<GpigLogS> gpig_logS_cache;

static void gpig_logS_row(double a, double b, double c, bool want_grad,
                          double& lS, double& ga, double& gb, double& gc)
{
  ga = gb = gc = 0.0;
  if (gpig_method != GPIG_SADDLE && gpig_method != GPIG_HYBRID) { lS = 0.0; return; }
  if (!(c > 0 && c < 1)) Rcpp::stop("GPIG lattice normalization requires 0 < c < 1");
  const int N = gpig_nquad;
  std::vector<double> key = {a, b, c, (double)N, (double)gpig_order, (double)gpig_method,
    (double)gpig_yswitch, (double)gpig_ymax, gpig_sptol, gpig_eps};
  for (size_t i = 0; i < gpig_logS_cache.size(); i++) {
    const GpigLogS& e = gpig_logS_cache[i];
    if (e.key == key && (e.grad || !want_grad)) {
      lS = e.lS; ga = e.ga; gb = e.gb; gc = e.gc;
      return;
    }
  }
  GpigLadder* L = NULL;
  double sum = 0.0, sa = 0.0, sb = 0.0, sc = 0.0, previous = INFINITY, ratio = (1.0 + c) / 2.0;
  int small = 0;
  for (int k = 0; k <= 100000; k++) {
    GpigPt pt;
    gpig_point(k, a, b, c, want_grad, L, pt);
    if (!std::isfinite(pt.lp))
      Rcpp::stop("GPIG lattice normalization encountered unresolved mass");
    double pk = std::exp(pt.lp);
    sum += pk;
    if (want_grad) { sa += pk * pt.ga; sb += pk * pt.gb; sc += pk * pt.gc; }
    if (k > std::max(gpig_yswitch, N) && pk <= previous * ratio &&
        pk / (1.0 - ratio) < 1e-11 * sum) small++; else small = 0;
    if (small >= 32) {
      if (!(sum > 0.0 && std::isfinite(sum)))
        Rcpp::stop("GPIG lattice normalization failed");
      lS = std::log(sum);
      if (want_grad) { ga = sa / sum; gb = sb / sum; gc = sc / sum; }
      if (gpig_logS_cache.size() >= 64) gpig_logS_cache.erase(gpig_logS_cache.begin());
      GpigLogS e; e.key = key; e.lS = lS; e.ga = ga; e.gb = gb; e.gc = gc; e.grad = want_grad;
      gpig_logS_cache.push_back(e);
      return;
    }
    previous = pk;
  }
  Rcpp::stop("GPIG lattice normalization did not converge within 100000 counts");
}

static inline bool gpig_normalising() {
  return gpig_normalize && (gpig_method == GPIG_SADDLE || gpig_method == GPIG_HYBRID);
}

// Evaluate one parameter row at the counts ys[idx[0..n-1]], writing res[idx[i]].
static void gpig_eval_row(double a, double b, double c, const std::vector<int>& ys,
                          const int* idx, int n, bool want_grad, std::vector<GpigPt>& res)
{
  GpigLadder* L = NULL;
  for (int i = 0; i < n; i++) gpig_point(ys[idx[i]], a, b, c, want_grad, L, res[idx[i]]);
  if (gpig_normalising()) {
    double lS, ga, gb, gc;
    gpig_logS_row(a, b, c, want_grad, lS, ga, gb, gc);
    for (int i = 0; i < n; i++) {
      GpigPt& o = res[idx[i]];
      o.lp -= lS;
      if (want_grad) { o.ga -= ga; o.gb -= gb; o.gc -= gc; }
    }
  }
  gpig_ladder_trim(NULL);
}

// Test hook: the normalising constant under the current options.
// [[Rcpp::export]]
double gpig_logS_cpp(double a, double b, double c) {
  double lS, ga, gb, gc;
  gpig_logS_row(a, b, c, false, lS, ga, gb, gc);
  return lS;
}

// Test/validation hook: log pmf and d log pmf / d(a,b,c) under the current method, and
// |corr| of the order-2 saddlepoint at y (NA at y = 0).
// [[Rcpp::export]]
Rcpp::NumericVector gpig_logp_cpp(int y, double a, double b, double c) {
  std::vector<int> ys(1, y);
  std::vector<GpigPt> res(1);
  int idx0 = 0;
  gpig_eval_row(a, b, c, ys, &idx0, 1, true, res);
  double acorr = NA_REAL;
  if (y > 0) { GpigSP s; gpig_sp(y, a, b, c, 2, false, s); acorr = std::fabs(s.corr); }
  return Rcpp::NumericVector::create(Rcpp::Named("logp") = res[0].lp,
                                     Rcpp::Named("da") = res[0].ga,
                                     Rcpp::Named("db") = res[0].gb,
                                     Rcpp::Named("dc") = res[0].gc,
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
  // A clamped shape/down-weight has no derivative with respect to its link
  // predictor, including its indirect contribution through b on the mean scale.
  if (a == GPIG_A_LO || a == GPIG_A_HI) ge[param == 0 ? 0 : 1] = 0.0;
  if (c == GPIG_C_LO || c == GPIG_C_HI) ge[2] = 0.0;
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

// ------------------------------------------------------ rows of parameters ----
//
// Map every observation to (a, b, c), group identical rows, and evaluate each row once
// (one ladder, one normalising constant). NaN rows sort last and are never grouped
// with each other; their log pmf is NaN, which the fitter treats as a failed point.
static inline bool gpig_lt(double x, double y) {
  if (std::isnan(x)) return false;
  if (std::isnan(y)) return true;
  return x < y;
}

static void gpig_rows(int param, const arma::vec& e1, const arma::vec& e2,
                      const arma::vec& e3, const arma::vec& yvec, bool want_grad,
                      std::vector<double>& A, std::vector<double>& B,
                      std::vector<double>& C, std::vector<GpigPt>& res)
{
  int n = yvec.n_elem;
  A.resize(n); B.resize(n); C.resize(n); res.resize(n);
  std::vector<int> ys(n), idx(n);
  for (int j = 0; j < n; j++) {
    eta_to_abc(param, e1[j], e2[j], e3[j], A[j], B[j], C[j]);
    ys[j] = (int)(yvec[j] + 0.5);
    idx[j] = j;
  }
  std::sort(idx.begin(), idx.end(), [&](int i, int j) {
    if (gpig_lt(A[i], A[j])) return true;
    if (gpig_lt(A[j], A[i])) return false;
    if (gpig_lt(B[i], B[j])) return true;
    if (gpig_lt(B[j], B[i])) return false;
    return gpig_lt(C[i], C[j]);
  });
  int s = 0;
  while (s < n) {
    int i0 = idx[s], e = s + 1;
    while (e < n && A[idx[e]] == A[i0] && B[idx[e]] == B[i0] && C[idx[e]] == C[i0]) e++;
    gpig_eval_row(A[i0], B[i0], C[i0], ys, &idx[s], e - s, want_grad, res);
    s = e;
  }
}

static void gpig_eta3(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2,
                      const arma::mat& X3, const arma::uvec& dupid, int dcate,
                      const Rcpp::List& offsets, arma::vec& e1, arma::vec& e2, arma::vec& e3)
{
  e1 = X1 * Rcpp::as<arma::vec>(pars[0]);
  e2 = X2 * Rcpp::as<arma::vec>(pars[1]);
  e3 = X3 * Rcpp::as<arma::vec>(pars[2]);
  if (dcate == 1) { e1 = e1.elem(dupid); e2 = e2.elem(dupid); e3 = e3.elem(dupid); }
  arma::vec o0 = Rcpp::as<arma::vec>(offsets[0]); if (o0.n_elem > 0) e1 += o0;
  arma::vec o1 = Rcpp::as<arma::vec>(offsets[1]); if (o1.n_elem > 0) e2 += o1;
  arma::vec o2 = Rcpp::as<arma::vec>(offsets[2]); if (o2.n_elem > 0) e3 += o2;
}

// ----------------------------------------------------------- GPIG (3 par) -----

// Shared d0 for both parameterisations (param flag).
static double gpig_d0_impl(int param, const Rcpp::List& pars,
                           const arma::mat& X1, const arma::mat& X2,
                           const arma::mat& X3, arma::vec yvec,
                           const arma::uvec& dupid, int dcate,
                           const Rcpp::List& offsets)
{
  arma::vec e1, e2, e3;
  gpig_eta3(pars, X1, X2, X3, dupid, dcate, offsets, e1, e2, e3);
  std::vector<double> A, B, C; std::vector<GpigPt> res;
  gpig_rows(param, e1, e2, e3, yvec, false, A, B, C, res);
  double nllh = 0.0;
  for (size_t j = 0; j < res.size(); j++) nllh += -res[j].lp;
  return nllh;
}

static arma::mat gpig_d12_impl(int param, const Rcpp::List& pars,
                               const arma::mat& X1, const arma::mat& X2,
                               const arma::mat& X3, arma::vec yvec,
                               const arma::uvec& dupid, int dcate,
                               const Rcpp::List& offsets)
{
  arma::vec e1, e2, e3;
  gpig_eta3(pars, X1, X2, X3, dupid, dcate, offsets, e1, e2, e3);
  std::vector<double> A, B, C; std::vector<GpigPt> res;
  gpig_rows(param, e1, e2, e3, yvec, true, A, B, C, res);
  int nobs = yvec.n_elem;
  arma::mat out(nobs, 9);
  for (int j = 0; j < nobs; j++) {
    double ge[3];
    map_grad(param, A[j], B[j], C[j], res[j].ga, res[j].gb, res[j].gc, ge);
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
  arma::vec e1, e2, e3;
  gpig_eta3(pars, X1, X2, X3, dupid, dcate, offsets, e1, e2, e3);
  arma::vec e4 = X4 * Rcpp::as<arma::vec>(pars[3]);
  if (dcate == 1) e4 = e4.elem(dupid);
  { arma::vec o3 = Rcpp::as<arma::vec>(offsets[3]); if (o3.n_elem > 0) e4 += o3; }
  std::vector<double> A, B, C; std::vector<GpigPt> res;
  gpig_rows(param, e1, e2, e3, yvec, false, A, B, C, res);
  int nobs = yvec.size();
  double nllh = 0.0;
  for (int j = 0; j < nobs; j++) {
    double eta_pi = e4[j];
    double log_pi = -R::log1pexp(-eta_pi);
    double log1m_pi = -R::log1pexp(eta_pi);
    int y = (int)(yvec[j] + 0.5);
    double lpy = res[j].lp, logf;
    if (y == 0) {
      double b_term = log1m_pi + lpy;
      if (log_pi == b_term) {
        logf = log_pi + M_LN2;
      } else {
        double m = (log_pi > b_term) ? log_pi : b_term;
        logf = std::isfinite(m) ? (m + R::log1pexp(-std::fabs(log_pi - b_term))) : m;
      }
    } else {
      logf = log1m_pi + lpy;
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
  arma::vec e1, e2, e3;
  gpig_eta3(pars, X1, X2, X3, dupid, dcate, offsets, e1, e2, e3);
  arma::vec e4 = X4 * Rcpp::as<arma::vec>(pars[3]);
  if (dcate == 1) e4 = e4.elem(dupid);
  { arma::vec o3 = Rcpp::as<arma::vec>(offsets[3]); if (o3.n_elem > 0) e4 += o3; }
  std::vector<double> A, B, C; std::vector<GpigPt> res;
  gpig_rows(param, e1, e2, e3, yvec, true, A, B, C, res);
  int nobs = yvec.size();
  arma::mat out(nobs, 14);
  for (int j = 0; j < nobs; j++) {
    double a = A[j], b = B[j], c = C[j];
    double eta_pi = e4[j];
    double log_pi = -R::log1pexp(-eta_pi);
    double log1m_pi = -R::log1pexp(eta_pi);
    double pi = std::exp(log_pi);
    double one_m_pi = std::exp(log1m_pi);
    int y = (int)(yvec[j] + 0.5);
    double lpy = res[j].lp, dla = res[j].ga, dlb = res[j].gb, dlc = res[j].gc;
    double g[4];

    if (y > 0) {
      double ge[3];
      map_grad(param, a, b, c, dla, dlb, dlc, ge);
      g[0] = -ge[0]; g[1] = -ge[1]; g[2] = -ge[2];
      g[3] = pi;                                        // -(score_pi = -pi)
    } else {
      double delta = eta_pi - lpy;
      double log_w0 = -R::log1pexp(delta);
      double log_w1 = -R::log1pexp(-delta);
      double w0 = (delta >= 0.0) ? std::exp(log_w0) : (1.0 - std::exp(log_w1));
      double w1 = (delta >= 0.0) ? (1.0 - w0) : std::exp(log_w1);
      double ge[3];
      map_grad(param, a, b, c, w0 * dla, w0 * dlb, w0 * dlc, ge);
      double score_pi = (eta_pi >= 0.0) ? (one_m_pi - w0) : (w1 - pi);
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
  if (n < 0 || !(a > 0 && a <= 1 && b > 0 && c >= 0 && c <= 1) ||
      !std::isfinite(b)) Rcpp::stop("invalid GPIG parameters or support");
  Rcpp::NumericVector out(n + 1);
  Rcpp::NumericVector logs(n + 1);
  if (a == 1 || c == 0) {
    for (int k = 0; k <= n; k++) logs[k] = R::dpois(k, b*c, true);
  } else {
    double lp, da, db, dc;
    std::vector<double> ladder;
    if (!gpig_recur(n, a, b, c, false, lp, da, db, dc, false, gpig_eps, &ladder))
      Rcpp::stop("GPIG exact recursion cannot resolve this support");
    for (int k = 0; k <= n; k++) logs[k] = ladder[k];
  }
  for (int k = 0; k <= n; k++) out[k] = std::exp(logs[k]);
  out.attr("logpmf") = logs;
  return out;
}
