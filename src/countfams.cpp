// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>

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
// The Poisson-lognormal integral is evaluated by ADAPTIVE Gauss-Kronrod (G7-K15) quadrature of
// its log-concave integrand, between the points where the integrand has fallen e^-50 below
// its mode, with interval bisection until the estimated error is below 1e-13 of the integral.
// A FIXED Gauss-Hermite grid fails outright -- the Poisson kernel is a sharp spike when
// mu ~ 5e3, exactly the regime of the largest Tycho series, and fixed nodes then miss the
// mass and return zero. A mode-centred Gauss-Hermite rule fixes that, but not the opposite
// shape: at y = 0 with a large sigma the integrand is a wide Gaussian cut off sharply on
// one side by exp(-e^z), which 21 Hermite nodes resolved only to ~1e-3 in log P(Y = 0)
// (sum(dpoislnorm(0:2e5, 20, 3)) was 1 - 2e-4). Its score is exact, not a finite
// difference: d log p / d(muz, sigma) are posterior moments of z, which the same nodes
// already supply.
//
// The generalized Waring log pmf and score are assembled from cancellation-free first and
// second differences of lgamma and digamma (series for small increments, term-by-term
// differences of the asymptotic expansions otherwise) rather than from lgamma(a + rho) -
// lgamma(a) and the like, which lose ~1e-6 when a ~ 1e8 and everything as a -> Inf in the
// negative-binomial limit. The GAM likelihood and dgwaring() share that code.
//
// Exact tails. The standalone survivors that cannot be reached by summing the pmf use
// mixture integrals with the same quadrature: the Poisson-lognormal survivor is
// E[P(Pois(e^Z) > q)], and the generalized Waring is the beta-negative-binomial
// Y | p ~ NB(k, p), p ~ Beta(rho, a), so P(Y > q) = E[I_{1-p}(q+1, k)]. Both integrands are
// log-concave, so the same mode-and-limits construction applies.
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
  if (std::exp(e1) < GP_MU_LO || std::exp(e1) > GP_MU_HI) g[0] = 0.0;
  double rawlam = 1.0 / (1.0 + std::exp(-e2));
  if (rawlam < GP_LAM_LO || rawlam > GP_LAM_HI) g[1] = 0.0;
}

// ------------------------------------------ cancellation-free log-gamma differences ---
// The generalized Waring needs lgamma and digamma at arguments that may be ~1e15 or more
// times their differences (in its negative-binomial and Poisson limits rho, a and k grow
// without bound), so its log pmf and score are assembled from differences that are evaluated
// directly, never as differences of large lgamma or digamma values:
//
//   cf_digamma_diff(z, d)   psi(z + d) - psi(z)
//   cf_lgamma_dd(z, s, t)   lgamma(z+s+t) - lgamma(z+s) - lgamma(z+t) + lgamma(z)   (>= 0)
//   cf_lnbcoef(q, y)        lgamma(q + y) - lgamma(q) - lgamma(y + 1)
//
// The first two are Taylor series when the (smaller) increment is below 1e-4 z; otherwise the
// recurrences psi(x) = psi(x+1) - 1/x and lgamma(x) = lgamma(x+1) - log(x) lift z to >= 10
// (adding only positive terms), and the asymptotic expansions are differenced term by term,
// in forms whose terms are all of the size of the result. cf_lnbcoef is
// -log(q + y) - lbeta(q, y + 1); R's lbeta is assembled in the same way.
static inline double cf_stirling_corr(double x) {          // lgamma(x) - Stirling, x >= 10
  double r = 1.0 / x, r2 = r * r;
  return r * (1.0 / 12 - r2 * (1.0 / 360 - r2 * (1.0 / 1260 - r2 * (1.0 / 1680 -
         r2 * (1.0 / 1188 - r2 * (691.0 / 360360 - r2 * (1.0 / 156 - r2 * 3617.0 / 122400)))))));
}
// With S(x) = sum_n B_2n / (2n x^2n), so that psi(x) = log(x) - 1/(2x) - S(x) for x >= 10,
// S(z + d) - S(z) term by term as B_2n/(2n z^2n) expm1(-2n log1p(d/z)), which keeps its
// relative accuracy however small d is.
static inline double cf_digamma_dS(double z, double d) {
  static const double c[8] = {1.0 / 12, -1.0 / 120, 1.0 / 252, -1.0 / 240, 1.0 / 132,
                              -691.0 / 32760, 1.0 / 12, -3617.0 / 8160};
  const double l1 = std::log1p(d / z), r2 = 1.0 / (z * z);
  double zp = 1.0, out = 0.0;
  for (int n = 1; n <= 8; n++) { zp *= r2; out += c[n - 1] * zp * std::expm1(-2.0 * n * l1); }
  return out;
}
static double cf_digamma_diff(double z, double d) {
  if (d == 0.0) return 0.0;
  double acc = 0.0;
  if (!(d < 1e-4 * z))
    for (; z < 10.0; z += 1.0) acc += (d / (z + d)) / z;  // 1/z - 1/(z + d) > 0
  if (d < 1e-4 * z)
    return acc + d * (R::trigamma(z) + d * (0.5 * R::tetragamma(z) +
                 d * (R::pentagamma(z) / 6.0 + d * R::psigamma(z, 4.0) / 24.0)));
  return acc + std::log1p(d / z) + (d / (z + d)) / (2.0 * z) - cf_digamma_dS(z, d);
}
// (w + s) log1p(s/w) - s for 0 <= s <= w, without cancellation.
static inline double cf_bd0(double w, double s) {
  double r = s / w;
  return s * r + (w + s) * R::log1pmx(r);
}
static double cf_lgamma_dd(double z, double s, double t) {
  if (s > t) std::swap(s, t);
  if (s == 0.0) return 0.0;
  double acc = 0.0;
  if (!(s < 1e-4 * z))
    for (; z < 10.0; z += 1.0) acc += std::log1p((s / z) * (t / (z + s + t)));
  if (s < 1e-4 * z) {
    // sum_m s^m / m! [psi^(m-1)(z + t) - psi^(m-1)(z)]; the leading difference is the only
    // one that needs care, the later ones being smaller by s/z relative to it.
    double zt = z + t;
    return acc + s * (cf_digamma_diff(z, t) + s * (0.5 * (R::trigamma(zt) - R::trigamma(z)) +
           s * ((R::tetragamma(zt) - R::tetragamma(z)) / 6.0 +
           s * (R::pentagamma(zt) - R::pentagamma(z)) / 24.0)));
  }
  // z >= 10 and t >= s >= 1e-4 z: Stirling's series term by term. Its linear and constant
  // terms cancel exactly, and (x - 1/2) log(x) leaves the integrals of 1/x (I1) and 1/(2x^2)
  // (I2) over the square [0, s] x [0, t] shifted by z, both positive. I1 is written so
  // that its terms are of its own size: through cf_bd0 when s <= z, else directly.
  double I1 = (s <= z)
    ? s * std::log1p(t / z) + cf_bd0(z + t, s) - cf_bd0(z, s)
    : (z + t) * std::log1p(s / (z + t)) + (z + s) * std::log1p(t / (z + s)) -
      z * std::log1p((s + t) / z);
  double I2 = 0.5 * std::log1p((s / z) * (t / (z + s + t)));
  return acc + I1 + I2 + ((cf_stirling_corr(z + s + t) - cf_stirling_corr(z + t)) -
                          (cf_stirling_corr(z + s) - cf_stirling_corr(z)));
}
static inline double cf_lnbcoef(double q, double y) {
  return (y == 0.0) ? 0.0 : -std::log(q + y) - R::lbeta(q, y + 1.0);
}

// Generalized Waring at natural (a, k, rho), with c = a + k + rho:
//   log p(0) = lgamma(a+rho) + lgamma(k+rho) - lgamma(rho) - lgamma(c) = -DD(rho; a, k)
//   log p(y) = log p(0) + [lgamma(a+y) - lgamma(a)] + [lgamma(k+y) - lgamma(k)]
//              - [lgamma(c+y) - lgamma(c)] - lgamma(y+1)
//            = log p(0) - DD(hi; y, lo + rho) + lnbcoef(lo, y),  hi/lo = max/min(a, k).
// Pairing the larger of a and k with c differences exactly the two terms that grow without
// bound in the negative-binomial limit (rho -> Inf with a = mu (rho-1)/k), and the remaining
// negative-binomial coefficient is bounded there; for y >> a, k, rho every term is
// O(log y) rather than O(y log y). a is formed identically by the GAM likelihood and the
// standalone dgwaring(), so the two agree to the last bit inside the clamps.
static inline double gwaring_a(double mu, double k, double rho) {
  return mu * ((rho - 1.0) / k);
}
static inline double gwaring_logp0(double a, double k, double rho) {
  return -cf_lgamma_dd(rho, a, k);
}
static inline double gwaring_logpmf_one(double y, double a, double k, double rho,
                                        double l0) {
  if (y == 0.0) return l0;
  double hi = std::max(a, k), lo = std::min(a, k);
  return l0 - cf_lgamma_dd(hi, y, lo + rho) + cf_lnbcoef(lo, y);
}

// -------------------------------------------------- Generalized Waring -------------
// rho > 1 keeps the mean finite; rho > 2 is needed for a finite variance, but that is a
// property of the fit, not a constraint to impose -- a fitted rho in (1,2] is a meaningful
// statement (infinite variance) and is reported rather than forbidden.
//
// The clamps sit where the likelihood is numerically flat, so that they cannot distort a
// fit. rho -> Inf at fixed k is the negative-binomial limit NB(size = k, mu), and rho >> k
// -> Inf the Poisson limit: the optimiser escapes towards them on light-tailed data, and
// the generalized Waring differs from them by O(1/rho). At GW_RHO_HI = 1e15 the bound is
// therefore invisible for practical purposes, and the Poisson limit is approached to a squared
// coefficient of variation of the Poisson rate of 2/sqrt(mu rho), at k = sqrt(mu rho), well
// inside GW_K_HI. (The former bounds, rho <= 1e6 and k <= 1e8, cut the likelihood off
// before it had flattened: on Poisson(30) samples of 400 the maximum lies at rho ~ 3e6-6e6.)
// The other bounds are not plateaus: mu -> Inf, rho -> 1, k -> 0 or k -> Inf at fixed rho
// drive a or k to 0 or Inf and the likelihood of any positive count to -Inf; only mu -> 0,
// on all-zero data, is flat. Beyond a bound the likelihood is evaluated AT the bound, so it
// is flat there with a zero score, and predict() reports the clamped values (cf_bounds_cpp),
// at which dgwaring() reproduces the fitted likelihood. a itself is not clamped: inside the
// box a >= 1e-31 and the log pmf is exact for tiny a.
static const double GW_RHO_LO = 1.0 + 1e-6, GW_RHO_HI = 1e15;
static const double GW_K_LO   = 1e-8,  GW_K_HI  = 1e15;
static const double GW_MU_LO  = 1e-10, GW_MU_HI = 1e12;

static void gwaring_lpg(double y, double e1, double e2, double e3, double& lp, double* g)
{
  double mu  = cf_clamp(std::exp(e1), GW_MU_LO,  GW_MU_HI);
  double k   = cf_clamp(std::exp(e2), GW_K_LO,   GW_K_HI);
  double rho = cf_clamp(std::exp(e3), GW_RHO_LO, GW_RHO_HI);
  double a   = gwaring_a(mu, k, rho);

  lp = gwaring_logpmf_one(y, a, k, rho, gwaring_logp0(a, k, rho));
  if (!g) return;

  // d log p / d(a, k, rho) at fixed (a, k, rho), as digamma differences evaluated directly:
  //   la = Dpsi(a, y) - Dpsi(a+rho, k+y),  lk = Dpsi(k, y) - Dpsi(k+rho, a+y),
  //   lr = Dpsi(rho, lo) - Dpsi(hi+rho, lo+y),
  // and a = mu (rho-1)/k gives g = (a la, k lk - a la, rho lr + rho (mu/k) la). In the
  // negative-binomial and Poisson limits the two terms of each difference, and then the
  // terms of g, agree to O(1/rho), so there (arguments >= 10, where Stirling's series holds,
  // and |x| <= z/2) the leading logarithms are merged analytically: with
  // x = (y rho - a k)/(c + y),
  //   log1p(y/a) - log1p((k+y)/(a+rho)) = log1p(x/a),   likewise x/k for lk,
  //   log1p(lo/rho) - log1p((lo+y)/(hi+rho)) = log1p(-x/rho),
  // and the O(x) parts of g1 and g3 cancel exactly, leaving log1pmx terms.
  const double r = rho - 1.0, cy = a + k + rho + y, m = mu / k;
  const double mr = mu * r, mr_err = std::fma(mu, r, -mr);    // a k = mu (rho - 1), exactly
  const double x = (std::fma(y, rho, -mr) - mr_err) / cy;
  const double hi = std::max(a, k), lo = std::min(a, k), ax = std::fabs(x);
  const bool Am = a >= 10.0 && ax <= 0.5 * a, Km = k >= 10.0 && ax <= 0.5 * k,
             Rm = rho >= 10.0 && ax <= 0.5 * rho;
  // Dpsi(z, d) = log1p(d/z) + rest(z, d) for z >= 10
  auto rest = [](double z, double d) { return (d / (z + d)) / (2.0 * z) - cf_digamma_dS(z, d); };
  double ra = 0.0, rk = 0.0, rr = 0.0, la, lk, lr;
  if (Am) { ra = rest(a, y) - rest(a + rho, k + y); la = std::log1p(x / a) + ra; }
  else la = cf_digamma_diff(a, y) - cf_digamma_diff(a + rho, k + y);
  if (Km) { rk = rest(k, y) - rest(k + rho, a + y); lk = std::log1p(x / k) + rk; }
  else lk = cf_digamma_diff(k, y) - cf_digamma_diff(k + rho, a + y);
  if (Rm) { rr = rest(rho, lo) - rest(hi + rho, lo + y); lr = std::log1p(-x / rho) + rr; }
  else lr = cf_digamma_diff(rho, lo) - cf_digamma_diff(hi + rho, lo + y);

  g[0] = la * a;                                        // d/d eta1 (log mu)
  g[1] = (Am && Km)                                     // d/d eta2 (log k)
    ? k * R::log1pmx(x / k) - a * R::log1pmx(x / a) + (k * rk - a * ra)
    : lk * k - la * a;
  g[2] = (Am && Rm)                                     // d/d eta3 (log rho)
    ? x / r + rho * R::log1pmx(-x / rho) + rho * m * R::log1pmx(x / a) + rho * (rr + m * ra)
    : (lr + la * m) * rho;
  if (std::exp(e1) < GW_MU_LO || std::exp(e1) > GW_MU_HI) g[0] = 0.0;
  if (std::exp(e2) < GW_K_LO || std::exp(e2) > GW_K_HI) g[1] = 0.0;
  if (std::exp(e3) < GW_RHO_LO || std::exp(e3) > GW_RHO_HI) g[2] = 0.0;
}

// ------------------------------------------- adaptive Gauss-Kronrod quadrature ------
// G7-K15 on [a, b] with bisection of the worst interval (QUADPACK's error heuristic) until
// the total error estimate is below rtol times the integral. The integrand is supplied on
// the log scale and integrated as exp(g(z) - gref); the first two moments of z - c0 under
// the normalised integrand are returned on request (used by the Poisson-lognormal score).
static const double GK_XK[8] = {
  0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
  0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
  0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
  0.207784955007898467600689403773245, 0.0};
static const double GK_WK[8] = {
  0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
  0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
  0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
  0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
static const double GK_WG[4] = {
  0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
  0.381830050505118944950369775488975, 0.417959183673469387755102040816327};

struct GkSeg { double a, b, i0, i1, i2, err; };

template <class G>
static GkSeg gk_seg(G& g, double a, double b, double gref, double c0, bool mom)
{
  GkSeg o; o.a = a; o.b = b;
  double c = 0.5 * (a + b), h = 0.5 * (b - a);
  double fv[15], zv[15];
  zv[7] = c; fv[7] = std::exp(g(c) - gref);
  for (int j = 0; j < 7; j++) {
    double dx = h * GK_XK[j];
    zv[j] = c - dx;      fv[j] = std::exp(g(c - dx) - gref);
    zv[14 - j] = c + dx; fv[14 - j] = std::exp(g(c + dx) - gref);
  }
  double rk = GK_WK[7] * fv[7], rg = GK_WG[3] * fv[7], r1 = 0.0, r2 = 0.0;
  if (mom) { double d = c - c0; r1 = GK_WK[7] * fv[7] * d; r2 = r1 * d; }
  for (int j = 0; j < 7; j++) {
    double w = GK_WK[j], f1 = fv[j], f2 = fv[14 - j];
    rk += w * (f1 + f2);
    if (j % 2 == 1) rg += GK_WG[(j - 1) / 2] * (f1 + f2);
    if (mom) {
      double d1 = zv[j] - c0, d2 = zv[14 - j] - c0;
      r1 += w * (f1 * d1 + f2 * d2);
      r2 += w * (f1 * d1 * d1 + f2 * d2 * d2);
    }
  }
  double half = 0.5 * rk, asc = GK_WK[7] * std::fabs(fv[7] - half);
  for (int j = 0; j < 7; j++)
    asc += GK_WK[j] * (std::fabs(fv[j] - half) + std::fabs(fv[14 - j] - half));
  o.i0 = rk * h; o.i1 = r1 * h; o.i2 = r2 * h;
  double err = std::fabs((rk - rg) * h);
  asc *= h;
  if (asc != 0.0 && err != 0.0) err = asc * std::min(1.0, std::pow(200.0 * err / asc, 1.5));
  o.err = err;
  return o;
}

// Integral of exp(g - gref) over consecutive intervals given by `brk` (ascending).
template <class G>
static double gk_adapt(G& g, const std::vector<double>& brk, double gref, double c0,
                       double* m1, double* m2, double rtol = 1e-13, int maxseg = 400)
{
  bool mom = (m1 != NULL);
  std::vector<GkSeg> seg;
  for (size_t i = 0; i + 1 < brk.size(); i++)
    if (brk[i + 1] > brk[i]) seg.push_back(gk_seg(g, brk[i], brk[i + 1], gref, c0, mom));
  for (;;) {
    double I = 0.0, E = 0.0; size_t worst = 0;
    for (size_t i = 0; i < seg.size(); i++) {
      I += seg[i].i0; E += seg[i].err;
      if (seg[i].err > seg[worst].err) worst = i;
    }
    if (!(E > rtol * I) || (int)seg.size() >= maxseg || seg.empty()) {
      if (mom) {
        double s1 = 0.0, s2 = 0.0;
        for (size_t i = 0; i < seg.size(); i++) { s1 += seg[i].i1; s2 += seg[i].i2; }
        *m1 = s1 / I; *m2 = s2 / I;
      }
      return I;
    }
    GkSeg w = seg[worst];
    double mid = 0.5 * (w.a + w.b);
    if (!(mid > w.a && mid < w.b)) {                    // cannot split further
      seg[worst].err = 0.0;
      continue;
    }
    seg[worst] = gk_seg(g, w.a, mid, gref, c0, mom);
    seg.push_back(gk_seg(g, mid, w.b, gref, c0, mom));
  }
}

// For a concave log-integrand g with maximum gmax at zm, the points left and right of zm
// where g has fallen by T (located to a small fraction of `scale`).
template <class G>
static void lc_limits(G& g, double zm, double gmax, double scale, double T,
                      double& zl, double& zr)
{
  for (int side = -1; side <= 1; side += 2) {
    double inner = zm, d = scale, outer = zm + side * d;
    int it = 0;
    while (g(outer) > gmax - T && it < 200) { inner = outer; d *= 2.0; outer = zm + side * d; it++; }
    for (int b = 0; b < 60; b++) {
      double mid = 0.5 * (inner + outer);
      if (std::fabs(outer - inner) <= 1e-2 * scale) break;
      if (g(mid) > gmax - T) inner = mid; else outer = mid;
    }
    if (side < 0) zl = outer; else zr = outer;
  }
}

// Width 1/sqrt(-g''(zm)) of a concave log-integrand at its maximum, by a central
// difference whose step is adapted until the drop over it is moderate.
template <class G>
static double lc_width(G& g, double zm, double h0)
{
  double h = h0, g0 = g(zm);
  for (int i = 0; i < 60; i++) {
    double drop = g0 - 0.5 * (g(zm + h) + g(zm - h));
    if (drop > 1.0) { h *= 0.25; continue; }
    if (!(drop >= 1e-4)) { h *= 4.0; continue; }
    return h / std::sqrt(2.0 * drop);
  }
  return h0;
}

// Breakpoints zl < ... < zm < ... < zr spaced w, 2w, 4w, ... either side of the mode, so
// that the adaptive rule starts from intervals matched to the integrand's local scale
// (a single interval across a sharp shoulder can fool its error estimate).
static std::vector<double> lc_breaks(double zl, double zm, double zr, double w)
{
  std::vector<double> left, out;
  for (double d = w; zm - d > zl; d *= 2.0) left.push_back(zm - d);
  out.push_back(zl);
  for (int i = (int)left.size() - 1; i >= 0; i--) out.push_back(left[i]);
  out.push_back(zm);
  for (double d = w; zm + d < zr; d *= 2.0) out.push_back(zm + d);
  out.push_back(zr);
  return out;
}

// Maximum of a concave function by bracket expansion and golden-section search.
template <class G>
static double lc_argmax(G& g, double z0, double scale)
{
  double a = z0 - scale, b = z0, c = z0 + scale;
  double fa = g(a), fb = g(b), fc = g(c);
  int it = 0;
  while (!(fb >= fa && fb >= fc) && it < 200) {
    if (fa > fb) { c = b; fc = fb; b = a; fb = fa; a = b - 2.0 * (c - b); fa = g(a); }
    else         { a = b; fa = fb; b = c; fb = fc; c = b + 2.0 * (b - a); fc = g(c); }
    it++;
  }
  const double r = 0.6180339887498949;
  double x1 = c - r * (c - a), x2 = a + r * (c - a), f1 = g(x1), f2 = g(x2);
  for (int i = 0; i < 200 && (c - a) > 1e-10 * (1.0 + std::fabs(b)); i++) {
    if (f1 >= f2) { c = x2; x2 = x1; f2 = f1; x1 = c - r * (c - a); f1 = g(x1); }
    else          { a = x1; x1 = x2; f1 = f2; x2 = a + r * (c - a); f2 = g(x2); }
  }
  return (f1 >= f2) ? x1 : x2;
}

// -------------------------------------------------- Poisson-lognormal -------------
static const double PLN_MU_LO = 1e-10, PLN_MU_HI = 1e12;
static const double PLN_SG_LO = 1e-6,  PLN_SG_HI = 1e3;
static const double PLN_T = 50.0;                // integrand truncated e^-50 below its mode
static double pln_lcdf(double q, double mu, double sg, bool lower, double* score = NULL);

// For either a count or a tail event, differentiate the normal mixing law
// under the integral. m1 and m2 are conditional moments of z - muz.
static void pln_score(double s2, double m1, double m2, double* score) {
  score[0] = m1 / s2;
  score[1] = m2 / s2 - 1.0 - m1;  // muz = log(mu) - s2/2
}

// e^d - 1 - d without cancellation.
static inline double cf_em1d(double d) {
  if (std::fabs(d) < 1e-2) {
    double d2 = d * d;
    return d2 * (0.5 + d * (1.0 / 6 + d * (1.0 / 24 + d * (1.0 / 120 + d * (1.0 / 720 +
           d * (1.0 / 5040 + d / 40320))))));
  }
  return std::expm1(d) - d;
}

// log p(y) and optionally its log-parameter scores. The integrand is written about its
// mode z*, h(z* + d) - h(z*) = -e^{z*} (e^d - 1 - d) + d h'(z*) - d^2/(2 s^2), so it is
// evaluated to full relative accuracy even when h itself is of order y log y; and
// h(z*) - log y! is taken from the Poisson log pmf, which is computed without cancellation.
static double pln_lpmf(double y, double mu, double sg, double* score)
{
  double s2 = sg * sg;
  // P(Y>0) <= E[sqrt(Lambda)] = sqrt(mu)*exp(-s2/8). When small,
  // integrate the positive tail rather than trying to resolve log P(0)
  // from a quadrature result indistinguishable from unit total mass.
  if (y == 0.0 && 0.5 * std::log(mu) - s2 / 8.0 < -4.61) {
    double ls = pln_lcdf(0.0, mu, sg, false, score);
    double lp = std::log1p(-std::exp(ls));
    if (score) {
      double w = -std::exp(ls - lp);
      score[0] *= w; score[1] *= w;
    }
    return lp;
  }
  double muz = std::log(mu) - 0.5 * s2;
  auto hp = [&](double z) { return -std::exp(z) + y - (z - muz) / s2; };
  // h is strictly concave: bracket the root of h' and refine by safeguarded Newton.
  double z = (y > 0.0) ? std::max(muz, std::log(y)) : std::min(muz, 0.0);
  if (z > 700.0) z = 700.0;
  double lo = z, hi = z, st = 1.0;
  for (int i = 0; i < 2000 && hp(lo) <= 0.0; i++) { lo -= st; st *= 2.0; }
  st = 1.0;
  for (int i = 0; i < 2000 && hp(hi) >= 0.0; i++) { hi += st; st *= 2.0; }
  for (int it = 0; it < 200; it++) {
    double f = hp(z);
    if (f > 0.0) lo = z; else hi = z;
    double zn = z - f / (-std::exp(z) - 1.0 / s2);
    if (!(zn > lo && zn < hi)) zn = 0.5 * (lo + hi);
    if (std::fabs(zn - z) <= 1e-14 * (1.0 + std::fabs(z))) { z = zn; break; }
    z = zn;
  }
  double Ez = std::exp(z), hps = hp(z), dz = z - muz;
  auto g = [&](double d) {
    // At d>50, (1+d)/exp(d) is below working precision. Combining
    // exponents avoids multiplying an underflowed exp(z) by an infinity.
    double bend = d > 50.0 ? std::exp(z + d) : Ez * cf_em1d(d);
    return -bend + d * hps - d * d / (2.0 * s2);
  };
  double sc = 1.0 / std::sqrt(Ez + 1.0 / s2), dl, dr;
  lc_limits(g, 0.0, 0.0, sc, PLN_T, dl, dr);
  std::vector<double> brk = lc_breaks(dl, 0.0, dr, sc);
  double e1 = 0.0, e2 = 0.0;
  double I = gk_adapt(g, brk, 0.0, 0.0, score ? &e1 : NULL, score ? &e2 : NULL);
  if (score) pln_score(s2, dz + e1, dz * dz + 2.0 * dz * e1 + e2, score);
  return R::dpois(y, Ez, 1) - dz * dz / (2.0 * s2) - 0.5 * std::log(2.0 * M_PI) -
    std::log(sg) + std::log(I);
}

// log P(Y <= q) (lower) or log P(Y > q) for the Poisson-lognormal: the Poisson cdf in
// pgamma form, P(Pois(l) <= q) = P(Gamma(q + 1) > l), averaged over the lognormal.
static double pln_lcdf(double q, double mu, double sg, bool lower, double* score)
{
  double s2 = sg * sg;
  double muz = std::log(mu) - 0.5 * s2;
  auto g = [&](double z) {
    double d = z - muz;
    double lp;
    if (z < -700.0) {
      // The leading Poisson upper-tail term retains log(lambda) even
      // when lambda itself is unrepresentable; relative corrections are O(lambda).
      double ls = (q + 1.0) * z - R::lgammafn(q + 2.0);
      lp = lower ? std::log1p(-std::exp(ls)) : ls;
    } else {
      lp = R::pgamma(std::exp(z), q + 1.0, 1.0, lower ? 0 : 1, 1);
    }
    return lp - d * d / (2.0 * s2);
  };
  double z0 = lower ? std::min(muz, std::log(q + 1.0)) : std::max(muz, std::log(q + 1.0));
  double zm = lc_argmax(g, z0, std::max(sg, 1e-3));
  double gm = g(zm), zl, zr;
  double w = lc_width(g, zm, std::max(std::min(sg, 1.0), 1e-3));
  lc_limits(g, zm, gm, w, PLN_T, zl, zr);
  std::vector<double> brk = lc_breaks(zl, zm, zr, w);
  double m1, m2;
  double I = gk_adapt(g, brk, gm, muz, score ? &m1 : NULL, score ? &m2 : NULL);
  if (score) pln_score(s2, m1, m2, score);
  double out = gm + std::log(I) - 0.5 * std::log(2.0 * M_PI) - std::log(sg);
  return std::min(out, 0.0);
}

static void plnorm_lpg(double y, double e1, double e2, double& lp, double* g)
{
  double mu = cf_clamp(std::exp(e1), PLN_MU_LO, PLN_MU_HI);
  double sg = cf_clamp(std::exp(e2), PLN_SG_LO, PLN_SG_HI);
  lp = pln_lpmf(y, mu, sg, g);
  if (!g) return;
  if (std::exp(e1) < PLN_MU_LO || std::exp(e1) > PLN_MU_HI) g[0] = 0.0;
  if (std::exp(e2) < PLN_SG_LO || std::exp(e2) > PLN_SG_HI) g[1] = 0.0;
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
// Exact (unclamped) natural-parameter log pmfs, vectorised over non-negative integer y.
// Used by the standalone d/p/q/r wrappers (whereas the GAM likelihoods above enforce
// the cf_clamp bounds on eta).

static inline double gpois_logpmf_exact_one(double y, double mu, double lambda) {
  if (lambda == 0.0) return R::dpois(y, mu, 1);
  double th = mu * (1.0 - lambda);
  if (y == 0.0) return -th;
  double A = th + lambda * y;
  if (th > 0.0) return -std::log1p((lambda * y) / th) + R::dpois(y, A, 1);
  return std::log(mu) + std::log1p(-lambda) + (y - 1.0) * std::log(A) - A - R::lgammafn(y + 1.0);
}

static inline double plnorm_logpmf_exact_one(double y, double mu, double sg) {
  if (sg * (1.0 + mu + y) < 1e-7) {
    double corr = 0.5 * sg * sg * ((y - mu) * (y - mu) - y);
    return R::dpois(y, mu, 1) + corr;
  }
  return pln_lpmf(y, mu, sg, NULL);
}

// [[Rcpp::export]]
Rcpp::NumericVector gpois_logpmf_cpp(Rcpp::NumericVector y, double mu, double lambda) {
  int n = y.size();
  Rcpp::NumericVector out(n);
  for (int i = 0; i < n; i++) out[i] = gpois_logpmf_exact_one(y[i], mu, lambda);
  return out;
}
// [[Rcpp::export]]
Rcpp::NumericVector gwaring_logpmf_cpp(Rcpp::NumericVector y, double mu, double k, double rho) {
  int n = y.size();
  Rcpp::NumericVector out(n);
  double a = gwaring_a(mu, k, rho);
  double l0 = gwaring_logp0(a, k, rho);
  for (int i = 0; i < n; i++) out[i] = gwaring_logpmf_one(y[i], a, k, rho, l0);
  return out;
}
// [[Rcpp::export]]
Rcpp::NumericVector plnorm_logpmf_cpp(Rcpp::NumericVector y, double mu, double sigma) {
  int n = y.size();
  Rcpp::NumericVector out(n);
  for (int i = 0; i < n; i++) out[i] = plnorm_logpmf_exact_one(y[i], mu, sigma);
  return out;
}

// Poisson-lognormal log cdf (lower) or log survivor at integer q >= 0.
// [[Rcpp::export]]
Rcpp::NumericVector plnorm_lcdf_cpp(Rcpp::NumericVector q, double mu, double sigma, bool lower) {
  int n = q.size();
  Rcpp::NumericVector out(n);
  for (int i = 0; i < n; i++) out[i] = pln_lcdf(q[i], mu, sigma, lower);
  return out;
}

// Generalized Waring log cdf (lower) or log survivor at integer q >= 0, as the
// beta-negative-binomial mixture integral over w ~ Beta(rho, a):
//   P(Y > q) = E[ P(Beta(k, q+1) > w) ],  P(Y <= q) = E[ P(Beta(k, q+1) <= w) ],
// integrated over x = logit(w), where the integrand is log-concave. In the negative-binomial
// and Poisson limits a and rho are huge and the mixing density is a spike of width
// ~1/sqrt(a + rho), so its logarithm h(x) = rho log w + a log(1 - w) - lbeta(rho, a) is never
// formed from terms of size a or rho (which lost 1.5% of P(Y > q) at rho = 1e14). Points
// are parametrised by t = e^-x about a centre t_c (w_c = 1/(1+t_c), u_c = 1 - w_c); with
// E = e^d - 1, the drop h(x_c + d) - h(x_c) is
//   -E h1 - rho u_c E^2/(1+E) - a log1pmx(w_c E) - rho log1pmx(-u_c E/(1+E)),
// where h1 = (a - rho t_c)/(1 + t_c), evaluated by fma, is the only first-order term; and
// the level at the mode t* = a/rho is Stirling's
//   (1/2) log(rho a / ((rho + a) 2 pi)) - [corr(rho) + corr(a) - corr(rho + a)].
// For |d| > 1 the drop is -a log(u_c + w_c e^d) - rho log(w_c + u_c e^-d). The first sum
// is formed as 1 + w_c E, except below 2^-26, where it is added from its positive terms:
// for t_c < 2^-53, w_c rounds to 1 and w_c E to -1 once e^d < 2^-54, and the drop was +Inf
// (for a < ~1e-16 rho the search then missed the integrand's mode, ~log(rho/a) below x*,
// and P(Y > q) was returned as 1). The second is its mirror image, formed as
// 1 + u_c (e^-d - 1) except below 2^-26 (u_c rounds to 1 for t_c > 2^53; F(q) was returned as
// 1 for mu >~ 1e33 with k >~ 1e17).
// Small a also gives Beta(rho, a) a tail of mass ~e^(-a x) beyond x, which the quadrature
// cannot follow past |d| ~ 709, where e^d overflows. The upper-tail integrand vanishes there
// like e^(-(q+1) x); the lower-tail one does not, and lost e^(-709 a), which is below rounding
// only for a > 0.056 (log P(Y <= 0) was 0.49 too low at a = 1.3e-3). Below GW_A_LOWER = 1/16
// the lower tail is therefore log(1 - P(Y > q)).
// Continued fraction for the regularised incomplete beta I_x(a, b) (modified Lentz, as in
// Numerical Recipes' betacf), I_x(a, b) = x^a (1-x)^b / (a B(a, b)) * betacf; it converges
// quickly for x < (a + 1)/(a + b + 2).
static double cf_betacf(double a, double b, double x) {
  const double tiny = 1e-300;
  const double qab = a + b, qap = a + 1.0, qam = a - 1.0;
  double c = 1.0, d = 1.0 - qab * x / qap;
  if (std::fabs(d) < tiny) d = tiny;
  d = 1.0 / d;
  double h = d;
  for (int m = 1; m <= 100000; m++) {
    const double m2 = 2.0 * m;
    double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1.0 + aa * d; if (std::fabs(d) < tiny) d = tiny;
    c = 1.0 + aa / c; if (std::fabs(c) < tiny) c = tiny;
    d = 1.0 / d; h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1.0 + aa * d; if (std::fabs(d) < tiny) d = tiny;
    c = 1.0 + aa / c; if (std::fabs(c) < tiny) c = tiny;
    d = 1.0 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1.0) < 1e-16) break;
  }
  return h;
}
static inline double cf_log1mexp(double l) {                // log(1 - e^l), l <= 0
  return (l > -M_LN2) ? std::log(-std::expm1(l)) : std::log1p(-std::exp(l));
}
// log I_x(a, n + 1) = log P(NB(a, x) <= n), for an integer n >= 0 with (a + n - 1) y > n, as the
// finite sum of its negative-binomial terms, from the largest, T_n, down; T_{j-1}/T_j =
// j/((a + j - 1) y). In terms of y = 1 - x (lx = log x, accurate on both sides),
//   log T_n = log Pois(n; a y) + a [log(1 - y) + y] + [lgamma(a + n) - lgamma(a) - n log a],
// the last bracket being bd0(a, n) - log1p(n/a)/2 + corr(a + n) - corr(a), for a >= 10.
static double cf_lnbcdf(double a, double n, double y, double lx) {
  double l = R::dpois(n, a * y, 1) + a * ((y <= 0.5) ? R::log1pmx(-y) : lx + y);
  if (n > 0.0)
    l += cf_bd0(a, n) - 0.5 * std::log1p(n / a) + (cf_stirling_corr(a + n) - cf_stirling_corr(a));
  double s = 1.0, r = 1.0;
  for (double j = n; j > 0.0; j -= 1.0) {
    r *= j / ((a + j - 1.0) * y);
    s += r;
    if (r < 1e-20 * s) break;
  }
  return l + std::log(s);
}
// log P(Beta(a, b) <= x) (lower) or log P(Beta(a, b) > x), from x and y = 1 - x, both
// accurate. R's pbeta(log.p = TRUE) underflows to -Inf, with a warning, in far tails when a
// or b is large (bpser) -- P(Y > q) was then returned as 1 -- so a tail below e^-40 is
// evaluated in log space from the continued fraction, on whichever side it converges.
// The continued fraction for I_x(a, b) sees only x, so near x = 1 it loses the relative accuracy
// of y: its error is ~eps a / ((a + b) y - b + 1). For the negative-binomial size a = k this
// matters beyond CF_NB_HUGE: the region y ~ b/a puts x within ~b/a of 1, and the tails of the
// generalized Waring had relative errors of 1e-5 at k = 1e15 (P(NB <= q) 8e-4 off) and were NaN
// from k ~ 1e18, where x rounds to 1. There, with b = q + 1 an integer, the branch is chosen on
// y and the far lower tail is the finite sum above.
static const double CF_NB_HUGE = 1e12;
static double cf_lpbeta(double x, double y, double a, double b, bool lower) {
  if (!(x > 0.0)) return lower ? R_NegInf : 0.0;
  if (!(y > 0.0)) return lower ? 0.0 : R_NegInf;
  const double lx = (x < 0.5) ? std::log(x) : std::log1p(-y);
  const double ly = (y < 0.5) ? std::log(y) : std::log1p(-x);
  const double lpre = a * lx + b * ly - R::lbeta(a, b);
  const bool huge = a > CF_NB_HUGE && b >= 1.0 && b == std::floor(b) && b < a;
  if (huge ? y > (b + 1.0) / (a + b + 2.0) : x < (a + 1.0) / (a + b + 2.0)) {
    const double li = lpre - std::log(a);                  // log I_x(a, b), up to betacf
    if (li < -40.0) {
      const double l = huge ? cf_lnbcdf(a, b - 1.0, y, lx) : li + std::log(cf_betacf(a, b, x));
      return lower ? l : cf_log1mexp(l);
    }
  } else {
    const double lj = lpre - std::log(b);                  // log I_y(b, a), up to betacf
    if (lj < -40.0) {
      const double l = lj + std::log(cf_betacf(b, a, y));
      return lower ? cf_log1mexp(l) : l;
    }
  }
  return (x <= 0.5) ? R::pbeta(x, a, b, lower ? 1 : 0, 1) : R::pbeta(y, b, a, lower ? 0 : 1, 1);
}
struct GwMix { double a, rho, t, w, u, h1; };
static GwMix gw_mix(double a, double rho, double t) {
  GwMix m; m.a = a; m.rho = rho; m.t = t;
  m.w = 1.0 / (1.0 + t); m.u = t / (1.0 + t);
  m.h1 = std::fma(-rho, t, a) / (1.0 + t);
  return m;
}
static double gw_mix_c(const GwMix& m, double d) {
  if (d == 0.0) return 0.0;
  double E = std::expm1(d);
  if (std::fabs(d) > 1.0) {                  // the drop is then large: no cancellation,
    // except in 1 + w E = u + w e^d (d < 0) when t is tiny, and its mirror image
    // 1 + u (e^-d - 1) = w + u e^-d (d > 0) when t is huge: these are then summed from their terms
    const double wE = m.w * E, uF = m.u * std::expm1(-d);
    const double l1 = (1.0 + wE > 0x1p-26) ? std::log1p(wE) : std::log(m.u + m.w * std::exp(d));
    const double l2 = (1.0 + uF > 0x1p-26) ? std::log1p(uF) : std::log(m.w + m.u * std::exp(-d));
    return -m.a * l1 - m.rho * l2;
  }
  double F = -E / (1.0 + E);
  return -E * m.h1 - m.rho * m.u * (E * E / (1.0 + E)) - m.a * R::log1pmx(m.w * E) -
    m.rho * R::log1pmx(m.u * F);
}
static inline void gw_mix_wu(const GwMix& m, double d, double& w, double& u) {
  double v = m.t * std::exp(-d);
  if (v > 1.0) { double iv = 1.0 / v; u = 1.0 / (1.0 + iv); w = iv / (1.0 + iv); }
  else         { w = 1.0 / (1.0 + v);   u = v / (1.0 + v); }
}
static const double GW_A_LOWER = 0.0625;
static double gw_lcdf(double q, double a, double k, double rho, bool lower)
{
  if (lower && a < GW_A_LOWER) return cf_log1mexp(gw_lcdf(q, a, k, rho, false));
  // P(Beta(k, q+1) > w) (<= w for the lower tail) = P(NB(k, w) > q)
  auto lpb = [&](double w, double u) { return cf_lpbeta(w, u, k, q + 1.0, lower); };
  const double ts = a / rho;                          // mode of the mixing density
  const GwMix m0 = gw_mix(a, rho, ts);
  double w, u;
  gw_mix_wu(m0, 0.0, w, u);
  const double lp0 = lpb(w, u);
  auto g0 = [&](double d) {
    double wd, ud; gw_mix_wu(m0, d, wd, ud);
    return (lpb(wd, ud) - lp0) + gw_mix_c(m0, d);
  };
  const double sd0 = 1.0 / std::sqrt(std::max((a + rho) * m0.w * m0.u, 1e-300));
  const double ds = lc_argmax(g0, 0.0, std::min(1.0, 3.0 * sd0));
  const GwMix mc = gw_mix(a, rho, ts * std::exp(-ds));   // recentred at the integrand's mode
  gw_mix_wu(mc, 0.0, w, u);
  const double lpc = lpb(w, u);
  auto gc = [&](double d) {
    double wd, ud; gw_mix_wu(mc, d, wd, ud);
    return (lpb(wd, ud) - lpc) + gw_mix_c(mc, d);
  };
  double lmode;                                       // log mixing density at its mode
  if (std::min(a, rho) >= 10.0)
    lmode = 0.5 * std::log(rho * (a / (rho + a)) / (2.0 * M_PI)) -
      ((cf_stirling_corr(rho) + cf_stirling_corr(a)) - cf_stirling_corr(rho + a));
  else
    lmode = -rho * std::log1p(a / rho) - a * std::log1p(rho / a) - R::lbeta(rho, a);
  const double rt = (ts - mc.t) / mc.t;                    // x_c - x* = log(ts / t_c)
  const double lmix = lmode + gw_mix_c(m0, std::fabs(rt) <= 0.5 ? std::log1p(rt)
                                                                : std::log(ts / mc.t));
  double sc = lc_width(gc, 0.0, 1.0 / std::sqrt(std::max((a + rho) * mc.w * mc.u, 1e-300)));
  sc = std::min(std::max(sc, 1e-10), 10.0);
  double dl, dr;
  lc_limits(gc, 0.0, 0.0, sc, 50.0, dl, dr);
  std::vector<double> brk = lc_breaks(dl, 0.0, dr, sc);
  double I = gk_adapt(gc, brk, 0.0, 0.0, (double*)NULL, (double*)NULL);
  return std::min(lpc + lmix + std::log(I), 0.0);
}

// [[Rcpp::export]]
Rcpp::NumericVector gwaring_lcdf_cpp(Rcpp::NumericVector q, double mu, double k, double rho,
                                     bool lower) {
  int n = q.size();
  Rcpp::NumericVector out(n);
  double a = gwaring_a(mu, k, rho);
  for (int i = 0; i < n; i++) out[i] = gw_lcdf(q[i], a, k, rho, lower);
  return out;
}

// Cumulative log-sum-exp of a log-probability ladder: forward, out[i] = log sum_{j<=i} p_j;
// reverse, out[i] = log( sum_{j>=i} p_j + exp(ltail) ). Stable for terms far below double
// range, unlike cumsum(exp(.)).
// [[Rcpp::export]]
Rcpp::NumericVector count_lcse_cpp(Rcpp::NumericVector lp, bool reverse, double ltail) {
  int n = lp.size();
  Rcpp::NumericVector out(n);
  auto add = [](double x, double y) {
    if (x == -INFINITY) return y;
    if (y == -INFINITY) return x;
    double m = std::max(x, y);
    return m + std::log1p(std::exp(std::min(x, y) - m));
  };
  if (reverse) {
    double acc = ltail;
    for (int i = n - 1; i >= 0; i--) { acc = add(acc, lp[i]); out[i] = acc; }
  } else {
    double acc = -INFINITY;
    for (int i = 0; i < n; i++) { acc = add(acc, lp[i]); out[i] = acc; }
  }
  return out;
}

// The parameter clamps, exported so that predict.egpd() and the tests read the SAME numbers
// the likelihood enforces, rather than duplicating them as literals that can silently drift.
// [[Rcpp::export]]
Rcpp::List cf_bounds_cpp() {
  return Rcpp::List::create(
    Rcpp::Named("gpois")   = Rcpp::NumericVector::create(
      Rcpp::Named("mu_lo") = GP_MU_LO, Rcpp::Named("mu_hi") = GP_MU_HI,
      Rcpp::Named("lambda_lo") = GP_LAM_LO, Rcpp::Named("lambda_hi") = GP_LAM_HI),
    Rcpp::Named("gwaring") = Rcpp::NumericVector::create(
      Rcpp::Named("mu_lo") = GW_MU_LO, Rcpp::Named("mu_hi") = GW_MU_HI,
      Rcpp::Named("k_lo") = GW_K_LO, Rcpp::Named("k_hi") = GW_K_HI,
      Rcpp::Named("rho_lo") = GW_RHO_LO, Rcpp::Named("rho_hi") = GW_RHO_HI),
    Rcpp::Named("plnorm")  = Rcpp::NumericVector::create(
      Rcpp::Named("mu_lo") = PLN_MU_LO, Rcpp::Named("mu_hi") = PLN_MU_HI,
      Rcpp::Named("sigma_lo") = PLN_SG_LO, Rcpp::Named("sigma_hi") = PLN_SG_HI));
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
