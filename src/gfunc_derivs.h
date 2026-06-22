#ifndef EGPD_GFUNC_DERIVS_H
#define EGPD_GFUNC_DERIVS_H
// Shared closed-form derivative building blocks for the truncated-normal (model 5)
// and truncated-beta (model 6) G-functions, used to give analytic / hybrid gradients
// and Hessians for the discrete families (DEGPD, ZIDEGPD). The same math is verified
// inline in degpd.cpp (degpd5d12 analytic, degpd6d12 hybrid) against the FD oracle.
#include <Rcpp.h>
#include <cmath>

// --- discrete GPD CDF F = 1 - (1+xi m/sigma)^(-1/xi) (log-xi link) -------------
// F and its first/second derivatives w.r.t. (lsigma, lxi). ok=false if support is violated.
struct FDeriv { double F, F_ls, F_lx, F_lsls, F_lslx, F_lxlx; bool ok; };

static inline FDeriv gpd_Fderiv(double m, double sigma, double xi) {
  FDeriv r;
  double q = 1.0 / xi, v = xi * m / sigma, t = 1.0 + v;
  if (t <= 0.0) { r.ok = false; r.F = r.F_ls = r.F_lx = r.F_lsls = r.F_lslx = r.F_lxlx = 0.0; return r; }
  double L = std::log(t), w = std::pow(t, -q);
  r.F = 1.0 - w;
  double e_ls = q * v / t;
  double e_lx = q * (L - v / t);
  double e_lsls = -q * (v / t - v * v / (t * t));
  double e_lslx = -q * v * v / (t * t);
  double e_lxlx = -q * (L - v / t - v * v / (t * t));
  r.F_ls = -w * e_ls;
  r.F_lx = -w * e_lx;
  r.F_lsls = -w * (e_ls * e_ls + e_lsls);
  r.F_lslx = -w * (e_ls * e_lx + e_lslx);
  r.F_lxlx = -w * (e_lx * e_lx + e_lxlx);
  r.ok = true;
  return r;
}

// G and its theta-derivatives (theta = lsigma, lxi, lkappa) at one F point.
struct GTheta { double G, ls, lx, lk, lsls, lslx, lslk, lxlx, lxlk, lklk; };

// --- truncated-normal G: kappa-side precompute (shared across both F endpoints) --
struct TNkappa { double kappa, sk, skp, skpp, Fmin, phisk, D, Dk, Dkk, Fmin_kk; };

static inline TNkappa tn_kappa(double kappa) {
  TNkappa k;
  k.kappa = kappa;
  k.sk = std::sqrt(kappa);
  k.skp = 0.5 / k.sk;
  k.skpp = -k.skp / (2.0 * kappa);
  k.Fmin = R::pnorm(-k.sk, 0.0, 1.0, 1, 0);
  k.phisk = R::dnorm(k.sk, 0.0, 1.0, 0);
  k.D = 0.5 - k.Fmin; if (k.D < 1e-300) k.D = 1e-300;
  k.Dk = k.skp * k.phisk;
  k.Dkk = k.phisk * (k.skpp - k.skp * k.skp * k.sk);
  k.Fmin_kk = k.phisk * (k.skp * k.skp * k.sk - k.skpp);
  return k;
}

// truncated-normal G and full theta-derivatives at an F point (fully closed form).
static inline GTheta tn_Gtheta(const FDeriv& f, const TNkappa& k) {
  GTheta g;
  double sk = k.sk, D = k.D, kappa = k.kappa;
  double z = sk * (f.F - 1.0);
  double phiz = R::dnorm(z, 0.0, 1.0, 0);
  double Phiz = R::pnorm(z, 0.0, 1.0, 1, 0);
  double N = Phiz - k.Fmin;
  double G = N / D;
  double G_F = phiz * sk / D;
  double G_FF = -sk * sk * z * phiz / D;
  double zk = k.skp * (f.F - 1.0);
  double zkk = k.skpp * (f.F - 1.0);
  double Nk = phiz * zk + k.skp * k.phisk;
  double G_K = (Nk * D - N * k.Dk) / (D * D);
  double Nkk = phiz * (zkk - z * zk * zk) - k.Fmin_kk;
  double G_KK = (Nkk * D - N * k.Dkk) / (D * D) - 2.0 * k.Dk * G_K / D;
  double G_FK = ((-z * phiz) * zk * sk + phiz * k.skp) / D - phiz * sk * k.Dk / (D * D);
  g.G    = G;
  g.ls   = G_F * f.F_ls;
  g.lx   = G_F * f.F_lx;
  g.lk   = G_K * kappa;
  g.lsls = G_FF * f.F_ls * f.F_ls + G_F * f.F_lsls;
  g.lslx = G_FF * f.F_ls * f.F_lx + G_F * f.F_lslx;
  g.lxlx = G_FF * f.F_lx * f.F_lx + G_F * f.F_lxlx;
  g.lslk = G_FK * kappa * f.F_ls;
  g.lxlk = G_FK * kappa * f.F_lx;
  g.lklk = G_KK * kappa * kappa + G_K * kappa;
  return g;
}

// --- truncated-beta G: sigma/xi block only (kappa handled numerically) ----------
struct GThetaSX { double G, ls, lx, lsls, lslx, lxlx; };

static inline GThetaSX tb_GthetaSX(const FDeriv& f, double kappa, double pb_min, double D) {
  const double c1 = 0.5 - 1.0 / 32.0, c2 = 1.0 / 32.0;
  double u = c1 * f.F + c2;
  double db = R::dbeta(u, kappa, kappa, 0);
  double G_F = c1 * db / D;
  double dbu = (kappa - 1.0) * db * (1.0 - 2.0 * u) / (u * (1.0 - u));
  double G_FF = c1 * c1 * dbu / D;
  GThetaSX g;
  g.G    = (R::pbeta(u, kappa, kappa, 1, 0) - pb_min) / D;
  g.ls   = G_F * f.F_ls;
  g.lx   = G_F * f.F_lx;
  g.lsls = G_FF * f.F_ls * f.F_ls + G_F * f.F_lsls;
  g.lslx = G_FF * f.F_ls * f.F_lx + G_F * f.F_lslx;
  g.lxlx = G_FF * f.F_lx * f.F_lx + G_F * f.F_lxlx;
  return g;
}

// --- truncated-beta G: FULL theta-derivatives (analytic kappa via quadrature) -----
// d/dkappa I(x;k,k) has no elementary form; use
//   dI/dk = K1 - I*b1,  d2I/dk2 = K2 - 2*K1*b1 + I*(b1^2 - b1'),
// Km = int_0^x dbeta(t;k,k) [ln(t(1-t))]^m dt, b1 = 2(psi(k)-psi(2k)),
// b1' = 2(psi'(k)-2 psi'(2k)). Km via substitution t = x z^(1/k) (removes the
// t^(k-1) endpoint) + adaptive Gauss-Kronrod 7-15. Validated to ~1e-9 vs pbeta FD.
static const double EGPD_XGK[8] = {0.991455371120813,0.949107912342759,0.864864423359769,
  0.741531185599394,0.586087235467691,0.405845151377397,0.207784955007898,0.0};
static const double EGPD_WGK[8] = {0.022935322010529,0.063092092629979,0.104790010322250,
  0.140653259715525,0.169004726639267,0.190350578064785,0.204432940075298,0.209482141084728};
static const double EGPD_WG[4] = {0.129484966168870,0.279705391489277,0.381830050505119,0.417959183673469};

// substituted integrand at z in (0,1): g1 = w*L, g2 = w*L^2  (L = ln(t(1-t)), t = x z^(1/k))
static inline void egpd_tb_intg(double z, double x, double kappa, double c, double lxkB,
                                double& g1, double& g2) {
  if (z <= 0.0) { g1 = g2 = 0.0; return; }
  double t = x * std::pow(z, c);
  double w = std::exp(lxkB + (kappa - 1.0) * std::log(1.0 - t));
  double L = std::log(x) + std::log(1.0 - t) + c * std::log(z);
  g1 = w * L; g2 = w * L * L;
}

static inline void egpd_tb_qk15(double a, double b, double x, double kappa, double c, double lxkB,
                                double& v1, double& v2, double& e1, double& e2) {
  double ctr = 0.5 * (a + b), hl = 0.5 * (b - a);
  double f1c, f2c; egpd_tb_intg(ctr, x, kappa, c, lxkB, f1c, f2c);
  double rk1 = EGPD_WGK[7] * f1c, rk2 = EGPD_WGK[7] * f2c;
  double rg1 = EGPD_WG[3] * f1c,  rg2 = EGPD_WG[3] * f2c;
  for (int j = 0; j < 7; j++) {
    double dx = hl * EGPD_XGK[j];
    double l1, l2, r1, r2;
    egpd_tb_intg(ctr - dx, x, kappa, c, lxkB, l1, l2);
    egpd_tb_intg(ctr + dx, x, kappa, c, lxkB, r1, r2);
    double s1 = l1 + r1, s2 = l2 + r2;
    rk1 += EGPD_WGK[j] * s1; rk2 += EGPD_WGK[j] * s2;
    if (j % 2 == 1) { rg1 += EGPD_WG[(j - 1) / 2] * s1; rg2 += EGPD_WG[(j - 1) / 2] * s2; }
  }
  v1 = rk1 * hl; v2 = rk2 * hl;
  e1 = std::fabs((rk1 - rg1) * hl); e2 = std::fabs((rk2 - rg2) * hl);
}

static inline void egpd_tb_adapt(double a, double b, double x, double kappa, double c, double lxkB,
                                 double tol, int depth, double& K1, double& K2) {
  double v1, v2, e1, e2; egpd_tb_qk15(a, b, x, kappa, c, lxkB, v1, v2, e1, e2);
  double s1 = std::fabs(v1) > 1e-3 ? std::fabs(v1) : 1e-3;
  double s2 = std::fabs(v2) > 1e-3 ? std::fabs(v2) : 1e-3;
  if ((e1 < tol * s1 && e2 < tol * s2) || (b - a) < 1e-9 || depth > 30) { K1 = v1; K2 = v2; return; }
  double m = 0.5 * (a + b), aK1, aK2, bK1, bK2;
  egpd_tb_adapt(a, m, x, kappa, c, lxkB, tol, depth + 1, aK1, aK2);
  egpd_tb_adapt(m, b, x, kappa, c, lxkB, tol, depth + 1, bK1, bK2);
  K1 = aK1 + bK1; K2 = aK2 + bK2;
}

// dI/dkappa and d2I/dkappa2 of the regularized incomplete beta I(x;kappa,kappa).
static inline void incbeta_kk_derivs(double x, double kappa, double& dI, double& d2I) {
  double c = 1.0 / kappa;
  double lxkB = kappa * std::log(x) - std::log(kappa) - R::lbeta(kappa, kappa);
  double K1, K2; egpd_tb_adapt(0.0, 1.0, x, kappa, c, lxkB, 1e-11, 0, K1, K2);
  double I = R::pbeta(x, kappa, kappa, 1, 0);
  double b1 = 2.0 * (R::digamma(kappa) - R::digamma(2.0 * kappa));
  double b1p = 2.0 * (R::trigamma(kappa) - 2.0 * R::trigamma(2.0 * kappa));
  dI = K1 - I * b1;
  d2I = K2 - 2.0 * K1 * b1 + I * (b1 * b1 - b1p);
}

// kappa-side precompute for truncated-beta G (shared across both F endpoints).
struct TBkappa { double kappa, c1, c2, I_c2, dI_c2, d2I_c2, D, Dk, Dkk, b1; };

static inline TBkappa tb_kappa(double kappa) {
  TBkappa k; k.kappa = kappa; k.c1 = 0.5 - 1.0 / 32.0; k.c2 = 1.0 / 32.0;
  double dI_h, d2I_h;
  incbeta_kk_derivs(k.c2, kappa, k.dI_c2, k.d2I_c2);
  incbeta_kk_derivs(0.5,  kappa, dI_h,    d2I_h);
  k.I_c2 = R::pbeta(k.c2, kappa, kappa, 1, 0);
  double I_h = R::pbeta(0.5, kappa, kappa, 1, 0);
  k.D = I_h - k.I_c2; if (k.D < 1e-300) k.D = 1e-300;
  k.Dk = dI_h - k.dI_c2;
  k.Dkk = d2I_h - k.d2I_c2;
  k.b1 = 2.0 * (R::digamma(kappa) - R::digamma(2.0 * kappa));
  return k;
}

// truncated-beta G and FULL theta-derivatives at an F point (kappa via quadrature).
static inline GTheta tb_Gtheta(const FDeriv& f, const TBkappa& k) {
  double kappa = k.kappa, D = k.D, c1 = k.c1, c2 = k.c2;
  double u = c1 * f.F + c2;
  double dI_u, d2I_u; incbeta_kk_derivs(u, kappa, dI_u, d2I_u);
  double I_u = R::pbeta(u, kappa, kappa, 1, 0);
  double N = I_u - k.I_c2, Nk = dI_u - k.dI_c2, Nkk = d2I_u - k.d2I_c2;
  double G = N / D;
  double G_K = (Nk * D - N * k.Dk) / (D * D);
  double G_KK = Nkk / D - 2.0 * Nk * k.Dk / (D * D) - N * k.Dkk / (D * D) + 2.0 * N * k.Dk * k.Dk / (D * D * D);
  double db = R::dbeta(u, kappa, kappa, 0);
  double G_F = c1 * db / D;
  double dbu = (kappa - 1.0) * db * (1.0 - 2.0 * u) / (u * (1.0 - u));
  double G_FF = c1 * c1 * dbu / D;
  double ddb = db * (std::log(u * (1.0 - u)) - k.b1);   // d(beta density)/dkappa
  double G_FK = c1 * (ddb * D - db * k.Dk) / (D * D);
  GTheta g;
  g.G    = G;
  g.ls   = G_F * f.F_ls;
  g.lx   = G_F * f.F_lx;
  g.lk   = G_K * kappa;
  g.lsls = G_FF * f.F_ls * f.F_ls + G_F * f.F_lsls;
  g.lslx = G_FF * f.F_ls * f.F_lx + G_F * f.F_lslx;
  g.lxlx = G_FF * f.F_lx * f.F_lx + G_F * f.F_lxlx;
  g.lslk = G_FK * kappa * f.F_ls;
  g.lxlk = G_FK * kappa * f.F_lx;
  g.lklk = G_KK * kappa * kappa + G_K * kappa;
  return g;
}

#endif
