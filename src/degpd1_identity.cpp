// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
//
// IDENTITY-LINK, bounded-support variant of DEGPD model 1.
//
// The standard degpd1d0/degpd1d12 parameterise the shape as lxi = log(xi), so
// xi = exp(.) > 0 and a bounded/short tail (xi < 0) cannot be represented.
// Here the shape uses an IDENTITY link: parameter 2 IS xi (real, may be < 0),
// giving a finite upper endpoint -sigma/xi. Observations beyond the endpoint
// have zero probability; we return a large finite penalty so the optimiser is
// pushed back inside the support.
//
// Parameters: (lsigma, xi, lkappa). sigma = exp(lsigma), kappa = exp(lkappa).
//   H(t) = 1 - (1 + xi t/sigma)^(-1/xi);  F(t) = H(t)^kappa;  P(X=y)=F(y+1)-F(y)
// The gradient and Hessian (in degpd1id12) were derived symbolically with SymPy
// (derive_degpd1_identity.py) and verified against finite differences.

// NOTE (xi ~ 0): the closed-form gradient/Hessian below match SymPy to ~1e-8 for
// |xi| >~ 0.05, but for |xi| < ~0.02 double-precision catastrophic cancellation
// (terms ~1/xi^2 that nearly cancel) degrades accuracy. That band is the light-tail
// boundary (the model -> geometric as xi -> 0), where the bounded-vs-light question
// is moot, so this is benign; guard_xi only prevents the exact 1/xi singularity.
static const double XI_EPS = 1e-6;          // avoid the 1/xi singularity at xi=0
static const double OOS_PEN = 1e20;         // out-of-support penalty for d0

static inline double guard_xi(double xi) {
  if (std::fabs(xi) < XI_EPS) return (xi < 0 ? -XI_EPS : XI_EPS);
  return xi;
}

// negative log-likelihood (scalar)
// [[Rcpp::export]]
double degpd1id0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2,
                 const arma::mat& X3, arma::vec yvec, const arma::uvec& dupid,
                 int dcate, const Rcpp::List& offsets, double xi_max)
{
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec xivec     = X2 * Rcpp::as<arma::vec>(pars[1]);   // identity link
  arma::vec lkappavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = yvec.size();
  if (dcate == 1) { lsigmavec = lsigmavec.elem(dupid); xivec = xivec.elem(dupid); lkappavec = lkappavec.elem(dupid); }
  { arma::vec o0 = Rcpp::as<arma::vec>(offsets[0]); if (o0.n_elem > 0) lsigmavec += o0;
    arma::vec o1 = Rcpp::as<arma::vec>(offsets[1]); if (o1.n_elem > 0) xivec     += o1;
    arma::vec o2 = Rcpp::as<arma::vec>(offsets[2]); if (o2.n_elem > 0) lkappavec += o2; }

  double nllh = 0.0;
  for (int j = 0; j < nobs; j++) {
    double y = yvec[j];
    double sigma = exp(lsigmavec[j]);
    double xi = guard_xi(xivec[j]);
    double kappa = exp(lkappavec[j]);
    double e2 = xi / sigma;
    double e4 = 1.0 + y * e2;            // 1 + xi*y/sigma
    double e3 = 1.0 + (y + 1.0) * e2;    // 1 + xi*(y+1)/sigma
    if (e4 <= 0.0) { nllh += OOS_PEN; continue; }                       // y beyond endpoint -> P=0
    double hi = (e3 <= 0.0) ? 1.0 : R_pow(1.0 - R_pow(1.0 / e3, 1.0 / xi), kappa); // F(y+1); =1 at endpoint cell
    double lo = R_pow(1.0 - R_pow(1.0 / e4, 1.0 / xi), kappa);          // F(y)
    double p = hi - lo;
    if (p <= 0.0 || !R_finite(p)) { nllh += OOS_PEN; continue; }
    nllh += -log(p);
  }
  return nllh;
}

// gradient (cols 0-2) + upper-tri Hessian (cols 3-8) w.r.t. (lsigma, xi, lkappa)
// [[Rcpp::export]]
arma::mat degpd1id12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3,
                     arma::vec yvec, const arma::uvec dupid, int dcate,
                     const Rcpp::List& offsets, double xi_max)
{
  arma::vec lsigmavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec xivec     = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec lkappavec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = yvec.size();
  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  if (dcate == 1) { lsigmavec = lsigmavec.elem(dupid); xivec = xivec.elem(dupid); lkappavec = lkappavec.elem(dupid); }
  { arma::vec o0 = Rcpp::as<arma::vec>(offsets[0]); if (o0.n_elem > 0) lsigmavec += o0;
    arma::vec o1 = Rcpp::as<arma::vec>(offsets[1]); if (o1.n_elem > 0) xivec     += o1;
    arma::vec o2 = Rcpp::as<arma::vec>(offsets[2]); if (o2.n_elem > 0) lkappavec += o2; }

  for (int j = 0; j < nobs; j++) {
    double y = yvec[j];
    double lsigma = lsigmavec[j];
    double xi = guard_xi(xivec[j]);
    double lkappa = lkappavec[j];
    // support check (e4, e3 = x2, x14 below)
    double e2 = xi / exp(lsigma);
    double e4g = 1.0 + y * e2, e3g = 1.0 + (y + 1.0) * e2;
    if (e4g <= 0.0) continue;                         // y beyond endpoint -> row stays 0
    double OUT0, OUT1, OUT2, OUT3, OUT4, OUT5, OUT6, OUT7, OUT8;
    if (y < 0.5) {
    // ---- AUTO-GENERATED y=0 cell: P=F(1) (F(0)=0 identically) ----
    double x0 = exp(lkappa);
    double x1 = exp(-lsigma);
    double x2 = x1*xi;
    double x3 = x2 + 1;
    double x4 = 1.0/x3;
    double x5 = x1*x4;
    double x6 = 1.0/xi;
    double x7 = R_pow(x3, -x6);
    double x8 = 1 - x7;
    double x9 = x7/x8;
    double x10 = x5*x9;
    double x11 = x0*x10;
    double x12 = log(x3);
    double x13 = -x12*x6 + x5;
    double x14 = x13*x6;
    double x15 = x14*x9;
    double x16 = -x0*x15;
    double x17 = -x0*log(x8);
    double x18 = R_pow(x13, 2)*x6;
    OUT0 = x11;
    OUT1 = x16;
    OUT2 = x17;
    OUT3 = x11*(x10 + x2*x4 + x5 - 1);
    OUT4 = -x11*(x14 + x15 + x5);
    OUT5 = x11;
    OUT6 = x0*x6*x9*(-2*x12/R_pow(xi, 2) + x18*x9 + x18 + 2*x5*x6 + exp(-2*lsigma)/R_pow(x3, 2));
    OUT7 = x16;
    OUT8 = x17;
    // ---- end generated ----
    } else if (e3g > 0.0) {
    // ---- AUTO-GENERATED (derive_degpd1_identity.py); interior cell P=F(y+1)-F(y) ----
    double x0 = exp(-lsigma);
    double x1 = x0*xi;
    double x2 = x1*y + 1;
    double x3 = y/x2;
    double x4 = 1.0/xi;
    double x5 = R_pow(x2, -x4);
    double x6 = 1 - x5;
    double x7 = exp(lkappa);
    double x8 = R_pow(x6, x7);
    double x9 = 1.0/x6;
    double x10 = x5*x9;
    double x11 = x10*x8;
    double x12 = x11*x3;
    double x13 = y + 1;
    double x14 = x1*x13 + 1;
    double x15 = x13/x14;
    double x16 = R_pow(x14, -x4);
    double x17 = 1 - x16;
    double x18 = R_pow(x17, x7);
    double x19 = x16/x17;
    double x20 = x18*x19;
    double x21 = x15*x20;
    double x22 = x12 - x21;
    double x23 = x7/(-x18 + x8);
    double x24 = x0*x23;
    double x25 = log(x2);
    double x26 = x0*x3;
    double x27 = -x25*x4 + x26;
    double x28 = log(x14);
    double x29 = x0*x15;
    double x30 = -x28*x4 + x29;
    double x31 = x18*x30;
    double x32 = x19*x31;
    double x33 = -x27*x5*x8*x9 + x32;
    double x34 = x23*x4;
    double x35 = log(x6);
    double x36 = x35*x8;
    double x37 = log(x17);
    double x38 = x18*x37;
    double x39 = x36 - x38;
    double x40 = x23*x39;
    double x41 = R_pow(y, 2)/R_pow(x2, 2);
    double x42 = x0*x41;
    double x43 = x11*x42;
    double x44 = 2*x4;
    double x45 = R_pow(x2, -x44)/R_pow(x6, 2);
    double x46 = x45*x8;
    double x47 = R_pow(x13, 2)/R_pow(x14, 2);
    double x48 = x0*x47;
    double x49 = x20*x48;
    double x50 = R_pow(x14, -x44)/R_pow(x17, 2);
    double x51 = x18*x50;
    double x52 = x18*x7;
    double x53 = x50*x52;
    double x54 = x7*x8;
    double x55 = x45*x54;
    double x56 = -x12 + x21;
    double x57 = x27*x8;
    double x58 = x3*x4;
    double x59 = x10*x57;
    double x60 = x15*x4;
    double x61 = x19*x38*x7;
    double x62 = x10*x36*x7;
    double x63 = R_pow(x27, 2)*x4;
    double x64 = R_pow(x30, 2)*x4;
    double x65 = 2/R_pow(xi, 2);
    double x66 = exp(-2*lsigma);
    OUT0 = x22*x24;
    OUT1 = x33*x34;
    OUT2 = -x40;
    OUT3 = x24*(x1*x11*x41 - x1*x20*x47 + R_pow(x22, 2)*x24 + x42*x46 - x42*x55 + x43 - x48*x51 + x48*x53 - x49 + x56);
    OUT4 = x24*(x22*x33*x34 + x27*x55*x58 - x30*x53*x60 + x31*x50*x60 + x32*x60 - x43 - x45*x57*x58 + x49 - x58*x59);
    OUT5 = -x24*(x15*x61 + x22*x40 - x3*x62 + x56);
    OUT6 = x34*(x11*x63 + x11*(-x25*x65 + x26*x44 + x41*x66) - x20*x64 - x20*(-x28*x65 + x29*x44 + x47*x66) + R_pow(x33, 2)*x34 + x46*x63 - x51*x64 + x53*x64 - x55*x63);
    OUT7 = -x34*(x27*x62 - x30*x61 - x32 + x33*x40 + x59);
    OUT8 = x23*(x23*R_pow(x39, 2) - R_pow(x35, 2)*x54 - x36 + R_pow(x37, 2)*x52 + x38);
    // ---- end generated ----
    } else {
    // ---- AUTO-GENERATED endpoint cell: P = 1 - F(y) (derive_degpd1_identity.py) ----
    double x0 = y*exp(-lsigma);
    double x1 = x0*xi;
    double x2 = x1 + 1;
    double x3 = 1.0/xi;
    double x4 = R_pow(x2, -x3);
    double x5 = 1 - x4;
    double x6 = exp(lkappa);
    double x7 = R_pow(x5, x6);
    double x8 = x7/(x7 - 1);
    double x9 = 1.0/x2;
    double x10 = x0*x9;
    double x11 = x4/x5;
    double x12 = x10*x11;
    double x13 = x12*x6;
    double x14 = x13*x8;
    double x15 = log(x2);
    double x16 = x10 - x15*x3;
    double x17 = x16*x3;
    double x18 = x11*x17;
    double x19 = x18*x6;
    double x20 = x19*x8;
    double x21 = x6*log(x5);
    double x22 = x21*x8;
    double x23 = -x22;
    double x24 = x21 + x23 + 1;
    double x25 = R_pow(x16, 2)*x3;
    double x26 = x11*x25;
    double x27 = x26*x6;
    OUT0 = x14;
    OUT1 = -x20;
    OUT2 = x23;
    OUT3 = x14*(x1*x9 + x10 + x12 - x13 + x14 - 1);
    OUT4 = -x14*(x10 + x17 + x18 - x19 + x20);
    OUT5 = x14*x24;
    OUT6 = x11*x3*x6*x8*(2*x10*x3 - 2*x15/R_pow(xi, 2) + x25 + x26 + x27*x8 - x27 + R_pow(y, 2)*exp(-2*lsigma)/R_pow(x2, 2));
    OUT7 = -x20*x24;
    OUT8 = -x22*x24;
    // ---- end generated ----
    }

    out(j, 0) = OUT0; out(j, 1) = OUT1; out(j, 2) = OUT2;
    out(j, 3) = OUT3; out(j, 4) = OUT4; out(j, 5) = OUT5;
    out(j, 6) = OUT6; out(j, 7) = OUT7; out(j, 8) = OUT8;
  }
  return out;
}
