#include <Rcpp.h>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <limits>

// Only used when ordinary probability arithmetic cannot separate a
// power-carrier zero mass from the empirical fraction. Work with the
// supplied binary parameters, not their printed decimal approximations.
// [[Rcpp::export(rng = false)]]
double profile_zero_odds_cpp(double zeros, double n, double sigma, double xi,
                             double kappa, double delta, double weight) {
  using Real = boost::multiprecision::cpp_bin_float_100;
  const Real s(sigma), x(xi);
  if (x <= -s) return R_NegInf;
  // For xi = +/-1 the GP mass is 1/d. Clearing the denominator keeps
  // rational equalities exact, e.g. one zero in 7^3 for d=7, kappa=3.
  if ((xi == 1 || xi == -1) && (weight == 1 || kappa == delta)) {
    const Real d = xi == 1 ? s+1 : s;
    const Real dk = pow(d, Real(kappa));
    const Real gap = zeros*dk-n;
    return gap > 0 ? static_cast<double>(log(gap/dk/(n-zeros))) : R_NegInf;
  }
  const Real h = xi == 0 ? 1/s : log1p(x/s)/x;
  const Real lu = h > 1 ? log1p(-exp(-h)) : log(-expm1(-h));
  const Real b = Real(weight)*exp(Real(kappa)*lu) + (1-Real(weight))*exp(Real(delta)*lu);
  const Real gap = zeros - n*b;
  // Rounding at working precision in the transcendental calculation.
  return gap > 32*std::numeric_limits<Real>::epsilon()*zeros ?
    static_cast<double>(log(gap/(n-zeros))) : R_NegInf;
}

// Two jump cells: the Bernoulli or thinned Poisson mean is affine in
// the mixture weight. Compute the two distances from the empirical mean
// separately, including when a positive interior weight is close to zero.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector two_cell_weight_cpp(double total, double n, double sigma,
                                        double xi, double kappa, double delta,
                                        double lambda, double pi, bool compound) {
  using Real = boost::multiprecision::cpp_bin_float_100;
  const Real x(xi), s(sigma);
  const Real h = log1p(x/s)/x;
  const Real lu = h > 1 ? log1p(-exp(-h)) : log(-expm1(-h));
  const Real scale = compound ? Real(lambda) : 1-Real(pi);
  const Real ns = n*scale;
  auto gap = [&](double power) -> Real {
    const Real l = Real(power)*lu;
    const bool small_zero = l < -log(Real(2));
    const Real a = small_zero ? Real(total)-ns : Real(total);
    const Real b = ns*(small_zero ? exp(l) : expm1(l));
    const Real g = a+b;
    return abs(g) <= 32*std::numeric_limits<Real>::epsilon()*(abs(a)+abs(b)) ? Real(0) : g;
  };
  const Real ga = gap(kappa), gb = gap(delta);
  const bool increasing = kappa > delta;
  if (increasing ? gb <= 0 : gb >= 0)
    return Rcpp::NumericVector::create(R_NegInf, NA_REAL);
  if (increasing ? ga >= 0 : ga <= 0)
    return Rcpp::NumericVector::create(R_PosInf, NA_REAL);
  const Real mu = Real(total)/n;
  const Real variance = compound ? mu : mu*(1-mu);
  return Rcpp::NumericVector::create(static_cast<double>(log(abs(gb))-log(abs(ga))),
                                     static_cast<double>(sqrt(n*variance)/abs(gb-ga)));
}
