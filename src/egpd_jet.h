#ifndef EGPD_JET_H
#define EGPD_JET_H
// Second-order forward-mode automatic differentiation ("Jets") shared by the
// stabilised discrete (DEGPD/ZIDEGPD) evaluators and the continuous-EGPD
// derivative fallback. Jet<N> carries a value, gradient and Hessian in N
// active directions. Jet<0> is value-only, so objective-only (d0) calls do not
// pay for derivative propagation; for N > 0 every operation performs exactly
// the arithmetic of the former fixed-size Jet, so derivatives are unchanged.
#include <cmath>
#include <Rcpp.h>

namespace egpd_ad {

template<int N>
struct Jet {
  static constexpr int M = N > 0 ? N : 1;
  double v, g[M], h[M][M];
  Jet(double value = 0.0) : v(value) {
    for (int i = 0; i < M; ++i) {
      g[i] = 0.0;
      for (int j = 0; j < M; ++j) h[i][j] = 0.0;
    }
  }
};

template<int N> inline Jet<N> variable(double value, int i) {
  Jet<N> x(value);
  if (i < N) x.g[i] = 1.0;
  return x;
}
template<int N> inline Jet<N> operator+(const Jet<N>& a, const Jet<N>& b) {
  Jet<N> c(a.v + b.v);
  for (int i = 0; i < N; ++i) {
    c.g[i] = a.g[i] + b.g[i];
    for (int j = 0; j < N; ++j) c.h[i][j] = a.h[i][j] + b.h[i][j];
  }
  return c;
}
template<int N> inline Jet<N> operator-(const Jet<N>& a) {
  Jet<N> c(-a.v);
  for (int i = 0; i < N; ++i) {
    c.g[i] = -a.g[i];
    for (int j = 0; j < N; ++j) c.h[i][j] = -a.h[i][j];
  }
  return c;
}
template<int N> inline Jet<N> operator-(const Jet<N>& a, const Jet<N>& b) { return a + (-b); }
template<int N> inline Jet<N> operator*(const Jet<N>& a, const Jet<N>& b) {
  Jet<N> c(a.v * b.v);
  for (int i = 0; i < N; ++i) {
    c.g[i] = a.g[i] * b.v + a.v * b.g[i];
    for (int j = 0; j < N; ++j)
      c.h[i][j] = a.h[i][j] * b.v + a.g[i] * b.g[j] +
        a.g[j] * b.g[i] + a.v * b.h[i][j];
  }
  return c;
}
template<int N> inline Jet<N> operator*(double s, const Jet<N>& a) {
  Jet<N> c(s * a.v);
  for (int i = 0; i < N; ++i) {
    c.g[i] = s * a.g[i];
    for (int j = 0; j < N; ++j) c.h[i][j] = s * a.h[i][j];
  }
  return c;
}
template<int N> inline Jet<N> operator*(const Jet<N>& a, double s) { return s * a; }
// Mixed Jet/double sums behave exactly like promoting the double to a constant Jet.
template<int N> inline Jet<N> operator+(const Jet<N>& a, double b) { return a + Jet<N>(b); }
template<int N> inline Jet<N> operator+(double a, const Jet<N>& b) { return Jet<N>(a) + b; }
template<int N> inline Jet<N> operator-(const Jet<N>& a, double b) { return a - Jet<N>(b); }
template<int N> inline Jet<N> operator-(double a, const Jet<N>& b) { return Jet<N>(a) - b; }
template<int N> inline Jet<N> unary(const Jet<N>& x, double value, double d1, double d2) {
  Jet<N> c(value);
  for (int i = 0; i < N; ++i) {
    c.g[i] = d1 * x.g[i];
    for (int j = 0; j < N; ++j)
      c.h[i][j] = d1 * x.h[i][j] + d2 * x.g[i] * x.g[j];
  }
  return c;
}
template<int N> inline Jet<N> reciprocal(const Jet<N>& x) {
  double r = 1.0 / x.v;
  Jet<N> c(r);
  for (int i = 0; i < N; ++i) c.g[i] = -r * (x.g[i] * r);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      c.h[i][j] = -r * (x.h[i][j] * r) + 2.0 * x.v * c.g[i] * c.g[j];
  return c;
}
template<int N> inline Jet<N> operator/(const Jet<N>& a, const Jet<N>& b) { return a * reciprocal(b); }
template<int N> inline Jet<N> operator/(const Jet<N>& a, double b) { return a * reciprocal(Jet<N>(b)); }
template<int N> inline Jet<N> operator/(double a, const Jet<N>& b) { return Jet<N>(a) * reciprocal(b); }
template<int N> inline Jet<N> jexp(const Jet<N>& x) { double e = std::exp(x.v); return unary(x, e, e, e); }
template<int N> inline Jet<N> jlog(const Jet<N>& x) {
  Jet<N> c(std::log(x.v));
  for (int i = 0; i < N; ++i) c.g[i] = x.g[i] / x.v;
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      c.h[i][j] = x.h[i][j] / x.v - c.g[i] * c.g[j];
  return c;
}
template<int N> inline Jet<N> jlog1p(const Jet<N>& x) {
  double r = 1.0 / (1.0 + x.v);
  Jet<N> c(std::log1p(x.v));
  for (int i = 0; i < N; ++i) c.g[i] = x.g[i] * r;
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      c.h[i][j] = x.h[i][j] * r - c.g[i] * c.g[j];
  return c;
}
template<int N> inline Jet<N> jexpm1(const Jet<N>& x) {
  double e = std::exp(x.v);
  return unary(x, std::expm1(x.v), e, e);
}
// (exp(z) - 1) / z with smooth Taylor continuation across z = 0.
template<int N> inline Jet<N> expm1_over_z(const Jet<N>& z) {
  if (std::abs(z.v) < 0.1) {
    // At |z| <= 0.1, 17-term Horner truncation error in d2 is < 7e-27.
    Jet<N> r(1.0);
    for (int k = 16; k >= 1; --k)
      r = Jet<N>(1.0) + (1.0 / static_cast<double>(k + 1)) * (z * r);
    return r;
  }
  return jexpm1(z) / z;
}
// log(1-exp(x)), x < 0, without losing the small complement or overflowing 1/s^2.
template<int N> inline Jet<N> log1mexp(const Jet<N>& x) {
  double s = -std::expm1(x.v), e = std::exp(x.v);
  Jet<N> c(x.v < -0.6931471805599453 ? std::log1p(-e) : std::log(s));
  for (int i = 0; i < N; ++i) c.g[i] = -e * (x.g[i] / s);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      c.h[i][j] = -e * (x.h[i][j] / s) + c.g[i] * (x.g[j] - c.g[j]);
  return c;
}
template<int N> inline Jet<N> logadd(const Jet<N>& a, const Jet<N>& b) {
  if (ISNAN(a.v) || ISNAN(b.v)) return Jet<N>(R_NaN);
  if (a.v == R_NegInf) return b;
  if (b.v == R_NegInf) return a;
  return a.v > b.v ? a + jlog1p(jexp(b-a)) : b + jlog1p(jexp(a-b));
}
template<int N> inline Jet<N> logdiff(const Jet<N>& a, const Jet<N>& b) {
  if (b.v == R_NegInf) return a;
  if (!(b.v < a.v)) return Jet<N>(-1e20);
  return a + log1mexp(b-a);
}

} // namespace egpd_ad

#endif
