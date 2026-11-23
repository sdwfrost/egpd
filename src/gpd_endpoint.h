#ifndef EGPD_GPD_ENDPOINT_H
#define EGPD_GPD_ENDPOINT_H
#include <cmath>

// Evaluate (sigma + xi*x)/sigma in scale-normalized units. Near an endpoint
// the sign of the fused residual, not the rounded quotient -sigma/xi,
// determines whether a cell has positive probability.
inline double egpd_endpoint_gap(double x, double sigma, double xi) {
  int ex, es, ea;
  const double mx = std::frexp(x, &ex), ms = std::frexp(sigma, &es);
  const double ma = std::frexp(xi, &ea);
  return std::fma(ma, std::ldexp(mx, ex + ea - es), ms) / ms;
}
#endif
