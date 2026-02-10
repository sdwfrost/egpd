# Density (PMF) of Compound Poisson-Discrete EGPD

Computes the probability mass function of the Compound Poisson-Discrete
EGPD distribution using Panjer recursion on the exact discrete EGPD
severity PMF (no discretization step).

The model is \\S = X_1 + \cdots + X_N\\ where \\N \sim
\mathrm{Poisson}(\lambda)\\ and \\X_i \sim
\mathrm{Discrete\mbox{-}EGPD}(\sigma, \xi, \kappa, \ldots)\\. Since the
severity is already integer-valued, the Panjer recursion produces an
exact compound distribution with no discretization error.

## Usage

``` r
dcpdegpd(
  x,
  lambda,
  prob = NA,
  kappa = NA,
  delta = NA,
  sigma,
  xi,
  type = 1,
  K = NULL,
  log = FALSE
)
```

## Arguments

- x:

  numeric vector of non-negative integer quantiles

- lambda:

  Poisson rate parameter (\> 0)

- prob:

  G-transformation parameter (type 6 only)

- kappa:

  G-transformation parameter

- delta:

  G-transformation parameter

- sigma:

  GPD scale parameter (\> 0)

- xi:

  GPD shape parameter

- type:

  integer 1-6 specifying G-transformation type

- K:

  grid size for Panjer recursion (default: automatically determined from
  data range)

- log:

  logical; if TRUE, return log-probabilities

## Value

numeric vector of (log-)probabilities

## See also

[`pcpdegpd`](https://sdwfrost.github.io/egpd/reference/pcpdegpd.md),
[`qcpdegpd`](https://sdwfrost.github.io/egpd/reference/qcpdegpd.md),
[`rcpdegpd`](https://sdwfrost.github.io/egpd/reference/rcpdegpd.md),
[`fitegpd`](https://sdwfrost.github.io/egpd/reference/fitegpd.md),
[`dcpegpd`](https://sdwfrost.github.io/egpd/reference/dcpegpd.md)
