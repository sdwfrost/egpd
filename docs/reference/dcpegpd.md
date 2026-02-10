# Density (PMF) of Compound Poisson-EGPD

Computes the probability mass function of the Compound Poisson-EGPD
distribution using Panjer recursion on a discretized EGPD severity.

The model is \\S = X_1 + \cdots + X_N\\ where \\N \sim
\mathrm{Poisson}(\lambda)\\ and \\X_i \sim \mathrm{EGPD}(\sigma, \xi,
\kappa, \ldots)\\.

## Usage

``` r
dcpegpd(
  x,
  lambda,
  prob = NA,
  kappa = NA,
  delta = NA,
  sigma,
  xi,
  type = 1,
  h = 0.2,
  K = NULL,
  log = FALSE
)
```

## Arguments

- x:

  numeric vector of quantiles

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

- h:

  bin width for discretization (default 0.2)

- K:

  number of bins (default: automatically determined from data range)

- log:

  logical; if TRUE, return log-probabilities

## Value

numeric vector of (log-)probabilities

## References

Ailliot, P., Gaetan, C. and Naveau, P. (2025). A parsimonious tail
compliant multiscale statistical model for aggregated rainfall. *arXiv
preprint* arXiv:2601.08350. <https://arxiv.org/abs/2601.08350>

## See also

[`pcpegpd`](https://sdwfrost.github.io/egpd/reference/pcpegpd.md),
[`qcpegpd`](https://sdwfrost.github.io/egpd/reference/qcpegpd.md),
[`rcpegpd`](https://sdwfrost.github.io/egpd/reference/rcpegpd.md),
[`fitegpd`](https://sdwfrost.github.io/egpd/reference/fitegpd.md)
