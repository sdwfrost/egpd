# Density (PMF) of Compound Poisson-EGPD

Computes the probability mass function of the Compound Poisson-EGPD
distribution using Panjer recursion on a discretized EGPD severity.

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

  number of bins (default: ceiling(max(x)/h) + 50)

- log:

  logical; if TRUE, return log-probabilities

## Value

numeric vector of (log-)probabilities
