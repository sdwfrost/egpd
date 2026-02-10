# Distribution function of Compound Poisson-EGPD

Computes the cumulative distribution function of the Compound
Poisson-EGPD distribution using Panjer recursion on a discretized EGPD
severity.

## Usage

``` r
pcpegpd(
  q,
  lambda,
  prob = NA,
  kappa = NA,
  delta = NA,
  sigma,
  xi,
  type = 1,
  h = 0.2,
  K = NULL
)
```

## Arguments

- q:

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

## Value

numeric vector of cumulative probabilities

## See also

[`dcpegpd`](https://sdwfrost.github.io/egpd/reference/dcpegpd.md),
[`qcpegpd`](https://sdwfrost.github.io/egpd/reference/qcpegpd.md),
[`rcpegpd`](https://sdwfrost.github.io/egpd/reference/rcpegpd.md)
