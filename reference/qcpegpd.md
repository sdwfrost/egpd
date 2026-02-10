# Quantile function of Compound Poisson-EGPD

Computes the quantile function of the Compound Poisson-EGPD distribution
by inverting the Panjer-based CDF.

## Usage

``` r
qcpegpd(
  p,
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

- p:

  numeric vector of probabilities

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

  number of bins (default: 5000)

## Value

numeric vector of quantiles

## See also

[`dcpegpd`](https://sdwfrost.github.io/egpd/reference/dcpegpd.md),
[`pcpegpd`](https://sdwfrost.github.io/egpd/reference/pcpegpd.md),
[`rcpegpd`](https://sdwfrost.github.io/egpd/reference/rcpegpd.md)
