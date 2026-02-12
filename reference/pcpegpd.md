# Distribution function of Compound Poisson-EGPD

Distribution function of Compound Poisson-EGPD

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

  number of bins (default: ceiling(max(x)/h) + 50)

## Value

numeric vector of cumulative probabilities
