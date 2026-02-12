# Random generation from Compound Poisson-EGPD

Direct simulation: draw N ~ Poisson(lambda), then sum N iid EGPD draws.

## Usage

``` r
rcpegpd(n, lambda, prob = NA, kappa = NA, delta = NA, sigma, xi, type = 1)
```

## Arguments

- n:

  number of observations

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

## Value

numeric vector of length n
