# Random generation from Compound Poisson-Discrete EGPD

Direct simulation: draw \\N \sim \mathrm{Poisson}(\lambda)\\, then sum
\\N\\ i.i.d. Discrete-EGPD draws.

## Usage

``` r
rcpdegpd(n, lambda, prob = NA, kappa = NA, delta = NA, sigma, xi, type = 1)
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

integer vector of length `n`

## See also

[`dcpdegpd`](https://sdwfrost.github.io/egpd/reference/dcpdegpd.md),
[`pcpdegpd`](https://sdwfrost.github.io/egpd/reference/pcpdegpd.md),
[`qcpdegpd`](https://sdwfrost.github.io/egpd/reference/qcpdegpd.md)
