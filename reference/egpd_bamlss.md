# bamlss family for continuous EGPD

Creates a `family.bamlss` object for fitting continuous Extended
Generalized Pareto Distribution models with `bamlss()`.

## Usage

``` r
egpd_bamlss(m = 1, ...)
```

## Arguments

- m:

  integer 1–4 selecting the G transformation:

  1

  :   Power: \\G(u) = u^\kappa\\

  2

  :   Mixture: \\G(u) = p u^{\kappa} + (1-p) u^{\delta}\\

  3

  :   Incomplete beta

  4

  :   Power-beta

- ...:

  arguments passed to link specification

## Value

An object of class `family.bamlss`
