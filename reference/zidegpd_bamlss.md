# bamlss family for zero-inflated discrete EGPD (ZIDEGPD)

Creates a `family.bamlss` object for fitting zero-inflated discrete
Extended Generalized Pareto Distribution models with `bamlss()`.

## Usage

``` r
zidegpd_bamlss(m = 1, ...)
```

## Arguments

- m:

  integer 1–4 selecting the G transformation (see
  [`egpd_bamlss`](https://sdwfrost.github.io/egpd/reference/egpd_bamlss.md)
  for details)

- ...:

  arguments passed to link specification

## Value

An object of class `family.bamlss`
