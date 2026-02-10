# Predictions from a fitted `egpd` object

Predictions from a fitted `egpd` object

## Usage

``` r
# S3 method for class 'egpd'
predict(
  object,
  newdata,
  type = "link",
  prob = NULL,
  se.fit = FALSE,
  marginal = TRUE,
  trace = 0,
  ...
)
```

## Arguments

- object:

  a fitted `egpd` object

- newdata:

  a data frame

- type:

  a character string: "link", "response", "lpmatrix", or "quantile"

- prob:

  a scalar or vector of probabilities for quantile estimation

- se.fit:

  logical: should standard errors be returned? Defaults to FALSE

- marginal:

  logical: should uncertainty integrate out smoothing parameter
  uncertainty? Defaults to TRUE

- trace:

  an integer controlling output verbosity

- ...:

  unused

## Value

A data frame, list, or design matrix depending on `type`
