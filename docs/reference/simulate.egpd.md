# Simulations from a fitted `egpd` object

Simulations from a fitted `egpd` object

## Usage

``` r
# S3 method for class 'egpd'
simulate(
  object,
  nsim = 1000,
  seed = NULL,
  newdata,
  type = "link",
  probs = NULL,
  threshold = 0,
  marginal = TRUE,
  ...
)
```

## Arguments

- object:

  a fitted `egpd` object

- nsim:

  an integer giving the number of simulations

- seed:

  an integer giving the seed for simulations

- newdata:

  a data frame

- type:

  a character string: "link" or "response"

- probs:

  a scalar or vector of probabilities

- threshold:

  a scalar added to simulations

- marginal:

  logical: should uncertainty integrate out smoothing parameter
  uncertainty?

- ...:

  additional arguments

## Value

Simulations of parameters
