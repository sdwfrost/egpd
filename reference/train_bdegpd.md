# \\Experimental\\ Train a neural estimator for bivariate discrete EGPD

Trains Neural Posterior Estimator (NPE) and/or Neural Bayesian Estimator
(NBE) networks for bivariate discrete EGPD inference (BDEGPD or
BZIDEGPD). Requires Julia (\>= 1.11) with NeuralEstimators.jl and
Flux.jl installed, plus the R packages JuliaConnectoR and
NeuralEstimators.

## Usage

``` r
train_bdegpd(
  savepath = "inst/models",
  family = c("bdegpd", "bzidegpd"),
  estimator = c("both", "npe", "nbe"),
  K = NULL,
  m = NULL,
  epochs = NULL,
  stopping_epochs = 10,
  quick = FALSE,
  mc.cores = parallel::detectCores() - 1L,
  seed = 1L,
  verbose = TRUE
)
```

## Arguments

- savepath:

  character; directory to save trained model files (.bson). Created if
  it does not exist.

- family:

  character; `"bdegpd"` (6 parameters) or `"bzidegpd"` (7 parameters,
  adds pi0).

- estimator:

  character; which estimator(s) to train: `"both"` (default), `"npe"`,
  or `"nbe"`.

- K:

  integer; number of training parameter sets to sample from the prior.

- m:

  integer vector; range of sample sizes for simulated datasets.

- epochs:

  integer; maximum number of training epochs.

- stopping_epochs:

  integer; early stopping patience.

- quick:

  logical; if `TRUE`, uses reduced training settings.

- mc.cores:

  integer; number of parallel cores for data simulation.

- seed:

  integer or `NULL`; random seed for reproducibility.

- verbose:

  logical; if `TRUE`, print progress messages.

## Value

A named list of file paths to saved .bson model files (invisible).

## Details

This function is **experimental** and its interface may change.

The training procedure mirrors
[`train_begpd`](https://sdwfrost.github.io/egpd/reference/train_begpd.md)
but uses
[`rbdegpd`](https://sdwfrost.github.io/egpd/reference/rbdegpd.md) or
[`rbzidegpd`](https://sdwfrost.github.io/egpd/reference/rbzidegpd.md)
for data simulation. For BZIDEGPD, the seventh parameter (pi0) uses a
logit transform.

## See also

[`train_begpd`](https://sdwfrost.github.io/egpd/reference/train_begpd.md),
[`rbdegpd`](https://sdwfrost.github.io/egpd/reference/rbdegpd.md),
[`rbzidegpd`](https://sdwfrost.github.io/egpd/reference/rbzidegpd.md)

## Examples

``` r
if (FALSE) { # \dontrun{
paths <- train_bdegpd(savepath = tempdir(), family = "bdegpd", quick = TRUE)
} # }
```
