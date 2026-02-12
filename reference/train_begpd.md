# Train a neural estimator for bivariate BEGPD

Trains Neural Posterior Estimator (NPE) and/or Neural Bayesian Estimator
(NBE) networks for bivariate BEGPD inference. Requires Julia (\>= 1.11)
with NeuralEstimators.jl and Flux.jl installed, plus the R packages
JuliaConnectoR and NeuralEstimators.

## Usage

``` r
train_begpd(
  savepath = "inst/models",
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

- estimator:

  character; which estimator(s) to train: `"both"` (default), `"npe"`,
  or `"nbe"`.

- K:

  integer; number of training parameter sets to sample from the prior.

- m:

  integer vector; range of sample sizes for simulated datasets (e.g.,
  `1000:4000`).

- epochs:

  integer; maximum number of training epochs.

- stopping_epochs:

  integer; early stopping patience (epochs without improvement).

- quick:

  logical; if `TRUE`, uses reduced training settings for quick testing
  (K=20000, m=500:1000, epochs=10).

- mc.cores:

  integer; number of parallel cores for data simulation. Forced to 1 on
  Windows.

- seed:

  integer or `NULL`; random seed for reproducibility.

- verbose:

  logical; if `TRUE`, print progress messages.

## Value

A named list of file paths to saved .bson model files (invisible).

## Details

The training procedure:

1.  Samples parameters from uniform priors.

2.  Simulates bivariate BEGPD data using
    [`rbegpd`](https://sdwfrost.github.io/egpd/reference/rbegpd.md).

3.  Applies variance-stabilizing (signed log) and log parameter
    transforms.

4.  Trains the neural network(s) using
    [`NeuralEstimators::train()`](https://rdrr.io/pkg/NeuralEstimators/man/train.html).

5.  Saves trained states as .bson files.

## References

Sainsbury-Dale, M., Zammit-Mangion, A., and Huser, R. (2024).
Likelihood-free parameter estimation with neural Bayes estimators. *The
American Statistician*, 78(1), 1–14.

## Examples

``` r
if (FALSE) { # \dontrun{
paths <- train_begpd(savepath = tempdir(), quick = TRUE)
} # }
```
